#This is a test file for butane calculations based on energy_model:

import sys, os
sys.path.append("..") #this is to add the avobe folder to the package directory
import nexus_tool
from nexus_tool.weap_tools import create_learning_curve
import pandas as pd
import numpy as np

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
#demand_data = str(snakemake.input.demand_data)
#wwtp_inflow = str(snakemake.input.wwtp_inflow)
output = str(snakemake.output)

# load_folder = os.path.join('Data', 'Processed results', scenario, climate)
# results_folder = os.path.join('..', 'Morocco dashboard', 'data', scenario, climate)
# results_folder = output.split('results.gz')[0]

results_folder = output.split(os.path.basename(output))[0]

os.makedirs(results_folder, exist_ok = True)


#Define the path to read the cropland and builtup are data and reads it in
folder_path = os.path.join('data')
#cropland_path = os.path.join(folder_path,'Cropland and Builtarea' ,'cropland.gz')
cropland_path = str(snakemake.input.cropland)
cropland = nexus_tool.read_csv(cropland_path)
weap_path = str(snakemake.input.results)
weap = pd.read_csv(weap_path)
weap.loc[weap['Province']=='Chtouka-Ait Baha', 'Province'] = 'Chtouka-Aït Baha'
weap.loc[weap['Province']=='Inezgane-Ait Melloul', 'Province'] = 'Inezgane-Aït Melloul'
weap.loc[weap['Province']=='Agadir-Ida-Ou-Tanane', 'Province'] = 'Agadir-Ida ou Tanane'

#weap: the results of the enegy model based on weap demand nodes
cropland = cropland.loc[cropland.Date.isin(weap.Date.unique())] #this should be the results dataframe

summary_provinces_agri = weap.loc[weap['type'].str.contains('Agriculture')].groupby(['Province', 'Date'])[['sswd', 'pwd', 'swpa_e', 'swpp_e']].sum()
temp_cropland_provinces = cropland[['province', 'Date']].copy()
temp_cropland_provinces.loc[temp_cropland_provinces['province']=='Inezgane-Aït Melloul', 'province'] = 'Taroudannt'
for feature in list(summary_provinces_agri):
    cropland[feature] = temp_cropland_provinces.set_index(['province', 'Date']).index.map(summary_provinces_agri[feature]) * cropland['area_share']

cropland.rename(columns={'sswd':'water_demand(m3)',
                'pwd':'peak_water(m3)',
                'swpp_e':'peak_load(KW)',
                'swpa_e':'energy_demand(KWh)'}, inplace=True)


pv_year = 2040

butane_phaseout = str(snakemake.params.butane_phaseout)
pv_level = int(snakemake.params.pv_level)

discount_rate = 0.05

# Step1: Reading the input data file
#         file_path = os.path.join(load_folder, 'cropland.csv')
#         df = pd.read_csv(file_path)
df = cropland

# Step2: Modeling Butane Phase out rate
def create_tech_dist(iyear, eyear, share1, rates, method, order=2):
    tech_dist = pd.DataFrame({'Year': range(iyear, eyear+1)})
    if method!='step':
        for year, rate in rates.items():
            tech_dist.loc[tech_dist.Year==year, 'rate'] = rate
        tech_dist['share'] = (1-tech_dist['rate']) * share1
        tech_dist['share'] = tech_dist.share.interpolate(method=method, order=order)
    else:
        for year, rate in rates.items():
            tech_dist.loc[tech_dist.Year>=year, 'rate'] = rate
        tech_dist['share'] = (1-tech_dist['rate']) * share1
    return tech_dist.set_index('Year').share

                
#         butane_change ={'bau':0,'late_po':1,'early_po':1 }
if butane_phaseout!='None':
    butane_share = 1
    butane_phaseout = int(butane_phaseout)
    year = butane_phaseout
else:
    butane_share = 0
    year = 2050
        

# #         end_year ={'bau':2050,'late_po':2040,'early_po':2030 }
#         year = end_year[scenario]

#         pv_year ={'bau':2040,'late_po':2040,'early_po':2040}

#         pv_change ={'low_PV':0, 'mid_PV':-1, 'high_PV':-4} #The negative sign means increase in share. 
#         pv_share = pv_change[pv_level]
pv_share = 1 - pv_level / 10
        
butane_share1 = 0.2
pv_share1 = 0.05
df.start_year = 2020
df.end_year = year
#         pv_year = pv_year[scenario]
df.mid_year = ((df.start_year+df.end_year)/2)
bt_rates = {df.start_year: 0, df.mid_year: butane_share*0.5, df.end_year: butane_share}
pv_rates = {df.start_year: 0, 2030: pv_share*0.5, pv_year: pv_share, 2050: pv_share*1.2}
butane_dist = create_tech_dist(df.start_year, df.end_year, butane_share1, bt_rates, 'polynomial')
pv_dist = create_tech_dist(df.start_year, pv_year, pv_share1, pv_rates, 'polynomial')

#Step3: PV learning curve and crop in PV CAPEX
def create_learning_curve(iyear, eyear, cc, rates, method, order=2):
    learning_curve = pd.DataFrame({'Year': range(iyear, eyear+1)})
    if method!='step':
        for year, rate in rates.items():
            learning_curve.loc[learning_curve.Year==year, 'rate'] = rate
        learning_curve['capital_cost'] = (1-learning_curve['rate']) * cc
        learning_curve['capital_cost'] = learning_curve.capital_cost.interpolate(method=method, order=order)
    else:
        for year, rate in rates.items():
            learning_curve.loc[learning_curve.Year>=year, 'rate'] = rate
        learning_curve['capital_cost'] = (1-learning_curve['rate']) * cc
    return learning_curve.set_index('Year').capital_cost

pv_rate = 0.4 #Assuming 40% drop in PV capex by 2040.
cc = 7 #current capital cost of pv in MAD/wat 
#pv_cf = 0.25
pv_life=25
df.start_year=2020
df.end_year= 2050
#capital_cost = 1000
rates = {2020: 0, 2030: pv_rate*0.3, 2040: pv_rate, 2050: pv_rate*1.3}
cc_solar_curve = create_learning_curve(df.start_year, df.end_year, cc, rates, 'polynomial')
om_cost=0.01 
efficiency=1

dist = pd.DataFrame({'butane_share':pd.Series(butane_dist),
            'pv_share':pd.Series(pv_dist),
            'pv_cc':pd.Series(cc_solar_curve)})

dist['butane_share'].fillna(method='ffill',inplace=True) #This method fills any missing data with the previous value.
dist['pv_share'].fillna(method='ffill',inplace=True)
dist['pv_cc'].fillna(method='ffill',inplace=True)
dist['grid_share']= 1-(dist['pv_share']+dist['butane_share'])


souss_massa = pd.merge(df, dist, on='Year')

butane_req = 1/12.58 # LHV = 12.58 KWh/Kg (amount of butane (kg) required to produce 1 KWh) 
butane_em = 6.67 #kg CO"/Gallon
gallon2liter = 3.78541 #l/Gallon
butane_density = 573 #kg/m3
butane_ef = butane_em / (gallon2liter/1000 * butane_density) #kgCO2/kg
butane_em_fac = butane_ef * 1000 ##kgCO2/ton
#         butane_em_fac = 3079.5
bpump_eff = 0.2 #efficiency of butane 
epump_eff = 0.45 #assumption
conv_fac = 1000000000 # to convert emissions from kgCO2 to Million meteric tons of CO2 MtCO2

souss_massa['pv_load(KW)'] = souss_massa['peak_load(KW)']*souss_massa['pv_share']
souss_massa['pv_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['pv_share']
souss_massa['butane_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['butane_share']
souss_massa['grid_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['grid_share']


souss_massa['pv_elec(KWh)'] = souss_massa['pv_demand(KWh)']/epump_eff
souss_massa['grid_elec(KWh)'] = souss_massa['grid_demand(KWh)']/epump_eff

souss_massa['butane_cons(kg)'] = (souss_massa['butane_demand(KWh)']*butane_req)/bpump_eff
souss_massa['butane_cons(tonnes)'] = souss_massa['butane_cons(kg)']/1000
souss_massa['butane_emissions(MtCO2)'] = souss_massa['butane_cons(tonnes)'] * butane_em_fac/conv_fac

souss_massa['butane_FARcost(mMAD)'] = (souss_massa['butane_cons(kg)']*(40/12))/1000000 #in million MAD, this is what farmers pay
souss_massa['butane_ACTcost(mMAD)'] = (souss_massa['butane_cons(kg)']*(120/12))/1000000 #in million MAD, this is what farmers pay
souss_massa['butane_Subsidy(mMAD)'] = (souss_massa['butane_cons(kg)']*(80/12))/1000000 #in million MAD, this is the total subsidy cost


# Estimating the electriicty from the grid and emissions
grid_em_fac = 0.7 #kgco2/kwh  , This is based on data from IEA, Morocco 2019 
grid_cost = 0.57 # MAD/KWh, Assumption and to be updated 

#scenario1['grid_elec (kwh)'] = scenario1['grid_demand_kwh']/epump_eff #already calculated above
souss_massa['grid_emissions(MtCO2)'] = souss_massa['grid_elec(KWh)']*grid_em_fac/conv_fac
souss_massa['grid_cost(mMAD)'] = (souss_massa['grid_elec(KWh)']*grid_cost)/1000000
souss_massa['Total_emissions(MtCO2)'] = (souss_massa['grid_emissions(MtCO2)'] + souss_massa['butane_emissions(MtCO2)'])

# Estimating the required monthly capacity of pv NOTE: 1 kWh = 3600 kJ

souss_massa['cf'] = souss_massa['srad'] / (24*60*60) #This will give the cf in solar rad: (kJ/(m2.day))*30.day/month*1h/(24h*60m*60s) =kWh/(m2.month)*30/(60*60)
souss_massa['cap_m(MW)'] = souss_massa['pv_load(KW)'] / souss_massa['cf']/1000   #to convert to MW, check the units




souss_massa1 = souss_massa.groupby(['Demand point','Year']).agg({'water_demand(m3)': 'sum','energy_demand(KWh)': 'sum', 
                                                                    'pv_elec(KWh)': 'sum', 'grid_elec(KWh)': 'sum', 'cap_m(MW)': 'max',
                                                                    'butane_cons(tonnes)': 'sum', 'butane_FARcost(mMAD)': 'sum', 
                                                                    'butane_ACTcost(mMAD)': 'sum','butane_Subsidy(mMAD)': 'sum',
                                                                    'butane_emissions(MtCO2)': 'sum','grid_emissions(MtCO2)': 'sum',
                                                                    'Total_emissions(MtCO2)': 'sum','grid_cost(mMAD)': 'sum',
                                                                    'pv_demand(KWh)': 'sum', 'butane_demand(KWh)': 'sum', 'grid_demand(KWh)': 'sum'})

pv_installed_cap = pd.Series(dtype=float) #inicialize a pandas series that will be populated with the cumulative max of the max capacity of each point group
for index, group in souss_massa1.reset_index().groupby('Demand point'): # loops through each demand point set of data
    group_pv_cap = pd.Series(group['cap_m(MW)'].cummax().values, index=group.reset_index().set_index(['Demand point','Year']).index)
    pv_installed_cap = pv_installed_cap.append(group_pv_cap) #, ignore_index=True) # calculates the cummmax() for the demand point and append the values to the pv_installed_capacity

#souss_massa1['Province'] = (souss_massa.groupby(['Demand point','Year'])['province'])
souss_massa1['GWdepth'] = (souss_massa.groupby(['Demand point','Year'])['wtd'].mean())
souss_massa1['srad'] = (souss_massa.groupby(['Demand point','Year'])['srad'].mean())
souss_massa1['wind'] = (souss_massa.groupby(['Demand point','Year'])['wind'].mean())
souss_massa1['cap_mean(MW)'] = (souss_massa.groupby(['Demand point','Year'])['cap_m(MW)'].mean())
souss_massa1['cap_max(MW)'] = (souss_massa.groupby(['Demand point','Year'])['cap_m(MW)'].max())
souss_massa1['pv_cc'] = (souss_massa.groupby(['Demand point','Year'])['pv_cc'].mean())



#Calculating PV installed capacity and new capacity in each year:
#souss_massa1.reset_index(inplace=True)     
souss_massa1['PV_installed_cap(MW)'] = pv_installed_cap
souss_massa1['PV_new_cap(MW)'] = souss_massa1['PV_installed_cap(MW)'] - souss_massa1['PV_installed_cap(MW)'].shift(1)
souss_massa1.reset_index(inplace=True)
souss_massa1.loc[souss_massa1['Year']==2020, 'PV_new_cap(MW)'] = 0


#Calculating the required area for PV installations:
souss_massa1['PV_area(ha)'] = souss_massa1['PV_installed_cap(MW)'] * 0.9  # Since the area required to install 1 MW of PV = 1 ha or 10m2/1KW 



#PV reinvestment calculations
#         souss_massa2 = pd.DataFrame()
#         for index, group in souss_massa1.groupby(['Demand point']):
#             dff = group.copy()
#             dff['reinv_cap(MW)'] = dff['PV_new_cap(MW)'].shift(pv_life).fillna(0)
#             souss_massa2 = souss_massa2.append(dff, ignore_index=True)
    
#             souss_massa2['PV_Capex(mMAD)']=(souss_massa2['PV_new_cap(MW)']+souss_massa2['reinv_cap(MW)'])*souss_massa2['pv_cc']
#             souss_massa2['PV_Opex(mMAD)']=(souss_massa2['PV_Capex(mMAD)']*om_cost)

souss_massa1['reinv_cap(MW)'] = souss_massa1['PV_new_cap(MW)'].shift(pv_life).fillna(0)
souss_massa1.loc[souss_massa1['Year']<(2020+pv_life), 'reinv_cap(MW)'] = 0
souss_massa1['PV_Capex(mMAD)']=(souss_massa1['PV_new_cap(MW)']+souss_massa1['reinv_cap(MW)'])*souss_massa1['pv_cc']
souss_massa1['PV_Opex(mMAD)']=(souss_massa1['PV_Capex(mMAD)']*om_cost)

#NPV calculations:
souss_massa1['time'] = souss_massa1['Year']-2020

souss_massa1['PV_Capex_NPV(mMAD)'] = souss_massa1['PV_Capex(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['PV_Opex_NPV(mMAD)'] = souss_massa1['PV_Opex(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['butane_Subsidy_NPV(mMAD)'] = souss_massa1['butane_Subsidy(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['grid_cost_NPV(mMAD)'] = souss_massa1['grid_cost(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['PV_Total_NPV(mMAD)'] = souss_massa1['PV_Capex_NPV(mMAD)'] + souss_massa1['PV_Opex_NPV(mMAD)']

souss_massa1_summary = souss_massa1.groupby(['Year'])[['water_demand(m3)','energy_demand(KWh)', 'cap_max(MW)', 'pv_elec(KWh)', 
                                                        'grid_elec(KWh)','butane_cons(tonnes)', 'butane_FARcost(mMAD)',
                                                        'PV_new_cap(MW)','PV_installed_cap(MW)','PV_area(ha)','reinv_cap(MW)',
                                                        'butane_emissions(MtCO2)','grid_emissions(MtCO2)','Total_emissions(MtCO2)',
                                                        'butane_Subsidy(mMAD)','butane_Subsidy_NPV(mMAD)','grid_cost(mMAD)','grid_cost_NPV(mMAD)',
                                                        'PV_Capex(mMAD)','PV_Capex_NPV(mMAD)','PV_Opex(mMAD)','PV_Opex_NPV(mMAD)', 'PV_Total_NPV(mMAD)',
                                                        'pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum()

# Step: Saving Results
souss_massa1_summary.to_csv(output)

