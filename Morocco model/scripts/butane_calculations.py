#!/usr/bin/env python
# coding: utf-8

# # NEXUS tool: case study for the Souss-Massa basin - Butane phase out strategies 
# In this notebook a case study for the Souss-Massa basin is covered using the `nexustool` and other packages. This part of the analysis will consider three phase out scenarios for Butane phaseout as follows:
# 
# 1. Business as usual scenario (None): No phase out which means the current share of butane use to continue into the future
# 2. Late Phase out (2040): which assumes compelte phase out of butane by 2040.
# 3. Early Phase out (2030): which assumes compelte phase out of butane by 2030.
# 
# For all the three main scenarios, we will explore three levels of PV share in the region.
# 
# a. Low PV Share (10): which assumes the current level of PV share of 10% to continue in the future
# b. Medium PV Share (20): which assumes an increase of PV share from 10% to 20% by 2040.
# c. High PV Share (50): which assumes an increase of PV share from 10% to 50% by 2040.
# 
# First import the package by running the following block:

# In[ ]:


import sys
sys.path.append("..") #this is to add the above folder to the package directory
import os
import nexustool
from nexustool.weap_tools import create_learning_curve, create_tech_dist
from nexustool.gis_tools import disaggregate_data
import pandas as pd
import numpy as np
from dashboard.scripts.plotting import pv_installed_capacity, energy_demand_ag, emissions_ag, costs_plot


# ## Step 1:Reading the input files
# 
# The analysis in this phase will depend on the outcomes of WEAP and the energy demand calculations. The main two input files are: 
# 1) Cropland file: which shows the spatial distribution of irrigated cropland area.
# 
# 2) WEAP file: which is shows the energy demand results fir each of the WEAP scenarios. In the following steps the naming of some of the provinces is also corrected. 

# In[ ]:


#Define the path to read the cropland and builtup are data and reads it in
maincropland_path = snakemake.input.cropland
maincropland = pd.read_csv(maincropland_path)
weap = pd.read_csv(str(snakemake.input.results))

output = str(snakemake.output)

results_folder = output.split(os.path.basename(output))[0]

os.makedirs(results_folder, exist_ok = True)


# ## Step2: Cleaning and reorganizing the cropland dataframe
# In this step we take the `WEAP` and distribute the energy scenario results into the `cropland` data. This way we can correlate agricultural pumping energy demans with spatially distributed irrigated area and solar irradiation (PV potential).

# In[ ]:


cropland = disaggregate_data(weap, maincropland)


# In[ ]:


cropland.rename(columns={'sswd':'water_demand(m3)',
                'pwd':'peak_water(m3)',
                'pp_e':'peak_load(KW)',
                'pa_e':'energy_demand(KWh)'}, inplace=True)


# ## Step3: Define butane phase out year and PV share
# 
# In this step we introduce a new function that will model the change in the technologies disctributions over years. To have a better understanding of this step, check the scenarios section under Butane calculations in the Energy Training material. 

# In[ ]:


df = cropland

butane_phaseout = str(snakemake.params.butane_phaseout)
pv_level = int(snakemake.params.pv_level)

pv_year = 2040
discount_rate = 0.05


# In[ ]:


if butane_phaseout:
    butane_share = 1
    year = butane_phaseout
else:
    butane_share = 0
    year = 2050


# ## Step4: Defining technologies distributions over coming years
# 
# In this step we introduce a new function that allows us to simulate the drop in PV capital cost due to learning curve. In this analysis we assume that the PV capital cost we drop to 40% of the currrent level in 2020. 

# In[ ]:


pv_share = 1 - pv_level / 10

om_cost = 0.01 
efficiency = 1
butane_share1 = 0.2
pv_share1 = 0.1
df.start_year = 2020
df.end_year = year
df.mid_year = ((df.start_year+df.end_year)/2)
bt_rates = {df.start_year: 0, df.mid_year: butane_share*0.5, df.end_year: butane_share}
pv_rates = {df.start_year: 0, 2030: pv_share*0.5, pv_year: pv_share, 2050: pv_share*1.2}
butane_dist = create_tech_dist(df.start_year, df.end_year, butane_share1, bt_rates, 'polynomial')
pv_dist = create_tech_dist(df.start_year, pv_year, pv_share1, pv_rates, 'polynomial')

pv_rate = 0.4 #Assuming 40% drop in PV capex by 2040.
cc = 7 #current capital cost of pv in MAD/wat. source: energy experts in SM
pv_life = 15
df.start_year = 2020
df.end_year= 2050

rates = {2020: 0, 2030: pv_rate*0.3, 2040: pv_rate, 2050: pv_rate*1.3}
cc_solar_curve = create_learning_curve(df.start_year, df.end_year, cc, rates, 'polynomial')


# Intorducing a new dataframe to show the change in butane share, pv share and pv capital cost (cc)
dist = pd.DataFrame({'butane_share':pd.Series(butane_dist),
                   'pv_share':pd.Series(pv_dist),
                    'pv_cc':pd.Series(cc_solar_curve)})

dist['butane_share'].fillna(method='ffill',inplace=True) #This method fills any missing data with the previous value.
dist['pv_share'].fillna(method='ffill',inplace=True)
dist['pv_cc'].fillna(method='ffill',inplace=True)
dist['grid_share']= 1-(dist['pv_share']+dist['butane_share'])


# In[ ]:


# merging the distribution dataframe with the original dataframe that includes weap outputs, the merging is based on 'Year'

souss_massa = pd.merge(df, dist, on='Year')


# ## Step5: Calculating the monthly electricity demand for each technology 
# 
# Three main technologies used for groundwater pumping: Butane pumps, electric pumps taking electricty from the national grid and solar pumps. Since we do not have detailed information on the exact location where each technology is being used, we assumed that in each area there is a mix of the three technologies and the total pumping demand is split between the three technologies based on the share of each technology (Grid:70%, Butane: 20% and PV:10%).
# In this step we calculate the monthly demand for each technology in each irrigated area

# In[ ]:


# Introducing technologies characteristics: such as costs, emission factors, pumps effiencies ..etc

butane_req = 1/12.58 # LHV = 12.58 KWh/Kg (amount of butane (kg) required to produce 1 KWh) 
butane_em = 6.67 #kg CO"/Gallon
gallon2liter = 3.78541 #l/Gallon
butane_density = 573 #kg/m3
butane_ef = butane_em / (gallon2liter/1000 * butane_density) #kgCO2/kg
butane_em_fac = butane_ef * 1000 ##kgCO2/ton
butane_em_fac = 3079.5
bpump_eff = 0.2 #efficiency of butane 
epump_eff = 0.45 #assumption
conv_fac = 1000000000 # to convert emissions from kgCO2 to Million meteric tons of CO2 MtCO2
grid_em_fac = 0.7 #kgco2/kwh  , This is based on data from IEA, Morocco 2019 
grid_cost = 0.57 # MAD/KWh, Assumption and to be updated 


# In[ ]:


souss_massa['pv_load(KW)'] = souss_massa['peak_load(KW)']*souss_massa['pv_share']
souss_massa['pv_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['pv_share']
souss_massa['butane_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['butane_share']
souss_massa['grid_demand(KWh)'] = souss_massa['energy_demand(KWh)']*souss_massa['grid_share']


souss_massa['pv_elec(KWh)'] = souss_massa['pv_demand(KWh)']/epump_eff
souss_massa['grid_elec(KWh)'] = souss_massa['grid_demand(KWh)']/epump_eff

souss_massa['butane_cons(kg)'] = (souss_massa['butane_demand(KWh)']*butane_req)/bpump_eff
souss_massa['butane_cons(tonnes)'] = souss_massa['butane_cons(kg)']/1000


# ## Step6: Calculating the costs and emissions for butane pumps and the grid
# 
# After calculating the share of each technology, we move to calculate the cost of using each technology. 
# For Butane we are calculating three types of costs: what farmers pay, the subsidy level and the toatl cost. The most important one is the second cost (subsidy) which will be used in the next steps.
# In this step we also compute the emissions released due to the use of butane or the grid. knowing that the national grid is heavily dependant on fossil fuels. 

# In[ ]:


# Estimating butance costs and emissions:

souss_massa['butane_emissions(MtCO2)'] = souss_massa['butane_cons(tonnes)'] * butane_em_fac/conv_fac
souss_massa['butane_FARcost(mMAD)'] = (souss_massa['butane_cons(kg)']*(40/12))/1000000 #in million MAD, this is what farmers pay
souss_massa['butane_ACTcost(mMAD)'] = (souss_massa['butane_cons(kg)']*(120/12))/1000000 #in million MAD, this is what farmers pay
souss_massa['butane_Subsidy(mMAD)'] = (souss_massa['butane_cons(kg)']*(80/12))/1000000 #in million MAD, this is the total subsidy cost


# Estimating grid costs and emissions:

souss_massa['grid_emissions(MtCO2)'] = souss_massa['grid_elec(KWh)']*grid_em_fac/conv_fac
souss_massa['grid_cost(mMAD)'] = (souss_massa['grid_elec(KWh)']*grid_cost)/1000000
souss_massa['Total_emissions(MtCO2)'] = (souss_massa['grid_emissions(MtCO2)'] + souss_massa['butane_emissions(MtCO2)'])


# ## Step7: Calculate the required PV capacity for pumping (for each location in each month) 
# 
# Now let us move to the PV part of the calculations. We first use the solar radiation map to estimate the capacity factor at each location in each month. This will then be used to calculate the required capacity of PV in each month of the year. 

# In[ ]:


# Estimating the required monthly capacity of pv NOTE: 1 kWh = 3600 kJ

souss_massa['cf'] = souss_massa['srad'] / (24*60*60) #This will give the cf in solar rad: (kJ/(m2.day))*30.day/month*1h/(24h*60m*60s) =kWh/(m2.month)*30/(60*60)
souss_massa['cap_m(MW)'] = souss_massa['pv_load(KW)'] / souss_massa['cf']/1000   #to convert to MW, check the units


# ## Step 8: Calculate the annual required PV capacity (new installation, reinvestment and total capacity)
# 
# When we talk about investments we are usually interested in annual investment not monthly. Therefore we aggregate the dataframe to show the annual values for each type of inputs. The `groupby` method is used to group everything for each `demand point` and for each `Year` as shown here:

# In[ ]:


souss_massa1 = souss_massa.groupby(['Demand point','Year']).agg({'water_demand(m3)': 'sum','energy_demand(KWh)': 'sum', 
                                    'pv_elec(KWh)': 'sum', 'grid_elec(KWh)': 'sum', 'cap_m(MW)': 'max',
                                     'butane_cons(tonnes)': 'sum', 'butane_FARcost(mMAD)': 'sum', 
                                     'butane_ACTcost(mMAD)': 'sum','butane_Subsidy(mMAD)': 'sum',
                                      'butane_emissions(MtCO2)': 'sum','grid_emissions(MtCO2)': 'sum',
                                      'Total_emissions(MtCO2)': 'sum','grid_cost(mMAD)': 'sum',
                                      'pv_demand(KWh)': 'sum', 'butane_demand(KWh)': 'sum', 'grid_demand(KWh)': 'sum'})


# In[ ]:


# Introducing additional attributes to the dataframe:

souss_massa1['GWdepth'] = (souss_massa.groupby(['Demand point','Year'])['wtd'].mean())
souss_massa1['srad'] = (souss_massa.groupby(['Demand point','Year'])['srad'].mean())
souss_massa1['cap_mean(MW)'] = (souss_massa.groupby(['Demand point','Year'])['cap_m(MW)'].mean())
souss_massa1['cap_max(MW)'] = (souss_massa.groupby(['Demand point','Year'])['cap_m(MW)'].max())
souss_massa1['pv_cc'] = (souss_massa.groupby(['Demand point','Year'])['pv_cc'].mean())

del souss_massa


# Here we set the annual capacity of PV to be the maximum monthly capacity at each demand point. 

# In[ ]:


pv_installed_cap = pd.Series(dtype=float) #inicialize a pandas series that will be populated with the cumulative max of the max capacity of each point group
for index, group in souss_massa1.reset_index().groupby('Demand point'): # loops through each demand point set of data
    group_pv_cap = pd.Series(group['cap_m(MW)'].cummax().values, index=group.reset_index().set_index(['Demand point','Year']).index)
    pv_installed_cap = pv_installed_cap.append(group_pv_cap) # calculates the cummmax() for the demand point and append the values to the pv_installed_capacity


# In this step we introduce the `new capcity` which is the additional capacity of PV required in each year compare to the previous year. Also we intorduce the `reinvest capacity` which is the second investment required in certain locations after the lifetime of the pannels

# In[ ]:


#souss_massa1.reset_index(inplace=True)     
souss_massa1['PV_installed_cap(MW)'] = pv_installed_cap
souss_massa1['PV_new_cap(MW)'] = souss_massa1['PV_installed_cap(MW)'] - souss_massa1['PV_installed_cap(MW)'].shift(1)
souss_massa1.reset_index(inplace=True)
souss_massa1.loc[souss_massa1['Year']==2020, 'PV_new_cap(MW)'] = 0
souss_massa1['reinv_cap(MW)'] = souss_massa1['PV_new_cap(MW)'].shift(pv_life).fillna(0)
souss_massa1.loc[souss_massa1['Year']<(2020+pv_life), 'reinv_cap(MW)'] = 0


# ## Step 9: Calculate PV CAPEX and OPEX
# 
# After calculating the required capacity of PV, it is time to calculate PV capital cost (CAPEX) which is based in the annual capacity and the annual price of PV pannels given the annual drop due to learning curve. The Operating and Maintenance cost of PV (OPEX) are simply a fraction of the annual CAPEX.

# In[ ]:


# Calculating PV CAPEX and OPEX:

souss_massa1['PV_Capex(mMAD)']=(souss_massa1['PV_new_cap(MW)']+souss_massa1['reinv_cap(MW)'])*souss_massa1['pv_cc']
souss_massa1['PV_Opex(mMAD)']=(souss_massa1['PV_Capex(mMAD)']*om_cost)


# In[ ]:


#Calculating the required area for PV installations:
souss_massa1['PV_area(ha)'] = souss_massa1['PV_installed_cap(MW)'] * 0.9  # Since the area required to install 1 MW of PV = 1 ha or 10m2/1KW 


# In[ ]:


#NPV calculations:
souss_massa1['time'] = souss_massa1['Year']-2020

souss_massa1['PV_Capex_NPV(mMAD)'] = souss_massa1['PV_Capex(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['PV_Opex_NPV(mMAD)'] = souss_massa1['PV_Opex(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['butane_Subsidy_NPV(mMAD)'] = souss_massa1['butane_Subsidy(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['grid_cost_NPV(mMAD)'] = souss_massa1['grid_cost(mMAD)']/((1+discount_rate)**(souss_massa1['time']))
souss_massa1['PV_Total_NPV(mMAD)'] = souss_massa1['PV_Capex_NPV(mMAD)'] + souss_massa1['PV_Opex_NPV(mMAD)']


# ## Step 10: Summary dataframe 
# 
# Here we construct a summary dataframe that includes the final results that we would like to explore. 

# In[ ]:


souss_massa1_summary = souss_massa1.groupby(['Year'])[['water_demand(m3)','energy_demand(KWh)', 'cap_max(MW)', 'pv_elec(KWh)', 
                        'grid_elec(KWh)','butane_cons(tonnes)', 'butane_FARcost(mMAD)',
                        'PV_new_cap(MW)','PV_installed_cap(MW)','PV_area(ha)','reinv_cap(MW)',
                        'butane_emissions(MtCO2)','grid_emissions(MtCO2)','Total_emissions(MtCO2)',
                        'butane_Subsidy(mMAD)','butane_Subsidy_NPV(mMAD)','grid_cost(mMAD)','grid_cost_NPV(mMAD)',
                        'PV_Capex(mMAD)','PV_Capex_NPV(mMAD)','PV_Opex(mMAD)','PV_Opex_NPV(mMAD)', 'PV_Total_NPV(mMAD)',
                        'pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum()


# ## Step 11: Save results files
# 
# The results can then be saved into a defined output folder `butane_results_folder`:

# In[ ]:


# Step: Saving Results

souss_massa1_summary.reset_index(inplace=True)
souss_massa1_summary.to_csv(output)
