import sys, os
sys.path.append("..") #this is to add the avobe folder to the package directory
import nexus_tool
from nexus_tool.weap_tools import create_learning_curve
import pandas as pd
import numpy as np

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
w_rate = str(snakemake.params.w_rate)
pv_rate = str(snakemake.params.pv_rate)
grid_rate = str(snakemake.params.grid_rate)
input_path = str(snakemake.input.results)
cropland_path = str(snakemake.input.cropland)
output_path = str(snakemake.output)

# output_main_folder = os.path.join('..', 'Morocco dashboard', 'data', scenario, climate)
# output_lcoe_folder = os.path.join(output_main_folder, f'W{w_rate}_PV{pv_rate}_Grid{grid_rate}')
output_lcoe_folder = output_path.split(os.path.basename(output_path))[0]
os.makedirs(output_lcoe_folder, exist_ok = True)

# file_path = os.path.join(output_main_folder, 'results.gz')
df = nexus_tool.read_csv(input_path)

#Define the path to read the cropland and builtup are data and reads it in
# folder_path = os.path.join('Data', 'Cropland and Builtarea')
# cropland_path = os.path.join(folder_path, 'cropland.gz')
cropland = nexus_tool.read_csv(cropland_path)
cropland = cropland.loc[cropland.Date.isin(df.Date.unique())]

summary_provinces_agri = df.loc[df['type'].str.contains('Agriculture')].groupby(['Province', 'Date'])[['sswd', 'pwd', 'swpa_e', 'swpp_e']].sum()
temp_cropland_provinces = cropland[['province', 'Date']].copy()
temp_cropland_provinces.loc[temp_cropland_provinces['province']=='Inezgane-Aït Melloul', 'province'] = 'Taroudannt'

for feature in list(summary_provinces_agri):
    cropland[feature] = temp_cropland_provinces.set_index(['province', 'Date']).index.map(summary_provinces_agri[feature]) * cropland['area_share']

sm_cropland = nexus_tool.Model(cropland)

sm_cropland.start_year = 2020
sm_cropland.end_year = 2050
sm_cropland.discount_rate = 0.05
sm_cropland.pwd = 'pwd'
sm_cropland.sswd = 'sswd'
sm_cropland.pd_e = 'swpp_e'
sm_cropland.swpp_e = 'swpp_e'
sm_cropland.swpa_e = 'swpa_e'

cc = 1250
rates = {2020: 0, 2030: float(w_rate)*0.3, 2040: float(w_rate), 2050: float(w_rate)*1.3}
cc_wind_curve = create_learning_curve(sm_cropland.start_year, sm_cropland.end_year, cc, rates, 'polynomial')

sm_cropland.create_wind_turbine('Wind power', life=20,
                                om_cost=0.01, capital_cost=cc_wind_curve,
                                efficiency=1)
sm_cropland.technologies['Wind power'].p_rated = 80
sm_cropland.technologies['Wind power'].p_curve = [0, 0, 0, 0, 2.9, 6, 11, 17.7, 27.3, 39.2, 51.4, 63.8, 74.2, 79.9, 82.2, 82.9, 83.3, 83.3, 83, 83, 83]
sm_cropland.technologies['Wind power'].u_arr = range(0, 21)
sm_cropland.technologies['Wind power'].z = 39

cc = 1000
rates = {2020: 0, 2030: float(pv_rate)*0.3, 2040: float(pv_rate), 2050: float(pv_rate)*1.3}
cc_solar_curve = create_learning_curve(sm_cropland.start_year, sm_cropland.end_year, cc, rates, 'polynomial')

sm_cropland.create_pv_system('Solar PV', life=15,
                             om_cost=0.01, capital_cost=cc_solar_curve,
                             efficiency=1)

cc = 0.12
rates = {2020: 0, 2030: float(grid_rate)*0.3, 2040: float(grid_rate), 2050: float(grid_rate)*1.3}
e_price_curve = create_learning_curve(sm_cropland.start_year, sm_cropland.end_year, cc, rates, 'step')

sm_cropland.create_standard_tech('Grid pump', life=15, om_cost=0.1,
                                 capital_cost=0, fuel_cost=e_price_curve,
                                 fuel_req=1, efficiency=0.85, cf = 0.8,
                                 emission_factor=1.76, env_cost=0)

sm_cropland.get_cf('all', axis=0)
sm_cropland.get_installed_capacity('all', axis=0)
sm_cropland.get_max_capacity('a', axis=0)
sm_cropland.get_lcoe(years='all', axis=0)
sm_cropland.get_least_cost(technologies='a', years='all', axis=0)

sm_cropland.lcoe.drop(columns=['water demand', 'required capacity', 'energy demand'], inplace=True)
sm_cropland.lcoe.reset_index().to_csv(output_path, index=False)