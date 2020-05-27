import sys, os
sys.path.append("..") #this is to add the avobe folder to the package directory
import nexus_tool
from nexus_tool.weap_tools import create_folder
from nexus_tool.weap_tools import create_learning_curve
import pandas as pd
import numpy as np

scenario = str(sys.argv[1])
climate = str(sys.argv[2])

load_folder = os.path.join('Data', 'Processed results')
results_folder = os.path.join('..', 'Morocco dashboard', 'data')
# results_folder = os.path.join('Results')
create_folder(results_folder)

scenario_folder = os.path.join(load_folder, scenario)
create_folder(os.path.join(results_folder, scenario))

sub_scenario_folder = os.path.join(scenario_folder, climate)
create_folder(os.path.join(results_folder, scenario, climate))
output_main_folder = os.path.join(results_folder, scenario, climate)

load_data = sub_scenario_folder
#Define the path to read the scenario input data and reads it in
file_path = os.path.join(load_data, 'demand_data.gz')
df = nexus_tool.read_csv(file_path)

#Creates the nexus model with the input dataframe
souss_massa = nexus_tool.Model(df)

#Define the diameter of the pipelines. The first is general for all and the following specific for each case
souss_massa.df['Pipe_diameter'] = 0.4
souss_massa.df.loc[souss_massa.df['type'].str.contains('GW'), 'Pipe_diameter'] = 1
souss_massa.df.loc[souss_massa.df['type'].str.contains('Pipeline'), 'Pipe_diameter'] = 1.2
souss_massa.df.loc[souss_massa.df['type'].str.contains('DS'), 'Pipe_diameter'] = 1

#Define the variable to take into account for elevation difference and lenght of pipelines
souss_massa.elevation = 'elevation_diff' #for the case of GW, the elevation_diff is set to be the wtd
souss_massa.L = 'distance' #for the case of GW, the distance is set to be the wtd

#Defines the name of the variable for Peak Water Demand and Seasonal Water demand (monthly)
souss_massa.pwd = 'pwd'
souss_massa.sswd = 'sswd'
souss_massa.df.rename(columns={'value': 'sswd'}, inplace=True)
souss_massa.peak_Q = souss_massa.pwd #This is a work around so we are able to calculat all energy demand with the SW_pumping method
souss_massa.avg_Q = souss_massa.sswd #This is a work around so we are able to calculat all energy demand with the SW_pumping method
souss_massa.swpp_e = 'swpp_e'
souss_massa.pd_e = 'swpp_e'
souss_massa.swpa_e = 'swpa_e'
souss_massa.pumping_hours_per_day = 10
souss_massa.SWpump_eff = 0.6

#Defines the PWD. It is defined as double the seasonal demand for agricultural sites
souss_massa.df[souss_massa.pwd] = souss_massa.df[souss_massa.sswd] / 3600 / 30 / souss_massa.pumping_hours_per_day #to convert to cubic meter per second [m3/s]
souss_massa.df.loc[souss_massa.df['type']=='Agriculture', souss_massa.pwd] *= 2

#Calculates some required parameters
souss_massa.get_A(inplace=True)
souss_massa.get_V(inplace=True, axis=0)
souss_massa.get_Re(inplace=True, axis=0)
souss_massa.get_f(inplace=True, axis=0)

souss_massa.get_sw_tdh(inplace = True, axis=0) #this is called sw but it calculets for gw too. I still need to change the names
souss_massa.get_SWpumping_energy(inplace = True, axis=0) #the same here

souss_massa.df.loc[souss_massa.df.swpp_e<0, souss_massa.swpp_e] = 0
souss_massa.df.loc[souss_massa.df.swpa_e<0, souss_massa.swpa_e] = 0

#Define energy intensity for seawater desalination project
desalination_energy_int = 3.31 # kWh/m3
#We compute the energy demand for deslination multiplying the monthly water requirement by the energy intensity, 
#and add it to the current energy requirements for desalinated water pumping
sm_desal = nexus_tool.Model(souss_massa.df.loc[souss_massa.df['type'].str.contains('DS')].copy())

sm_desal.df[souss_massa.swpa_e] = sm_desal.df[souss_massa.sswd] * desalination_energy_int
#Then we divide the total energy requiered for desalination by the daily pumping hours and the days of the month
sm_desal.df[souss_massa.swpp_e] = desalination_energy_int

#Here we load the WWTP inflow data
file_path = os.path.join(load_data, 'wwtp_inflow.gz')
df_wwtp = pd.read_csv(file_path)

#We define an energy intensity for wastewater treatment and compute the energy demand
wwtp_energy_int = 0.6 # kWh/m3
df_wwtp['swpa_e'] = df_wwtp.value * wwtp_energy_int

souss_massa.df.to_csv(os.path.join(output_main_folder, 'results.gz'), index=False)
sm_desal.df.to_csv(os.path.join(output_main_folder, 'desal_data.gz'), index=False)
df_wwtp.to_csv(os.path.join(output_main_folder, 'wwtp_data.gz'), index=False)