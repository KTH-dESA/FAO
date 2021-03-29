import sys, os
sys.path.append("..") #this is to add the avobe folder to the package directory
import nexustool
import pandas as pd

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
demand_data = str(snakemake.input.demand_data)
wwtp_inflow = str(snakemake.input.wwtp_inflow)
output = str(snakemake.output)

results_folder = output.split(os.path.basename(output))[0]

os.makedirs(results_folder, exist_ok = True)

#Define the path to read the scenario input data and reads it in
df = nexus_tool.read_csv(demand_data)

#Creates the nexus model with the input dataframe
souss_massa = nexus_tool.Model(df)

#Define the diameter of the pipelines. The first is general for all and the following specific for each case
souss_massa.df['Pipe_diameter'] = 1
souss_massa.df.loc[souss_massa.df['type'].str.contains('GW'), 'Pipe_diameter'] = 1000
souss_massa.df.loc[souss_massa.df['type'].str.contains('Pipeline'), 'Pipe_diameter'] = 1.2
souss_massa.df.loc[souss_massa.df['type'].str.contains('DS'), 'Pipe_diameter'] = 1

#Define the variable to take into account for elevation difference and lenght of pipelines
souss_massa.elevation = 'elevation_diff' #for the case of GW, the elevation_diff is set to be the wtd
souss_massa.L = 'distance' #for the case of GW, the distance is set to be the wtd

#Defines the name of the variable for Peak Water Demand and Seasonal Water demand (monthly)
souss_massa.elevation_diff = 'elevation_diff' #for the case of GW, the elevation_diff is set to be the wtd
souss_massa.L = 'distance' #for the case of GW, the distance is set to be the wtd

#Defines the name of the variable for Peak Water Demand and Seasonal Scheme Water demand (monthly)
souss_massa.pwd = 'pwd'
souss_massa.sswd = 'sswd'
souss_massa.df.rename(columns={'value': 'sswd'}, inplace=True) #Renames the name of the column value to sswd 
souss_massa.pp_e = 'pp_e' #
souss_massa.pa_e = 'pa_e'
souss_massa.pumping_hours_per_day = 10
souss_massa.SWpump_eff = 0.6

#Defines the PWD. It is defined as double the seasonal demand for agricultural sites
souss_massa.df[souss_massa.pwd] = souss_massa.df[souss_massa.sswd] / 30 / souss_massa.pumping_hours_per_day / 3600 #to convert to cubic meter per second [m3/s]
souss_massa.df.loc[souss_massa.df['type']=='Agriculture', souss_massa.pwd] *= 2

#Calculates some required parameters
souss_massa.pipe_area()
souss_massa.flow_velocity()
souss_massa.reynolds()
souss_massa.friction_factor()

souss_massa.get_tdh()
souss_massa.get_pumping_energy()

souss_massa.df.loc[souss_massa.df.pp_e<0, souss_massa.pp_e] = 0 # ensures no negative energy values are considered
souss_massa.df.loc[souss_massa.df.pa_e<0, souss_massa.pa_e] = 0 # ensures no negative power values are considered

#Define energy intensity for seawater desalination project
desal_energy_int = 3.31 # kWh/m3

#Create a new nexus Model with the data relevant to the desalination plant only, filtering by the key work DS (Desalination)
sm_desal = nexustool.Model(souss_massa.df.loc[souss_massa.df['type'].str.contains('DS')].copy())

#Multiply the sswd by the energy intensity for treatment
sm_desal.df[souss_massa.pa_e] = sm_desal.df[souss_massa.sswd] * desal_energy_int

#Here we load the WWTP inflow data
# file_path = os.path.join(load_folder, 'wwtp_inflow.gz')
df_wwtp = pd.read_csv(wwtp_inflow)

#We define an energy intensity for wastewater treatment and compute the energy demand
if 'Reuse' in scenario:
    wwtp_energy_int = 0.8 # kWh/m3
    df_wwtp['pa_e'] = df_wwtp.value * wwtp_energy_int
    df_wwtp.loc[df_wwtp['point']=='Agadir WWTP', 'pa_e'] = df_wwtp.loc[df_wwtp['point']=='Agadir WWTP'].value * 0.6
else:
    wwtp_energy_int = 0.1 # kWh/m3
    df_wwtp['pa_e'] = df_wwtp.value * wwtp_energy_int
    
souss_massa.df.loc[souss_massa.df['type'].str.contains('WWR'), 'pa_e'] = None
souss_massa.df.loc[souss_massa.df['type'].str.contains('SW'), 'pa_e'] = None
souss_massa.df.loc[souss_massa.df['Supply point'].str.contains('Complexe Aoulouz Mokhtar Soussi'), 'pa_e'] = None

souss_massa.df.to_csv(output, index=False)
sm_desal.df.to_csv(os.path.join(results_folder, 'desal_data.gz'), index=False)
df_wwtp.to_csv(os.path.join(results_folder, 'wwtp_data.gz'), index=False)