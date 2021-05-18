#!/usr/bin/env python
# coding: utf-8

# # NEXUS tool: case study for the Souss-Massa basin - energy demand  calculations
# In this notebook a case study for the Souss-Massa basin is covered using the `nexustool` package. The water requirements for agricultural irrigation and domestic use were previously calculated using the Water Evaluation and Planning System (WEAP) model. In this case study, the energy requirements for groundwater pumping, wastewater treatment, desalination of seawater and pumping energy for water conveyance are estimated.
# 
# First import the package by running the following block:

# In[ ]:


import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import os
import nexustool
import pandas as pd


# ## 1. Read scenario data
# After importing all required packages, the input GIS data is loaded into the variable `df`. Change the `data_folder`, `scenario`, `climate` and `level` variables to reflect the name and relative location of your data file. This dataset should already have the water demand for irrigation results.

# In[ ]:

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
demand_data = str(snakemake.input.demand_data)
wwtp_inflow = str(snakemake.input.wwtp_inflow)
output = str(snakemake.output)

results_folder = output.split(os.path.basename(output))[0]

os.makedirs(results_folder, exist_ok = True)


# ## 2. Create nexus model 
# To create a model simply create an instance of the `nexus_tool.Model()` class and store it in a variable name. The `nexus_tool.Model()` class requires a dataframe as input data. Several other properties and parameter values can be defined by explicitly passing values to them. To see a full list of parameters and their explaination refer to the documentation of the package. We wil create a model using the `pipelines_flow.csv` data:

# In[ ]:


#Define the path to read the scenario input data and reads it in
df = pd.read_csv(demand_data)

#Creates the nexus model with the input dataframe
souss_massa = nexustool.Model(df)


# ## 3. Define variable names
# The names of the properties of the model can be changed at any time. This is important for the model to know how each property is called withing your input data. To check the current property names run the `.print_properties()` method, a list with the names of each property and its current value will be displayed.
# 
# Then you can provide the right names for each property, calling them and assigning a value as:
# ```python
# souss_massa.elevation_diff = 'elevation_delta'
# souss_massa.gw_depth = 'name_of_ground_water_depth'
# ```
# 
# In this particular case we will need to change the following default values:

# In[ ]:


souss_massa.elevation_diff = 'elevation_diff' #for the case of GW, the elevation_diff is set to be the wtd
souss_massa.L = 'distance' #for the case of GW, the distance is set to be the wtd
souss_massa.D = 'Pipe_diameter'

#Defines the name of the variable for Peak Water Demand and Seasonal Scheme Water demand (monthly)
souss_massa.pwd = 'pwd' # Peak Water Demand
souss_massa.sswd = 'sswd' # Seassonal Scheme Water Demand
souss_massa.df.rename(columns={'value': 'sswd'}, inplace=True) #Renames the name of the column value to sswd 
souss_massa.pp_e = 'pp_e' # Peak Pumping Energy
souss_massa.pa_e = 'pa_e' # Pumping Average Energy


# ## 4. Define pipelines diameters and average pumping hours, pumping efficiency
# Now we need to define the specifications of the water network, giving pipeline / canal diameter values:

# In[ ]:


souss_massa.df['Pipe_diameter'] = 1
souss_massa.df.loc[souss_massa.df['type'].str.contains('GW'), 'Pipe_diameter'] = 1000 # in here we use a large number to disregar friction losses in groundwater pumping
souss_massa.df.loc[souss_massa.df['type'].str.contains('DS'), 'Pipe_diameter'] = 1.5
souss_massa.df.loc[souss_massa.df['type'].str.contains('Pipeline'), 'Pipe_diameter'] = 1.7

souss_massa.pumping_hours_per_day = 10
souss_massa.pump_eff = 0.6


# ## 5. Peak Water Demand (PWD)
# The $PWD$ is definfe as the daily peak cubic meters of water pumped per second withing the month. To accomplish that, the $SSWD$ (m<sup>3</sup>/month) is divided by 30 days per month, 3600 seconds per hour and the amount of average pumping hours in a day. This provides the $PWD$ in m<sup>3</sup>/s:
# 
# $$
# PWD\,(m^3/s) = \frac{SSWD\,(m^3/month)}{30\,(day/month)\cdot PumpHours\,(h/day)\cdot 3600\, (s/h)}
# $$
# 
# Moreover, the $PWD$ for agricultural irrigation is assumed as double the normal $PWD$. We make this calculations as per the following cell:

# In[ ]:


#Defines the PWD. It is defined as double the seasonal demand for agricultural sites
souss_massa.df[souss_massa.pwd] = souss_massa.df[souss_massa.sswd] / 30 / souss_massa.pumping_hours_per_day / 3600 #to convert to cubic meter per second [m3/s]
souss_massa.df.loc[souss_massa.df['type']=='Agriculture', souss_massa.pwd] *= 2


# ## 6. Calculate pumping energy requirements
# To estimate the pumping energy requirements for conveyance, first we need to calculate the Total Dinamic Head (TDH). This, is a measure in meters that accounts for the elevation difference between two points and the pressure loss in distribution.
# 
# For that, the area $A$ `.pipe_area()`, the velocity $V$ `.flow_velocity()`, the Reynolds number $Re$ `.reynolds()` and the friction factor $f$ `.friction_factor()` need to be estimated. The `nexustool` provides simple functions that allows us make an easy estimation of these variables, which have the following formulas implemented in the background:
# 
# $$
# A\,(m^2) = \pi\cdot \frac{D^2}{4}
# $$
# 
# $$
# V\,(m/s) = \frac{SSWD\,(m^3/month)}{PumpHours\,(h/day)\cdot 30\,(day/month)\cdot 3600\,(s/h)\cdot A\,(m^2)}
# $$
# 
# $$
# Re = \frac{V\,(m/s)\cdot D\,(m)}{v\,(m^2/s)}
# $$
# 
# Where $v$ is the kinematic viscosity of water at around 1.004e-06 m<sup>2</sup>/s. And the frction factor is estimated according to the Swamee–Jain equation:
# 
# $$
# f = \frac{0.25}{\left[log_{10}\left(\frac{\epsilon}{3.7D}+\frac{5.74}{Re^{0.9}}\right)\right]^2}
# $$
# 
# Where $\epsilon$ is the roughness of the material. 

# In[ ]:


souss_massa.pipe_area() 
souss_massa.flow_velocity()
souss_massa.reynolds()
souss_massa.friction_factor()


# Then, the TDH can be calculated by simply calling the `.get_tdh()` function.
# 
# $$
# TDH\,(m) = f\cdot \frac{L\,(m)}{D\,(m)}\cdot \frac{V(m/s)^2}{2\cdot g\,(m/s^2)}
# $$
# 
# Whereas the conveyance pumping energy requirements by calling the `.get_pumping_energy()` method. The equation used to calculate the Electricity Demand ($E_D$) for pumping is as follows:
# 
# $$
# E_D\,(kW_h) = \frac{SSWD\,(m^3)\cdot \rho\,(kg/m^3)\cdot g\,(m/s^2)\cdot TDH\,(m)}{PP_{eff}\,(\%)\cdot 3600\,(s/h)\cdot 1000\,(W/kW)}
# $$
# 
# The variable withing the Model for the $E_D$ is the `pa_e` or Pumping Average Electricity requirements.
# 
# Moreover, the Power Demand for pumping ($PD$) is denoted by the variable `pp_e` and calculated by the following formula:
# 
# $$
# PD\,(kW) = \frac{PWD\,(m^3/s)\cdot \rho\,(kg/m^3)\cdot g\,(m/s^2)\cdot TDH\,(m)}{PP_{eff}\,(\%)\cdot 1000\,(W/kW)}
# $$
# 
# The `.get_pumping_energy()` method calculates both the $E_D$ (`pa_e`) and $PD$ (`pp_e`).

# In[ ]:


souss_massa.get_tdh()
souss_massa.get_pumping_energy()

souss_massa.df.loc[souss_massa.df.pp_e<0, souss_massa.pp_e] = 0 # ensures no negative energy values are considered
souss_massa.df.loc[souss_massa.df.pa_e<0, souss_massa.pa_e] = 0 # ensures no negative power values are considered

# We exclude energy for pumping calculations done for the Complexe Aoulouz Mokhtar Soussi, 
# as this pipeline is known to be driven by gravity only
souss_massa.df.loc[souss_massa.df['Supply point'].str.contains('Complexe Aoulouz Mokhtar Soussi'), 'pa_e'] = None


# ## 7. Calculating desalination energy requirements
# Desalination energy requirements are estimated by multipliying the monthly average desalinated water (`sswd`), by an energy intensity factor (`desal_energy_int`) based on the characteristics of the desalination plant.

# In[ ]:


#Define energy intensity for seawater desalination project
desal_energy_int = 3.31 # kWh/m3

#Create a new nexus Model with the data relevant to the desalination plant only, filtering by the key work DS (Desalination)
sm_desal = nexustool.Model(souss_massa.df.loc[souss_massa.df['type'].str.contains('DS')].copy())

#Multiply the sswd by the energy intensity for treatment
sm_desal.df[souss_massa.pa_e] = sm_desal.df[souss_massa.sswd] * desal_energy_int


# ## 8. Calculating wastewater treatment energy requirements
# Wastewater treatment energy is dependent on the type of treatment required. Wastewater treatment can be subdivided into three stages: primary, secondary and tertiary. The treatment stages used, are then dependent on the final quality requirements of the treated wastewater. Thus, for wastewater that will be treated and returned to the ecosystem, often primary to secondary treatment is enough. On the other hand, treated wastewater intended for agricultural irrigation or drinking purposes, should go through secondary to terciary treatment to ensure proper desinfecton levels. 
# 
# Depending on the scenario run, we will need then to use a proper wastewater treatment energy intensity. In general, the higher the number of stages, the higher the energy requirements. In this model, we used an energy intensity of **0.1 kWh/m<sup>3</sup>** for treated wastewater that is not being reused, and **0.8 kWh/m<sup>3</sup>** for treated wastewater reused in agricultural irrigation.

# In[ ]:


#Here we load the WWTP inflow data
df_wwtp = pd.read_csv(wwtp_inflow)

#We define an energy intensity for wastewater treatment and compute the energy demand
if scenario in ['Reference Wastewater Reuse', 'Desalination Wastewater Reuse',
                'Integrated Strategies', 'Green Generation']:
    wwtp_energy_int = 0.8 # kWh/m3
else:
    wwtp_energy_int = 0.1 # kWh/m3
df_wwtp['pa_e'] = df_wwtp.value * wwtp_energy_int


# ## 9. Saving the results
# Finally, we save the resulting dataframes as `.gz` files, which is a compressed version of a `csv` file:

# In[ ]:


#Define and create the output folder
os.makedirs(results_folder, exist_ok=True)

#Save the results
souss_massa.df.to_csv(os.path.join(results_folder, 'results.gz'), index=False)
sm_desal.df.to_csv(os.path.join(results_folder, 'desal_data.gz'), index=False)
df_wwtp.to_csv(os.path.join(results_folder, 'wwtp_data.gz'), index=False)
