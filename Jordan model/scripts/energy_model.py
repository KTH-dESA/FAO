# NEXUS tool: case study for the Jordan - energy demand calculations

import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import os
import nexus_tool
import pandas as pd

## 1. Read scenario data

pipelines_flow = snakemake.input.pipelines_flow
input_folder = pipelines_flow.split(os.path.basename(pipelines_flow))[0]

## 2. Create nexus model 

file_path = pipelines_flow
df = nexus_tool.read_csv(file_path)
jordan = nexus_tool.Model(df)

## 3. Define variable names

jordan.elevation = 'elevation_delta'
jordan.L = 'segment_length'
jordan.peak_Q = 'value'
jordan.avg_Q = 'value'

## 4. Define pipelines diameters and pumping efficiency

jordan.df['Pipe_diameter'] = 1 #in m
jordan.df.loc[(df['pipeline']=='KAC'), 'Pipe_diameter'] = 4 #in m
jordan.df.loc[(df['pipeline']=='PL_Disi'), 'Pipe_diameter'] = 1.4 #in m
jordan.df.loc[(df['pipeline']=='PL_RedDead'), 'Pipe_diameter'] = 1.4 #in m
jordan.df.loc[(df['pipeline']=='PS_ZaraMain'), 'Pipe_diameter'] = 1.8 #in m
jordan.df.loc[(df['pipeline']=='PL_KACToZay'), 'Pipe_diameter'] = 1.4 #in m
jordan.df.loc[(df['pipeline']=='PL_Dab_AinGhazal'), 'Pipe_diameter'] = 1.4 #in m
jordan.df.loc[(df['pipeline']=='PL_ZaytoDabouq'), 'Pipe_diameter'] = 1.4 #in m

jordan.SWpump_eff = 1 # pumping efficiency

## 5. Calculate pumping energy requirements

jordan.get_A(inplace=True)
jordan.get_V(inplace=True, axis=0)
jordan.get_Re(inplace=True, axis=0)
jordan.get_f(inplace=True, axis=0)

jordan.get_sw_tdh(inplace = True, axis=0)
jordan.get_SWpumping_energy(inplace = True, axis=0)

## 6. Create nexus model for groundwater pumping

file_path = os.path.join(input_folder, 'groundwater_supply.gz')
df_groundwater = pd.read_csv(file_path)

jordan_gw = nexus_tool.Model(df_groundwater)

jordan_gw.df['Pipe_diameter'] = 0.4
jordan_gw.elevation = 'wtd_m'
jordan_gw.L = 'wtd_m'
jordan_gw.peak_Q = 'value'
jordan_gw.avg_Q = 'value'

## 7. Calcuculate groundwater pumping energy

jordan_gw.get_A(inplace=True)
jordan_gw.get_V(inplace=True, axis=0)
jordan_gw.get_Re(inplace=True, axis=0)
jordan_gw.get_f(inplace=True, axis=0)

jordan_gw.get_sw_tdh(inplace = True, axis=0)
jordan_gw.get_SWpumping_energy(inplace = True, axis=0)

## 8. Calculate wastewater treatment energy

file_path = os.path.join(input_folder, 'wwtp_inflow.gz')
df_wwtp = pd.read_csv(file_path)

wwtp_energy_int = 0.6 # kWh/m3
df_wwtp['SWPA_E_'] = df_wwtp.value * wwtp_energy_int

## 9. Calculate desalination energy

file_path = os.path.join(input_folder, 'desalination.gz')
df_desal = pd.read_csv(file_path)

red_dead_energy_int = 3.31 # kWh/m3
aqaba_energy_int = 5 # kWh/m3
df_desal['SWPA_E_'] = 0

df_desal.loc[df_desal.point=='RedDead', 'SWPA_E_'] = df_desal.loc[df_desal.point=='RedDead', 'value'] * red_dead_energy_int
df_desal.loc[df_desal.point=='Aqaba Desal', 'SWPA_E_'] = df_desal.loc[df_desal.point=='Aqaba Desal', 'value'] * aqaba_energy_int

## 10. Save result files

pipelines_data = snakemake.output.pipelines_data
results_folder = pipelines_data.split(os.path.basename(pipelines_data))[0]
os.makedirs(results_folder, exist_ok=True)

jordan.df.to_csv(pipelines_data, index=False)
jordan_gw.df.to_csv(os.path.join(results_folder, 'groundwater_pumping.gz'), index=False)
df_wwtp.to_csv(os.path.join(results_folder, 'wwtp_data.gz'), index=False)
df_desal.to_csv(os.path.join(results_folder, 'desal_data.gz'), index=False)