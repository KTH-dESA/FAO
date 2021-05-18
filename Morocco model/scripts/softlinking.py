import os
import pandas as pd
import geopandas as gpd
import numpy as np
from softlinking_functions import integrate_data, data_merging

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
data_file = str(snakemake.input.data)
demand_data_output = str(snakemake.output.demand_data)
wwtp_inflow_output = str(snakemake.output.wwtp_inflow)
production_data_output = str(snakemake.output.production_data)

demand_points_path = str(snakemake.input.demand_points)
spatial_folder = demand_points_path.split(os.path.basename(demand_points_path))[0]

demand_links = gpd.read_file(demand_points_path)
supply_links = gpd.read_file(os.path.join(spatial_folder, 'supply_points.gpkg'))
wwtp = gpd.read_file(os.path.join(spatial_folder, 'wwtp.gpkg'))
diversions = gpd.read_file(os.path.join(spatial_folder, 'pipelines.gpkg'))
groundwater = gpd.read_file(os.path.join(spatial_folder, 'groundwater.gpkg'))
all_points = gpd.read_file(os.path.join(spatial_folder, 'all_points.gpkg'))

output_folder = demand_data_output.split(os.path.basename(demand_data_output))[0]

## Defining the years of the simulation

init_year = 2020
end_year = 2050

## Merging WEAP results for water supply with spatial data (what is actually being delivered)

MerchidSudMoroc = 26192
    
data = pd.ExcelFile(data_file)

sheet_names = {'Desalination': 'DS Agriculture', 
               'GP Irrigation': 'GW Agriculture', 
               'GP Domestic': 'GW Domestic',
               'MAR': 'Aquifer recharge',
               'SW Irrigation': 'SW Agriculture',
               'SW Domestic': 'SW Domestic',
               'Wastewater Reuse AG': 'WWR Agriculture',
               'Wastewater Reuse Agadir': 'WWR Domestic'}

df = pd.DataFrame()
for sheet_name, category in sheet_names.items():
    df = df.append(integrate_data(data, sheet_name, category, {'Demand point': demand_links, 
                                                         'Supply point': supply_links}, 
                                  demand_links, init_year, end_year))

## Extracting water table depth data to groundwater supply points

df['wtd'] = df['Supply point'].map(groundwater.set_index('point')['wtd_m'])

wtd_change = integrate_data(data,'GW Change in Elev', 'GW wtd', {'GW': groundwater}, 
                            demand_links, init_year, end_year, 'point', 'point')
df.set_index(['Date', 'Supply point'], inplace=True)
df['wtd'] -= df.index.map(wtd_change.set_index(['Date', 'point']).value)
df.reset_index(inplace=True)

## Calculating elevation difference between demand and supply points

df['elevation_diff'] = df.links.map(demand_links.set_index('links').elevation) - \
                       df.links.map(supply_links.set_index('links').elevation)

dff = df.loc[df['Supply point'].isin(diversions['Demand point'].unique())].groupby(['Date','Supply point','Year','Month']).agg({'value': 'sum'}).reset_index()
dff.rename(columns={'Supply point': 'Demand point'}, inplace=True)
dff['Supply point'] = dff['Demand point'].map(diversions.set_index('Demand point')['Supply point'])
dff['elevation_diff'] = dff['Demand point'].map(diversions.set_index('Demand point')['elevation_diff'])
dff['Province'] = dff['Supply point'].map(supply_links.drop_duplicates('point').set_index('point')['Province'])
dff['type'] = 'Transmission Pipeline'

df = df.append(dff, sort=False, ignore_index=True)
df.loc[df.Province.isna(),'Province'] = df['Demand point'].map(demand_links.drop_duplicates('point').set_index('point')['Province'])

## Cleaning up the variables

supply = gpd.GeoDataFrame(geometry=df['Supply point'].map(all_points.set_index('point').geometry), crs='epsg:4326')

supply.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)
demand = gpd.GeoDataFrame(geometry=df['Demand point'].map(all_points.set_index('point').geometry), crs='epsg:4326')

demand.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)

df['distance'] = supply.distance(demand)

df.loc[df['type'].str.contains('GW'), 'distance'] = df.loc[df['type'].str.contains('GW'), 'wtd']
df.loc[df['type'].str.contains('GW'), 'elevation_diff'] = df.loc[df['type'].str.contains('GW'), 'wtd']

df.loc[df['type'].str.contains('SW'), 'distance'] = np.nan
df.loc[df['type'].str.contains('SW'), 'elevation_diff'] = np.nan

#define the desalinated water for agriculture and for domestic use
df.loc[(df['type']=='DS Agriculture')&(df['Demand point']=='Agadir'), 'type'] = 'DS Domestic'

## Merging WEAP results for wastewater treatment with spatial data

df_wwtp = integrate_data(data, 'WWTP Inflow', 'wwtp', {'WWTP': wwtp}, 
                         demand_links, init_year, end_year, 'point', 'point')
df_wwtp['Province'] = None
df_wwtp.loc[df_wwtp['point'].isin(['Agadir WWTP', 'Drargua WWTP']), 'Province'] = 'Agadir-Ida-Ou-Tanane'
df_wwtp.loc[df_wwtp['point'].isin(['Lqliaa WWTP', 'Biougra WWTP', 'Ait Baha WWTP']), 'Province'] = 'Chtouka-Ait Baha'
df_wwtp.loc[df_wwtp['point'].isin(['Ouled Teima WWTP', 'ELGuerdane WWTP',
                                   'Ait Iaaza WWTP', 'Oulad Berhil WWTP',
                                   'Aoulouz WWTP']), 'Province'] = 'Taroudannt'
df_wwtp.loc[df_wwtp['point'].isin(['Drargua WWTP']), 'Province'] = 'Inezgane-Ait Melloul'

## Merging WEAP results for water requirements with spatial data (what is needed)

sheet_names = {'AgWaterDemand': 'Agriculture', 
               'DomSupplyReq': 'Domestic'}

df_required = pd.DataFrame()
for sheet_name, category in sheet_names.items():
    df_required = df_required.append(integrate_data(data, sheet_name, category, {'Demand point': demand_links},
                                                    demand_links, init_year, end_year, 'point', 'point'))

df['water_required'] = df.set_index(['Date','Demand point']).index.map(df_required.set_index(['Date','point']).value)

df_unmet_month = 1 - (df.groupby(['Date', 'Demand point'])['value'].sum() / \
               df.groupby(['Date', 'Demand point'])['water_required'].mean())

df['unmet_demand_month'] = df.set_index(['Date','Demand point']).index.map(df_unmet_month)

water_req_year = df.groupby(['Year', 'Date', 'Demand point'])['water_required'].mean().reset_index().groupby(['Year', 'Demand point'])['water_required'].sum()

df_unmet_year = 1 - (df.groupby(['Year', 'Demand point'])['value'].sum() / \
                     water_req_year)

df['unmet_demand_year'] = df.set_index(['Year','Demand point']).index.map(df_unmet_year)


df.fillna({'unmet_demand_year': 0, 'unmet_demand_month': 0}, inplace=True)

## Merging WEAP results for crop productivity with spatial data

# Process crop production data
df_production = data.parse('Annual Crop Production', skiprows=3)
df_production.rename(columns={'Unnamed: 0': 'Year'}, inplace=True)
df_production.columns = df_production.columns.str.replace('"', '').str.strip()
df_production = df_production.loc[df_production.Year!='Sum']
df_production.drop(columns='Sum', inplace=True)
df_production = df_production.loc[(df_production.Year >= init_year) & (df_production.Year <= end_year)]
df_production = df_production.melt(id_vars=['Year'])
df_production['point'] = [row[0] for row in df_production['variable'].str.split('\\')]
df_production['crop'] = [row[-1] for row in df_production['variable'].str.split('\\')]
df_production.rename(columns={'value': 'production_kg'}, inplace=True)
df_production.drop(columns='variable', inplace=True)
df_production['crop'] = df_production['crop'].str.replace('_', ' ')

df_production.to_csv(production_data_output, index=False)
df.to_csv(demand_data_output, index=False)
df_wwtp.to_csv(wwtp_inflow_output, index=False)