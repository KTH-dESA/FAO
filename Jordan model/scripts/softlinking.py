import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
import nexus_tool.weap_tools as wp
import os
import scripts.softlinking_functions as sf

## Read processed schematic files

gis_folder = os.path.join('data', 'GIS', 'Processed layers')
demand = gpd.read_file(os.path.join(gis_folder, 'Demand_points.gpkg'))
supply = gpd.read_file(os.path.join(gis_folder, 'Supply_points.gpkg'))
pipelines = gpd.read_file(os.path.join(gis_folder, 'Pipelines.gpkg'))

## Read WEAP input data file

data_file = os.path.join('Data', 'WEAP Results', 'March 2021', 'WEAPResults - Reference.xlsx')
data = pd.ExcelFile(data_file)

## Read required water data sheets

sheet_names = {'Supply Requirement_AG': 'Agriculture', 
               'Supply Requirement_Muni': 'Municipality', 
               'Supply Requirement_Ind': 'Industry'}
required_demand = sf.get_data(sheet_names, data, demand, 'point', '{}')

## Read delivered water data sheets

sheet_names = {'Supply Delivered_AG': 'Agriculture', 
               'Supply Delivered_Muni': 'Municipality', 
               'Supply Delivered_Ind': 'Industry'}
delivered_demand = sf.get_data(sheet_names, data, demand, 'point', '{}')

delivered_demand['point'] = delivered_demand['point'].str.replace('BALQA', 'BQ')
delivered_demand['point'] = delivered_demand['point'].str.replace('MADABA', 'MD')
delivered_demand['point'] = delivered_demand['point'].str.replace('JERASH', 'JA')

governorates = {'AJ': 'Ajlun', 'AM': 'Amman', 'AQ': 'Aqaba', 'BQ': 'Balqa', 
                'IR': 'Irbid', 'JA': 'Jarash', 'KA': 'Karak', 'MA': 'Ma`an', 
                'MD': 'Madaba', 'MF': 'Mafraq', 'TA': 'Tafilah', 'ZA': 'Zarqa'}
delivered_demand['Governorate'] = [governorates[i[0:2]] for i in delivered_demand['point']]

## Mapping elevation from spatial layer to dataframe

delivered_demand['elevation_m'] = delivered_demand.point.map(demand.groupby('point')['elevation_m'].mean())

## Read desalination plants data sheets

#Red-Dead project
sheet_names = {'RedDead': 'Desalination'}
red_dead = sf.get_data(sheet_names, data, supply, 'point', '{}')
red_dead = red_dead.append(sf.get_data(sheet_names, data, demand, 'point', '{}'))
red_dead['point'] = 'RedDead'

#Aqaba plant
sheet_names = {'Aqaba Desal': 'Desalination'}
aqaba_desal = sf.get_data(sheet_names, data, demand, 'point', '{}')
aqaba_desal['point'] = 'Aqaba Desal'

#Merge both together
desalination = red_dead.append(aqaba_desal)
desalination['value'] = abs(desalination['value'])
desalination['variable'] = desalination.variable.str.strip()

## Read groundwater supply data sheet

sheet_names = {'GW Pumping': 'Groundwater supply'}
groundwater = supply.loc[supply['type']=='Groundwater supply']
gw_supply =  sf.get_data(sheet_names, data, groundwater, 'point', '{}')

## Read and process groundwater level change

sheet_names = {'Groundwater': 'Thickness'}
gw_thickness =  sf.get_data(sheet_names, data, supply, 'point', '{}')
gw_thickness.rename(columns={'value': 'thickness', 'units': 'thickness_units'}, inplace=True)
gw_thickness.drop(columns='type', inplace=True)

for point in gw_thickness.point.unique():
    _filter = (gw_thickness.point==point)
    init_thickness = gw_thickness.loc[_filter].iloc[0].thickness
    gw_thickness.loc[_filter, 'wtd_m'] = init_thickness - gw_thickness.loc[_filter, 'thickness'] + \
                                         float(groundwater.loc[groundwater.point==point,'wtd_m'].mean())

#Merge datasets on groundwater supply and level change (thickness)
gw_supply = gw_supply.merge(gw_thickness, on=['Year','Month','point','variable'])

## Read and process pipeline supply and wastewater treatment plants data sheets

sheet_names = {'Wadis': 'River/pipeline supply'}
surface_water =  sf.get_data(sheet_names, data, supply, 'point', '{}')

sheet_names = {'WWTP Inflow': 'WWTP'}
wwtp_inflow =  sf.get_data(sheet_names, data, supply, 'point', '{}')

## Pipeline flow processing

sheet_names = {'Pipelines': 'conveyance', 
               'PumpStations': 'conveyance'}
pl_flow = sf.get_data(sheet_names, data, pipelines, 'pipeline', '^{} [0-9]')

pl_flow = pl_flow.loc[~pl_flow.variable.str.contains('FR')]

pl_flow['point'] = np.nan
_vec = ~pl_flow.variable.isin(['Headflow','Reach'])
pl_flow.loc[_vec,'point'] = pl_flow.loc[_vec,'variable']

_df = pl_flow.loc[(pl_flow.Year==2020) & (pl_flow.Month==1)].groupby('pipeline').count()
idx = _df.loc[_df.point<3].index
_df = pl_flow.loc[pl_flow.pipeline.isin(idx)].copy()

for pipeline in _df.pipeline.unique():
    for row in pipelines.loc[pipelines.pipeline==pipeline].iterrows():
        if (row[1].point not in _df.loc[(_df.pipeline==pipeline), 'variable'].unique()) & (row[1].type=='Diversion Outflow'):
            pl_flow.loc[(pl_flow.pipeline==pipeline) & (pl_flow.variable=='Headflow'), 'point'] = row[1].point
        elif (row[1].point not in _df.loc[(_df.pipeline==pipeline), 'variable'].unique()) & (row[1].type=='Tributary inflow'):
                pl_flow.loc[(pl_flow.pipeline==pipeline) & (pl_flow.variable=='Reach'), 'point'] = row[1].point

_df = pipelines[['pipeline','point','index']].groupby(['pipeline','point']).mean()
_pl_flow = pl_flow.set_index(['pipeline','point'])
_index = _pl_flow.index.map(_df['index'].to_dict())
pl_flow['pl_index'] = [(i + 0.5) if i%1==0.5 else (i) for i in _index]
pl_flow['elevation'] = pl_flow.point.map(pipelines.groupby('point')['elevation_m'].mean().to_dict())
pl_flow['segment_length'] = pl_flow.pl_index.map(pipelines.groupby('index')['segment_length_m'].mean().to_dict())
pl_flow['pipeline_length'] = pl_flow.pipeline.map(pipelines.groupby('pipeline')['pl_length_m'].mean().to_dict())

pl_flow.dropna(subset=['point'], inplace=True)

sf.enumerate_segments(pl_flow)
sf.get_elevation_delta(pl_flow)

pl_flow.loc[(pl_flow.variable=='Reach') & (pl_flow.elevation_delta!=0), 'point'] = np.nan
pl_flow.dropna(subset=['point'], inplace=True)

sf.enumerate_segments(pl_flow)
sf.get_elevation_delta(pl_flow)

_point = supply.loc[supply['type']=='River/pipeline supply'].point.unique()

_pipe = pl_flow.pipeline.unique()
pl_flow['water_use'] = 0
for pipe in _pipe:
    n = pl_flow.loc[(pl_flow['pipeline']==pipe)].n.unique()
    for _n in range(1,n.max()+1):
        value2 = np.array(pl_flow.loc[(pl_flow['pipeline']==pipe) & (pl_flow['n']==_n), 'value'])
        value1 = np.array(pl_flow.loc[(pl_flow['pipeline']==pipe) & (pl_flow['n']==(_n-1)), 'value'])
        pl_flow.loc[(pl_flow['pipeline']==pipe) & (pl_flow['n']==_n), 'water_use'] = abs(value1 - value2)

    if pl_flow.loc[(pl_flow['pipeline']==pipe)].water_use.sum() == 0:
        pl_flow.loc[(pl_flow['pipeline']==pipe), 'water_use'] = \
                                        pl_flow.loc[(pl_flow['pipeline']==pipe), 'value'].mean()/\
                                        pl_flow.loc[(pl_flow['pipeline']==pipe), 'value'].count()

pl_flow['type'] = np.nan
pl_flow.loc[pl_flow['point'].isin(_point), 'type'] = 'River/pipeline supply'

## Read crop production data sheet

sheet_names = {'Production': 'Crop production'}
crop_production =  sf.get_data(sheet_names, data, demand, 'point', '^{}', {'Unnamed: 0': 'Year'}, ['Year'])

crop_production['variable'] = crop_production['variable'].str.replace('Summer ', '')
crop_production['variable'] = crop_production['variable'].str.replace('Winter ', '')
crop_production['Governorate'] = [governorates[i[0:2]] for i in crop_production['point']]

## Save the processed results

output_folder = os.path.join('Data', 'Processed results', 'Reference', 'Climate Change', 'level_1')
os.makedirs(output_folder, exist_ok=True)

sf.save_dataframe(desalination, 2020, 2050, output_folder, 'desalination.csv')
sf.save_dataframe(wwtp_inflow, 2020, 2050, output_folder, 'wwtp_inflow.csv')
# sf.save_dataframe(surface_water, 2020, 2050, output_folder, 'surface_water_supply.csv')
sf.save_dataframe(gw_supply, 2020, 2050, output_folder, 'groundwater_supply.csv')
sf.save_dataframe(pl_flow, 2020, 2050, output_folder, 'pipelines_flow.csv')

output_folder = os.path.join('..', 'Jordan dashboard', 'data_test', 'Reference', 'Climate Change', 'level_1')
os.makedirs(output_folder, exist_ok=True)

sf.save_dataframe(crop_production, 2020, 2050, output_folder, 'crop_production.csv')
sf.save_dataframe(delivered_demand, 2020, 2050, output_folder, 'water_delivered.csv')
sf.save_dataframe(required_demand, 2020, 2050, output_folder, 'water_requirements.csv')