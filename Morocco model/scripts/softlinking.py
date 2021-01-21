import os
import pandas as pd
import geopandas as gpd
import numpy as np

scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
data_file = str(snakemake.input.data)
demand_points = str(snakemake.input.demand_points)
demand_data_output = str(snakemake.output.demand_data)
wwtp_inflow_output = str(snakemake.output.wwtp_inflow)
production_data_output = str(snakemake.output.production_data)

# spatial_folder = f"{snakemake.params.dash_folder}/spatial_data"
spatial_folder = demand_points.split(os.path.basename(demand_points))[0]
output_folder = demand_data_output.split(os.path.basename(demand_data_output))[0]

# output_folder = output.split('demand_data.gz')[0]
# output_folder = os.path.join('Data', 'Processed results', scenario, climate)

def integrate_data(data, sheet_name, category, dff_dict, var_name='links', target='point'):
    df = data.parse(sheet_name, skiprows=3)
    df.rename(columns={'Unnamed: 0': 'Date'}, inplace=True)
    df.columns = df.columns.str.replace('"', '').str.strip()

    df.columns = df.columns.str.replace('Groundwater','GW')
    df.columns = df.columns.str.replace('Grounwater','GW')
    df.columns = df.columns.str.replace('GW of ','')
    df.columns = df.columns.str.replace('GW ','')
    df.columns = df.columns.str.replace('I_TRSPD','I_Traditional Rehabilite du Souss Perimetre Diffus')

    for link in demand_links.links:
        if np.array(df.columns[df.columns.str.contains(link)]).size > 0:
            df.rename(columns={df.columns[df.columns.str.contains(link)][0]: link}, inplace=True)

    df = df.loc[df.Date!='Sum']
    df.Date = pd.to_datetime(df.Date)
    df['Year'] = df.Date.dt.year
    df['Month'] = df.Date.dt.month
    
    drop_columns = []
    if 'Sum' in df.columns:
        drop_columns.append('Sum')
    df.drop(columns=drop_columns, inplace=True)

    df = df.melt(id_vars=['Date', 'Year', 'Month'])
    
    for name, dff in dff_dict.items():
        df_temp = dff.set_index(var_name)
        if var_name!=target:
            df[name] = df.variable.map(df_temp[target])
    
    df['type'] = category
    df.rename(columns={'variable': var_name}, inplace=True)
    if df.loc[~df[var_name].isin(df.dropna()[var_name].unique()),var_name].unique().size > 0:
        print("The following links were not found:")
        print(df.loc[~df[var_name].isin(df.dropna()[var_name].unique()),var_name].unique())
    return df

demand_links = gpd.read_file(demand_points)
supply_links = gpd.read_file(os.path.join(spatial_folder, 'Supply_points.gpkg'))
wwtp = gpd.read_file(os.path.join(spatial_folder, 'wwtp.gpkg'))
diversions = gpd.read_file(os.path.join(spatial_folder, 'Pipelines.gpkg'))
groundwater = gpd.read_file(os.path.join(spatial_folder, 'Groundwater.gpkg'))
all_points = gpd.read_file(os.path.join(spatial_folder, 'all_points.gpkg'))

MerchidSudMoroc = 26192

os.makedirs(output_folder, exist_ok = True)
    
# data_file = os.path.join('Data', 'WEAP Results', f'SoussMassa Results - {scenario} - {climate}.xlsx')
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
                                                         'Supply point': supply_links}))

df['wtd'] = df['Supply point'].map(groundwater.set_index('point')['wtd_m'])

wtd_change = integrate_data(data,'GW Change in Elev', 'GW wtd', {'GW': groundwater}, 'point', 'point')
df.set_index(['Date', 'Supply point'], inplace=True)
df['wtd'] -= df.index.map(wtd_change.set_index(['Date', 'point']).value)
df.reset_index(inplace=True)

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

supply = gpd.GeoDataFrame(geometry=df['Supply point'].map(all_points.set_index('point').geometry), crs='epsg:4326')
#             supply.dropna(inplace=True)
supply.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)
demand = gpd.GeoDataFrame(geometry=df['Demand point'].map(all_points.set_index('point').geometry), crs='epsg:4326')
#             demand.dropna(inplace=True)
demand.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)
df['distance'] = supply.distance(demand)

df.loc[df['type'].str.contains('GW'), 'distance'] = df.loc[df['type'].str.contains('GW'), 'wtd']
df.loc[df['type'].str.contains('GW'), 'elevation_diff'] = df.loc[df['type'].str.contains('GW'), 'wtd']

#temporal fix due to inacurate plecement of demand points
df.loc[df['type'].str.contains('SW'), 'distance'] = np.nan
df.loc[df['type'].str.contains('SW'), 'elevation_diff'] = np.nan

#define the delainated water for agriculture and for domestic use
df.loc[(df['type']=='DS Agriculture')&(df['Demand point']=='Agadir'), 'type'] = 'DS Domestic'

df_wwtp = integrate_data(data, 'WWTP Inflow', 'wwtp', {'WWTP': wwtp}, 'point', 'point')

sheet_names = {'AgWaterDemand': 'Agriculture', 
               'DomSupplyReq': 'Domestic'}

df_required = pd.DataFrame()
for sheet_name, category in sheet_names.items():
    df_required = df_required.append(integrate_data(data, sheet_name, category, {'Demand point': demand_links}, 'point', 'point'))

df['water_required'] = df.set_index(['Date','Demand point']).index.map(df_required.set_index(['Date','point']).value)

df_unmet_month = 1 - (df.groupby(['Date', 'Demand point'])['value'].sum() / \
               df.groupby(['Date', 'Demand point'])['water_required'].mean())

df['unmet_demand_month'] = df.set_index(['Date','Demand point']).index.map(df_unmet_month)

water_req_year = df.groupby(['Year', 'Date', 'Demand point'])['water_required'].mean().reset_index().groupby(['Year', 'Demand point'])['water_required'].sum()

df_unmet_year = 1 - (df.groupby(['Year', 'Demand point'])['value'].sum() / \
                     water_req_year)

df['unmet_demand_year'] = df.set_index(['Year','Demand point']).index.map(df_unmet_year)

# df.loc[df.unmet_demand_year<0, 'unmet_demand_year'] = 0
# df.loc[df.unmet_demand_month<0, 'unmet_demand_month'] = 0

df.fillna({'unmet_demand_year': 0, 'unmet_demand_month': 0}, inplace=True)

# Process crop production data
df_production = data.parse('Annual Crop Production', skiprows=3)
df_production.rename(columns={'Unnamed: 0': 'Year'}, inplace=True)
df_production.columns = df_production.columns.str.replace('"', '').str.strip()
df_production = df_production.loc[df_production.Year!='Sum']
df_production.drop(columns='Sum', inplace=True)
df_production = df_production.melt(id_vars=['Year'])
df_production['point'] = [row[0] for row in df_production['variable'].str.split('\\')]
df_production['crop'] = [row[-1] for row in df_production['variable'].str.split('\\')]
df_production['group'] = [row[1] for row in df_production['variable'].str.split('\\')]
df_production.rename(columns={'value': 'production_kg'}, inplace=True)
df_production.drop(columns='variable', inplace=True)
df_production['crop'] = df_production['crop'].str.replace('_', ' ')

df_production.to_csv(production_data_output, index=False)
df.to_csv(demand_data_output, index=False)
df_wwtp.to_csv(wwtp_inflow_output, index=False)

def data_merging(demand_points, supply_points, pipelines):
    df1 = demand_points.groupby('point').agg({'type': 'first',
                                              'geometry': 'first'}).reset_index()

    df2 = supply_points.groupby('point').agg({'type': 'first',
                                              'geometry': 'first'}).reset_index()

    df_pipelines = pipelines.groupby('diversion').agg({'geometry': 'first'}).reset_index()

    df = df1.append(df2, ignore_index=True)
    df['lon'] = [point.xy[0][0] for point in df.geometry]
    df['lat'] = [point.xy[1][0] for point in df.geometry]

    pipe_coords = pd.DataFrame({'lon': [], 'lat': []})
    for name, point in zip(df_pipelines.diversion, df_pipelines.geometry):
        lon = list(point.xy[0]) + [None]
        lat = list(point.xy[1]) + [None]
        df_temp = pd.DataFrame({'lon': lon, 'lat': lat})
        df_temp['name'] = name
        pipe_coords = pipe_coords.append(df_temp, ignore_index=True)

    pipe_coords['type'] = 'pipeline'
    return df, pipe_coords
    
points_coords, pipe_coords = data_merging(demand_links, supply_links, diversions)

points_coords.to_csv(os.path.join(spatial_folder, 'points_coords.csv'), index=False)
pipe_coords.to_csv(os.path.join(spatial_folder, 'pipe_coords.csv'), index=False)