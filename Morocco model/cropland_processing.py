import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
from nexus_tool.weap_tools import get_zonal_stats
import fiona
import itertools
import os
fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support
fiona.drvsupport.supported_drivers['KML'] = 'rw' # enable KML support

cropland = gpd.read_file(r"Data\GIS\Processed layers\cropland_2013_processed_20200227.csv", encoding='utf-8')
cropland.drop(columns=['WKT', 'gen2013', 'area'], inplace=True)
cropland.reset_index(inplace=True)
cropland.rename(columns={'area_ha': 'area_m2', 'index': 'Demand point', 'prov': 'province', 'wtd_mean': 'wtd'}, inplace=True)
cropland.columns = cropland.columns.str.replace('mean','')
cropland.columns = cropland.columns.str.replace('mea','')
for column in list(cropland)[1:-3]:
    cropland.loc[cropland[column]=='', column] = np.nan

builtarea = gpd.read_file(r"Data\GIS\Processed layers\builtup_2013_processed_20200218_clean.csv", encoding='utf-8')
builtarea.drop(columns=['WKT', 'fid', 'DN'], inplace=True)

provinces = gpd.read_file('Data/GIS/Admin/Provinces.gpkg', encoding='utf-8')
demand_sites = gpd.read_file('../Morocco dashboard/spatial_data/Demand_points.gpkg', encoding='utf-8')

MerchidSudMoroc = 26192
for gdf in [provinces, demand_sites]:
    gdf.to_crs(epsg=MerchidSudMoroc, inplace=True)

for column in list(cropland)[1:-3]:
    cropland[column] = cropland[column].astype('float')
    cropland[column] = cropland[column].interpolate()
cropland['wtd'] = cropland['wtd'].astype('float')
    
for column in list(builtarea)[0:-2]:
    builtarea[column] = builtarea[column].astype('float')

cropland.crs = f'epsg:{MerchidSudMoroc}'
cropland = cropland.loc[cropland.area_m2>=100] #choose

dff = cropland.loc[cropland.province=='Taroudannt', 'wtd']
cropland.loc[cropland.province=='Taroudannt', 'wtd'] *= 170 * dff.count() / dff.sum()
dff = cropland.loc[cropland.province=='Chtouka-Aït Baha', 'wtd']
cropland.loc[cropland.province=='Chtouka-Aït Baha', 'wtd'] *= 85 * dff.count() / dff.sum()

def create_time_data(data, first_year, end_year):
    df = pd.DataFrame()
    df['Year'] = list(itertools.chain.from_iterable([[year]*12*data.shape[0] for year in range(first_year,end_year+1)]))
    y_number = df.Year.unique().shape[0]
    dff = pd.DataFrame(data.filter(regex='_srad|province|Demand point|area|wtd').melt(id_vars=['Demand point','province','area_m2','wtd']))
    dff.rename(columns={'variable': 'Month', 'value': 'srad'}, inplace=True) 
    dff = dff.join(data.filter(like='_wind').melt())
    dff.rename(columns={'value': 'wind'}, inplace=True)
    dff.drop(columns=['variable'], inplace=True)
    dff['Month'] = dff['Month'].str.replace('_srad','').astype(int)
    dff.sort_values(['province','Month'], inplace=True)
    dff.reset_index(drop=True, inplace=True)
    df = df.join(pd.concat([dff]*y_number, ignore_index=True))
    df['Date'] = pd.to_datetime(df[['Year','Month']].join(pd.DataFrame({'day': [1]*df.shape[0]})))
    return df
    
df_cropland = create_time_data(cropland, 2019, 2050)

def get_area_share(df, by, value):
    dff = df.groupby(by)[value].sum().copy()
    return df[value] / df.set_index(by).index.map(dff)
    
cropland_temp = cropland.copy()
cropland_temp.loc[cropland_temp['province']=='Inezgane-Aït Melloul', 'province'] = 'Taroudannt'
df_cropland['area_share'] = df_cropland['Demand point'].map(get_area_share(cropland_temp, 'province', 'area_m2'))

output_folder = r'Data/Cropland and Builtarea'
cropland.to_crs('epsg:4326').to_file(r"Data\GIS\Processed layers/cropland.geojson", driver='GeoJSON')
df_cropland.to_csv(os.path.join(output_folder, 'cropland.gz'), index=False)