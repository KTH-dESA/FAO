import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
import os
from nexustool.gis_tools import download_data, create_time_data, get_area_share, get_zonal_stats
from nexustool.weap_tools import reproject_raster, sample_raster

## Downloading solar irradiation and water table depth data

url = 'https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_srad.zip'
file_path = os.path.join('data', 'gis', 'srad', 'wc2.1_30s_srad.zip')
download_data(url, file_path)

url = 'https://souss-massa-dev.s3.us-east-2.amazonaws.com/post_build/Africa_model_wtd_v2.nc'
file_path = os.path.join('data', 'gis', 'wtd', 'Africa_model_wtd_v2.nc')
download_data(url, file_path)

## Reading the input data

demand_path = str(snakemake.input.demand_points)
cropland_path = str(snakemake.input.cropland)

crop_df = pd.read_csv(cropland_path, encoding='utf-8')
geometry = crop_df['WKT'].map(shapely.wkt.loads)
cropland = gpd.GeoDataFrame(crop_df.drop(columns=['WKT']), crs="EPSG:26192", geometry=geometry)
provinces = gpd.read_file(os.path.join('data', 'gis', 'admin', 'provinces.gpkg'), encoding='utf-8')

output_file = str(snakemake.output)
output_folder = output_file.split(os.path.basename(output_file))[0]

## Convert coordenate reference system (crs)

MerchidSudMoroc = 26192
for gdf in [provinces, provinces]:
    gdf.to_crs(epsg=MerchidSudMoroc, inplace=True)

cropland = cropland.loc[cropland.area_m2>=100] #choose

## Solar irradiation zonal statistics
Loops through the 12 months of the year and gets the mean solar irradiation of each month within each cropland polygon

cropland.to_crs(epsg=4326, inplace=True)
for month in range(1, 13):
    cropland = get_zonal_stats(cropland, 
                               os.path.join('data', 'gis', 'srad', 
                                            f'wc2.1_30s_srad_{str(month).zfill(2)}.tif'), 
                               ['mean'], all_touched=True).rename(columns={'mean': f'srad{month}'})

## Water table depth zonal statistics

cropland.crs = 4326
cropland = get_zonal_stats(cropland, 
                           os.path.join('data', 'gis', 'wtd', 
                                        'Africa_model_wtd_v2.nc'), 
                           ['mean'], all_touched=True).rename(columns={'mean': 'wtd'})

cropland.crs = 4326
cropland.to_crs(epsg=MerchidSudMoroc, inplace=True)

## Creating time series data

df_cropland = create_time_data(cropland, 2019, 2050)


## Calculating the area share of each croplan area within each province

cropland.loc[cropland['province']=='Inezgane-Aït Melloul', 'province'] = 'Taroudannt' #Including Inezgane-Aït Melloul irrigated area into results from Taroudant due to lack of data for the former
cropland['area_share'] = get_area_share(cropland, 'province', 'area_m2')

df_cropland = pd.merge(df_cropland, cropland[['Demand point', 'area_share']], on='Demand point')

os.makedirs(output_folder, exist_ok = True)
df_cropland.to_csv(output_file, index=False)