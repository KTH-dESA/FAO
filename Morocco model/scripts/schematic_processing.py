import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import os
import geopandas as gpd
import rasterio
import nexustool.weap_tools as wp
from nexustool.gis_tools import download_data
from softlinking_functions import integrate_data, data_merging
import fiona
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw' # enable KML support which is disabled by default

provinces = gpd.read_file(str(snakemake.input.provinces), encoding='utf-8')
groundwater = gpd.read_file(str(snakemake.input.groundwater), encoding='utf-8')
wwtp = gpd.read_file(str(snakemake.input.wwtp), encoding='utf-8')
other_supply = gpd.read_file(str(snakemake.input.other_supply), encoding='utf-8')
river_withdrawals = gpd.read_file(str(snakemake.input.river_with), encoding='utf-8')
demand_sites = gpd.read_file(str(snakemake.input.demand_sites), encoding='utf-8')
catchments = gpd.read_file(str(snakemake.input.catchments), encoding='utf-8')
diversion = gpd.read_file(str(snakemake.input.diversions), encoding='utf-8')
reservoirs = gpd.read_file(str(snakemake.input.reservoirs), encoding='utf-8')
links = gpd.read_file(str(snakemake.input.links), encoding='utf-8')

gis_folder = os.path.join('data', 'gis')
output_demand_points = str(snakemake.output.demand_points)
output_folder = output_demand_points.split(os.path.basename(output_demand_points))[0]
spatial_output = str(snakemake.output.spatial_output)

## Downloading the GIS data

url = 'http://viewfinderpanoramas.org/dem3/I29.zip'
file_path = os.path.join('data', 'gis', 'dem', 'I29.zip')
download_data(url, file_path)

url = 'http://viewfinderpanoramas.org/dem3/H29.zip'
file_path = os.path.join('data', 'gis', 'dem', 'H29.zip')
download_data(url, file_path)

url = 'https://souss-massa-dev.s3.us-east-2.amazonaws.com/post_build/Africa_model_wtd_v2.nc'
file_path = os.path.join('data', 'gis', 'wtd', 'Africa_model_wtd_v2.nc')
download_data(url, file_path)

## Converting geometries and dropping unecesary columns

for gdf in [groundwater, wwtp, other_supply, river_withdrawals,
            demand_sites, catchments, diversion, reservoirs, links]:
    gdf['geometry'] = wp.polyz_to_poly(gdf)
    gdf.drop(gdf.columns.difference(['Name','geometry']), 1, inplace=True)

for gdf in [groundwater, wwtp, other_supply, river_withdrawals,
            reservoirs, demand_sites, catchments]:
    gdf.rename(columns={'Name': 'point'}, inplace=True)

links.rename(columns={'Name': 'links'}, inplace=True)
diversion.rename(columns={'Name': 'diversion'}, inplace=True)

## Defining types

groundwater['type'] = 'Groundwater supply'
wwtp['type'] = 'WWTP'
other_supply['type'] = 'Other supply'
river_withdrawals['type'] = 'Surfacewater withdrawal'
reservoirs['type'] = 'Reservoir supply'
catchments['type'] = 'Catchment'
demand_sites['type'] = 'Demand site'
diversion['type'] = 'Transmission Pipeline'
links['type'] = 'Transmission links'

## Cleaning up the data

supply_points = groundwater.append([other_supply, reservoirs, river_withdrawals, wwtp], ignore_index=True, sort=False)
demand_points = demand_sites.append(catchments, ignore_index=True, sort=False)
links['links'] = links.links.str.replace('Groundwater','GW')
links['links'] = links.links.str.replace('Grounwater','GW')
links['links'] = links.links.str.replace('GW of ','')
links['links'] = links.links.str.replace('GW ','')
links['links'] = links.links.str.strip()
links['links'] = links.links.str.replace(' El Guerdane',' I_El Guerdane')

## Spatially joining some layers

demand_links = gpd.sjoin(demand_points, links, how='inner', op='intersects')
supply_links = gpd.sjoin(supply_points, links, how='inner', op='intersects')

demand_links.rename(columns={'type_left': 'type'}, inplace=True)
demand_links.drop(columns=['type_right', 'index_right'], inplace=True)
supply_links.rename(columns={'type_left': 'type'}, inplace=True)
supply_links.drop(columns=['type_right', 'index_right'], inplace=True)

dff = gpd.sjoin(demand_links, provinces[['Province','geometry']], how='inner', op='within')
dff.drop(columns=['index_right'], inplace=True)
df_temp = demand_links.loc[~demand_links.point.isin(dff.point.unique())].copy()
df_temp['Province'] = 'Tiznit'
demand_links = dff.append(df_temp, sort=False)

dff = gpd.sjoin(supply_links, provinces[['Province','geometry']], how='inner', op='within')
dff.drop(columns=['index_right'], inplace=True)
df_temp = supply_links.loc[~supply_links.point.isin(dff.point.unique())].copy()
supply_links = dff.append(df_temp, sort=False)
df_temp['Province'] = 'Chtouka-Aït Baha'
supply_links = dff.append(df_temp, sort=False)

## Displaying the system

base = provinces.plot(color='white', edgecolor='black', figsize=(12, 12))
data = demand_links.append([supply_links, wwtp, diversion, links], ignore_index=True, sort=False)
data.plot(ax=base, column='type', cmap='Spectral_r', legend=True)

## Masking and displaying rasters

out_image, out_meta = wp.mask_raster(os.path.join(gis_folder, 'wtd', 'Africa_model_wtd_v2.nc'), 
                                  os.path.join(gis_folder, 'admin', 'provinces.gpkg'), 'EPSG:4326')
       
with rasterio.open(os.path.join(gis_folder, 'wtd', 'Souss-Massa WTD.tif'), 'w', **out_meta) as dest:
    dest.write(out_image)

dem_path = os.path.join(gis_folder, 'dem', 'Souss-Massa DEM.tif')
if not os.path.exists(dem_path):
    wp.merge_rasters(os.path.join(gis_folder, 'dem', '*', '*.hgt'), 
                     'EPSG:4326', 
                     dem_path)

## Sampling raster data

groundwater['wtd_m'] = wp.sample_raster(os.path.join(gis_folder, 'wtd', 'Africa_model_wtd_v2.nc'), 
                                        groundwater)
groundwater.loc[groundwater.point=='Souss GW','wtd_m'] = 170
groundwater.loc[groundwater.point=='Chtouka GW','wtd_m'] = 85

supply_links['elevation'] = wp.sample_raster(dem_path, supply_links)
demand_links['elevation'] = wp.sample_raster(dem_path, demand_links)
wwtp['elevation'] = wp.sample_raster(dem_path, wwtp)

## Spatially joining diversions with reservoirs and river withdrawals

diversions = gpd.sjoin(diversion, reservoirs, how='inner', op='intersects')
diversions.drop(columns=['index_right'], inplace=True)
diversions.rename(columns={'type_right': 'type_supply'}, inplace=True)
diversions = gpd.sjoin(diversions, river_withdrawals, how='inner', op='intersects')
diversions.drop(columns=['index_right'], inplace=True)
diversions.rename(columns={'point_left': 'Supply point', 'point_right': 'Demand point',
                           'type_left': 'type', 'type': 'type_demand'}, inplace=True)

## Sampling elevation data

dff1 = river_withdrawals.loc[river_withdrawals.point.isin(diversions['Demand point'])].copy()
dff1['diversion'] = dff1.point.map(diversions.set_index('Demand point').diversion)
dff1['elevation'] = wp.sample_raster(dem_path, dff1)

dff2 = reservoirs.loc[reservoirs.point.isin(diversions['Supply point'])].copy()
dff2['diversion'] = dff2.point.map(diversions.groupby('Supply point').agg({'diversion': 'first'})['diversion'])
dff2['elevation'] = wp.sample_raster(dem_path, dff2)
MerchidSudMoroc = 26192
dff2.crs = 4326
dff2.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)

## Calculating pipelines distance

dff1.crs = 4326
dff1.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)
dff1.set_index('diversion', inplace=True)
dff1['distance'] = dff1.distance(dff2.set_index('diversion')).copy()

## Calculating elevation difference

diversions['distance'] = diversions['Demand point'].map(dff1.set_index('point')['distance'])
diversions['elevation_diff'] = diversions['Demand point'].map(dff1.set_index('point')['elevation']) - \
                               diversions['Supply point'].map(dff2.set_index('point')['elevation'])

all_points = supply_links.append(demand_links, sort=False, ignore_index=True)
all_points.drop_duplicates(subset="point", inplace=True)

os.makedirs(output_folder, exist_ok = True)

demand_links.to_file(output_demand_points, driver="GPKG", OVERWRITE='YES')
supply_links.to_file(os.path.join(output_folder, 'supply_points.gpkg'), 
                     driver="GPKG", OVERWRITE='YES')
wwtp.to_file(os.path.join(output_folder, 'wwtp.gpkg'), 
             driver="GPKG", OVERWRITE='YES')
diversions.to_file(os.path.join(output_folder, 'pipelines.gpkg'), 
                   driver="GPKG", OVERWRITE='YES')
groundwater.to_file(os.path.join(output_folder, 'groundwater.gpkg'), 
                    driver="GPKG", OVERWRITE='YES')
all_points.to_file(os.path.join(output_folder, 'all_points.gpkg'), 
                   driver="GPKG", OVERWRITE='YES')
                   

points_coords, pipe_coords = data_merging(demand_links, supply_links, diversions)

spatial_output_folder = spatial_output.split(os.path.basename(spatial_output))[0]
os.makedirs(spatial_output_folder, exist_ok=True)

points_coords.to_csv(os.path.join(spatial_output_folder, 'points_coords.csv'), index=False)
pipe_coords.to_csv(os.path.join(spatial_output_folder, 'pipe_coords.csv'), index=False)