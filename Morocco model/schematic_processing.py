import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import os
import geopandas as gpd
import nexus_tool.weap_tools as wp

provinces = gpd.read_file('Data/GIS/Admin/Provinces.gpkg', encoding='utf-8')
groundwater = gpd.read_file('Data/Schematic/Groundwater.kml', encoding='utf-8')
wwtp = gpd.read_file('Data/Schematic/Wastewater Treatment Plants.kml', encoding='utf-8')
other_supply = gpd.read_file('Data/Schematic/Other Supplies.kml', encoding='utf-8')
river_withdrawals = gpd.read_file('Data/Schematic/River Withdrawals.kml', encoding='utf-8')
demand_sites = gpd.read_file('Data/Schematic/Demand Sites.kml', encoding='utf-8')
catchments = gpd.read_file('Data/Schematic/Catchments.kml', encoding='utf-8')
diversion = gpd.read_file('Data/Schematic/Diversions.kml', encoding='utf-8')
reservoirs = gpd.read_file('Data/Schematic/Reservoirs.kml', encoding='utf-8')
links = gpd.read_file('Data/Schematic/Transmission Links.kml', encoding='utf-8')

#Converting geometries and droping unecesary columns

for gdf in [groundwater, wwtp, other_supply, river_withdrawals,
            demand_sites, catchments, diversion, reservoirs, links]:
    gdf['geometry'] = wp.polyz_to_poly(gdf)
    gdf.drop(columns='Description', inplace=True)
    
for gdf in [groundwater, wwtp, other_supply, river_withdrawals,
            reservoirs, demand_sites, catchments]:
    gdf.rename(columns={'Name': 'point'}, inplace=True)

links.rename(columns={'Name': 'links'}, inplace=True)
diversion.rename(columns={'Name': 'diversion'}, inplace=True)

#Defining types

groundwater['type'] = 'Groundwater supply'
wwtp['type'] = 'WWTP'
other_supply['type'] = 'Other supply'
river_withdrawals['type'] = 'Surfacewater withdrawal'
reservoirs['type'] = 'Reservoir supply'
catchments['type'] = 'Catchment'
demand_sites['type'] = 'Demand site'
diversion['type'] = 'Transmission Pipeline'
links['type'] = 'Transmission links'

supply_points = groundwater.append([other_supply, reservoirs, river_withdrawals], ignore_index=True, sort=False)
demand_points = demand_sites.append(catchments, ignore_index=True, sort=False)
links['links'] = links.links.str.replace('Groundwater','GW')
links['links'] = links.links.str.replace('Grounwater','GW')
links['links'] = links.links.str.replace('GW of ','')
links['links'] = links.links.str.replace('GW ','')
links['links'] = links.links.str.strip()
links['links'] = links.links.str.replace(' El Guerdane',' I_El Guerdane')

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
    
groundwater['wtd_m'] = wp.sample_raster("Data/GIS/wtd/Africa_model_wtd_v2.nc", 
                                         groundwater)
groundwater.loc[groundwater.point=='Souss GW','wtd_m'] = 170
groundwater.loc[groundwater.point=='Chtouka GW','wtd_m'] = 85

if not os.path.exists("Data/GIS/DEM/Souss-Massa DEM.tif"):
    wp.merge_rasters('Data/GIS/DEM/*', 'EPSG:4326', 
                     "Data/GIS/DEM/Souss-Massa DEM.tif")
    
supply_links['elevation'] = wp.sample_raster("Data/GIS/DEM/Souss-Massa DEM.tif", 
                                               supply_links)
demand_links['elevation'] = wp.sample_raster("Data/GIS/DEM/Souss-Massa DEM.tif", 
                                               demand_links)
wwtp['elevation'] = wp.sample_raster("Data/GIS/DEM/Souss-Massa DEM.tif", wwtp)
                                       
diversions = gpd.sjoin(diversion, reservoirs, how='inner', op='intersects')
diversions.drop(columns=['index_right'], inplace=True)
diversions.rename(columns={'type_right': 'type_supply'}, inplace=True)
diversions = gpd.sjoin(diversions, river_withdrawals, how='inner', op='intersects')
diversions.drop(columns=['index_right'], inplace=True)
diversions.rename(columns={'point_left': 'Supply point', 'point_right': 'Demand point',
                           'type_left': 'type', 'type': 'type_demand'}, inplace=True)                           

dff1 = river_withdrawals.loc[river_withdrawals.point.isin(diversions['Demand point'])].copy()
dff1['diversion'] = dff1.point.map(diversions.set_index('Demand point').diversion)
dff1['elevation'] = wp.sample_raster("Data/GIS/DEM//Souss-Massa DEM.tif", dff1)

dff2 = reservoirs.loc[reservoirs.point.isin(diversions['Supply point'])].copy()
dff2['diversion'] = dff2.point.map(diversions.groupby('Supply point').agg({'diversion': 'first'})['diversion'])
dff2['elevation'] = wp.sample_raster("Data/GIS/DEM//Souss-Massa DEM.tif", 
                                       dff2)
MerchidSudMoroc = 26192
dff2.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)

dff1.to_crs(f'epsg:{MerchidSudMoroc}', inplace=True)
dff1.set_index('diversion', inplace=True)
dff1['distance'] = dff1.distance(dff2.set_index('diversion')).copy()

diversions['distance'] = diversions['Demand point'].map(dff1.set_index('point')['distance'])
diversions['elevation_diff'] = diversions['Demand point'].map(dff1.set_index('point')['elevation']) - \
                               diversions['Supply point'].map(dff2.set_index('point')['elevation'])

all_points = supply_links.append(demand_links, sort=False, ignore_index=True)
all_points.drop_duplicates(subset="point", inplace=True)

folder = r'../Morocco dashboard/spatial_data'
wp.create_folder(folder)
demand_links.to_file(os.path.join(folder, 'Demand_points.gpkg'), driver="GPKG")
supply_links.to_file(os.path.join(folder, 'Supply_points.gpkg'), driver="GPKG")
wwtp.to_file(os.path.join(folder, 'wwtp.gpkg'), driver="GPKG")
diversions.to_file(os.path.join(folder, 'Pipelines.gpkg'), driver="GPKG")
groundwater.to_file(os.path.join(folder, 'Groundwater.gpkg'), driver="GPKG")
all_points.to_file(os.path.join(folder, 'all_points.gpkg'), driver="GPKG")