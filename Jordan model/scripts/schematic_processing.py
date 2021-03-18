import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import nexus_tool.weap_tools as wp
import matplotlib.pyplot as plt
import shapely
import os

## Reading files

gis_folder = os.path.join('Data', 'GIS')

governorates = gpd.read_file(os.path.join(gis_folder, 'Admin', 'JOR_adm1.shp'))
links = gpd.read_file(str(snakemake.input.links))
groundwater = gpd.read_file(str(snakemake.input.groundwater))
river_withdrawals = gpd.read_file(str(snakemake.input.river_withdrawals))
wwtp = gpd.read_file(str(snakemake.input.wwtp))
other_supply = gpd.read_file(str(snakemake.input.other_supply))
demand_sites = gpd.read_file(str(snakemake.input.demand_sites))
tributary_inflows = gpd.read_file(str(snakemake.input.tributary_inflows))
diversion_outflows = gpd.read_file(str(snakemake.input.diversion_outflows))
diversion = gpd.read_file(str(snakemake.input.diversion))

output_demand_points = str(snakemake.output.demand_points)
output_folder = output_demand_points.split(os.path.basename(output_demand_points))[0]
output_points_coords = str(snakemake.output.points_coords)
dash_folder = output_points_coords.split(os.path.basename(output_points_coords))[0]

## Converting geometries and dropping unecesary columns

for gdf in [links, groundwater, river_withdrawals, wwtp, other_supply, 
            demand_sites, tributary_inflows, diversion, diversion_outflows]:
    gdf['geometry'] = wp.polyz_to_poly(gdf)
    gdf.drop(columns='Description', inplace=True)

## Defining supply, demand points and transmissions

for gdf in [groundwater, river_withdrawals, wwtp, other_supply, 
            diversion_outflows, demand_sites, tributary_inflows]:
    gdf.rename(columns={'Name': 'point'}, inplace=True)

links.rename(columns={'Name': 'links'}, inplace=True)
diversion.rename(columns={'Name': 'diversion'}, inplace=True)

## Reprojecting layers

PalestineBelt = 28192
for gdf in [links, groundwater, governorates, river_withdrawals, wwtp, 
            other_supply, demand_sites, tributary_inflows, diversion, 
            diversion_outflows]:
    gdf.crs = 4326
    gdf.to_crs(epsg=PalestineBelt, inplace=True)

## Defining types

groundwater['type'] = 'Groundwater supply'
river_withdrawals['type'] = 'River/pipeline supply'
wwtp['type'] = 'Wastewater plant'
other_supply['type'] = 'Other supply'
other_supply.loc[other_supply.point=='RedDead', 'type'] = 'Desalination'
other_supply.loc[other_supply.point=='Aqaba Desal', 'type'] = 'Desalination'
other_supply = other_supply.loc[other_supply['type']!='Other supply']
demand_sites['type'] = 'Municipality'
demand_sites.loc[demand_sites['point'].str.contains('Agri'), 'type'] = 'Agriculture'
demand_sites.loc[demand_sites['point'].str.contains('Ind'), 'type'] = 'Industry'
tributary_inflows['type'] = 'Tributary inflow'
diversion['type'] = 'Transmission Pipeline'
diversion_outflows['type'] = 'Diversion Outflow'

## Joining water transmition network points

points = tributary_inflows.append([diversion_outflows, river_withdrawals],
                                  ignore_index=True, sort=False)

## Calculate length of distribution and transmission links

links['length_m'] = links.length
diversion['pl_length_m'] = diversion.length

## Sampling raster data

groundwater['wtd_m'] = wp.sample_raster(os.path.join(gis_folder, 
                                                     'Water Table Depth', 
                                                     'Jordan_wtd_projected.tif'), 
                                        groundwater)

supply = groundwater.append([river_withdrawals, wwtp, other_supply], 
                            ignore_index=True, sort=False)
supply = gpd.sjoin(supply, links, how='inner', op='intersects')
supply.drop(columns='index_right', inplace=True)
demand = gpd.sjoin(demand_sites, links, how='inner', op='intersects')
demand.drop(columns='index_right', inplace=True)

points['elevation_m'] = wp.sample_raster(os.path.join(gis_folder, 
                                                      'DEM', 
                                                      'DEM_projected.tif'), 
                                         points)
demand['elevation_m'] = wp.sample_raster(os.path.join(gis_folder, 
                                                      'DEM', 
                                                      'DEM_projected.tif'), 
                                         demand)
supply['elevation_m'] = wp.sample_raster(os.path.join(gis_folder, 
                                                      'DEM', 
                                                      'DEM_projected.tif'), 
                                         supply)

## Intersect distribution links with supply and demand

distribution = gpd.sjoin(links, supply, how='inner', op='intersects')
distribution.drop(columns='index_right', inplace=True)
distribution = gpd.sjoin(distribution, demand, how='inner', op='intersects')
distribution.drop(columns='index_right', inplace=True)
distribution['type'] = 'Distribution link'

## Split segment in the pipeline system

points_coords = points.geometry.unary_union
list_gdf = []
i = 0

for pipeline in diversion.iterrows():
    split_pipeline = shapely.ops.split(pipeline[1].geometry, points_coords)
    segments = [feature for feature in split_pipeline]
    index = list(range(i, len(segments) + i))
    gdf_segments = gpd.GeoDataFrame(index, geometry=segments)
    gdf_segments.columns = ['index', 'geometry']
    gdf_segments['pipeline'] = pipeline[1].diversion
    gdf_segments['pl_length_m'] = pipeline[1].pl_length_m
    gdf_segments.crs = f'epsg:{PalestineBelt}'
    intersections = gpd.sjoin(gdf_segments, points, how='inner', op='intersects')
    list_gdf.append(intersections)
    i = list_gdf[-1]['index'].max() + 1
    
pipelines = gpd.GeoDataFrame(pd.concat(list_gdf, ignore_index=True))
pipelines.crs = supply.crs

## Calculate each segment lenght

pipelines['segment_length_m'] = pipelines.length

## Standardizing names

demand['point'] = demand['point'].str.replace('Agriculture', 'Agri')

## Save the spatial layers

os.makedirs(output_folder, exist_ok=True)
demand.to_file(os.path.join(output_folder, 'Demand_points.gpkg'), driver='GPKG')
supply.to_file(os.path.join(output_folder, 'Supply_points.gpkg'), driver='GPKG')
pipelines.to_file(os.path.join(output_folder, 'Pipelines.gpkg'), driver='GPKG')

## Save coordinates for visualization

all_points_coords = pd.DataFrame(demand.append(supply, sort=False, ignore_index=True).to_crs(epsg=4326))
all_points_coords.drop_duplicates(subset="point", inplace=True)
all_points_coords['lon'] = [point.xy[0][0] for point in all_points_coords.geometry]
all_points_coords['lat'] = [point.xy[1][0] for point in all_points_coords.geometry]
all_points_coords.drop(columns=['geometry', 'links', 'length_m', 'elevation_m', 'wtd_m'], inplace=True)

pipe_coords = pd.DataFrame({'lon': [], 'lat': []})
for name, point in zip(pipelines.pipeline, pipelines.to_crs(epsg=4326).geometry):
    lon = list(point.xy[0]) + [None]
    lat = list(point.xy[1]) + [None]
    df_temp = pd.DataFrame({'lon': lon, 'lat': lat})
    df_temp['name'] = name
    pipe_coords = pipe_coords.append(df_temp, ignore_index=True)
pipe_coords['type'] = 'pipeline'

os.makedirs(dash_folder, exist_ok=True)
all_points_coords.to_csv(output_points_coords, index=False)
pipe_coords.to_csv(os.path.join(dash_folder, 'pipe_coords.csv'), index=False)