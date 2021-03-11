import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import os
import geopandas as gpd
import nexus_tool.weap_tools as wp

#Reading files
gis_folder = os.path.join('Data', 'GIS')
sch_folder = os.path.join('Data', 'Schematic')

governorates = gpd.read_file(os.path.join(gis_folder, 'Admin', 'JOR_adm1.shp'))
links = gpd.read_file(os.path.join(sch_folder, 'Transmission Links.kml'))
groundwater = gpd.read_file(os.path.join(sch_folder, 'Groundwater.kml'))
river_withdrawals = gpd.read_file(os.path.join(sch_folder, 'River Withdrawals.kml'))
wwtp = gpd.read_file(os.path.join(sch_folder, 'Wastewater Treatment Plants.kml'))
other_supply = gpd.read_file(os.path.join(sch_folder, 'Other Supplies.kml'))
demand_sites = gpd.read_file(os.path.join(sch_folder, 'Demand Sites.kml'))
tributary_inflows = gpd.read_file(os.path.join(sch_folder, 'Tributary Inflows.kml'))
diversion_outflows = gpd.read_file(os.path.join(sch_folder, 'Diversion Outflows.kml'))
diversion = gpd.read_file(os.path.join(sch_folder, 'Diversions.kml'))

#Converting geometries and dropping unecesary columns
for gdf in [links, groundwater, river_withdrawals, wwtp, other_supply, 
            demand_sites, tributary_inflows, diversion, diversion_outflows]:
    gdf['geometry'] = wp.polyz_to_poly(gdf)
    gdf.drop(columns='Description', inplace=True)
    
#Defining supply, demand points and transmissions
for gdf in [groundwater, river_withdrawals, wwtp, other_supply, 
            diversion_outflows, demand_sites, tributary_inflows]:
    gdf.rename(columns={'Name': 'point'}, inplace=True)

links.rename(columns={'Name': 'links'}, inplace=True)
diversion.rename(columns={'Name': 'diversion'}, inplace=True)

#Reprojecting layers
PalestineBelt = 28192
for gdf in [links, groundwater, governorates, river_withdrawals, wwtp, 
            other_supply, demand_sites, tributary_inflows, diversion, 
            diversion_outflows]:
    gdf.crs = 4326
    gdf.to_crs(epsg=PalestineBelt, inplace=True)
    
#Defining types
groundwater['type'] = 'Groundwater supply'
river_withdrawals['type'] = 'River/pipeline supply'
wwtp['type'] = 'WWTP'
other_supply['type'] = 'Other supply'
demand_sites['type'] = 'Demand site'
tributary_inflows['type'] = 'Tributary inflow'
diversion['type'] = 'Transmission Pipeline'
diversion_outflows['type'] = 'Diversion Outflow'

#Joining water transmition network points
points = tributary_inflows.append([diversion_outflows, river_withdrawals],
                                  ignore_index=True, sort=False)
                                  
#Calculate length of distribution and transmission links
links['length_m'] = links.length
diversion['pl_length_m'] = diversion.length

#Sampling raster data
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

#Define desalination plants
supply.loc[supply.point=='RedDead', 'type'] = 'Desalination'
supply.loc[supply.point=='Aqaba Desal', 'type'] = 'Desalination'
                  
#Intersect distribution links with supply and demand
distribution = gpd.sjoin(links, supply, how='inner', op='intersects')
distribution.drop(columns='index_right', inplace=True)
distribution = gpd.sjoin(distribution, demand, how='inner', op='intersects')
distribution.drop(columns='index_right', inplace=True)
distribution['type'] = 'Distribution link'    

#Split segment in the pipeline system
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

#Calculate each segment lenght
pipelines = gpd.GeoDataFrame(pd.concat(list_gdf, ignore_index=True))
pipelines['segment_length_m'] = pipelines.length

#Save the spatial layers
#Define the output folder path
folder = os.path.join('..', 'Jordan dashboard', 'spatial_data')
demand.to_file(os.path.join(folder, 'Demand_points.gpkg'), driver='GPKG')
supply.to_file(os.path.join(folder, 'Supply_points.gpkg'), driver='GPKG')
pipelines.to_file(os.path.join(folder, 'Pipelines.gpkg'), driver='GPKG')
              