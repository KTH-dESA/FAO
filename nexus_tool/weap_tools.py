#Standard library imports
import os
import glob
import geopandas as gpd
import pandas as pd
import numpy as np
import fiona
import shapely
import rasterio
import rasterio.mask
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.fill import fillnodata
from rasterstats import zonal_stats
import scipy
fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support
fiona.drvsupport.supported_drivers['KML'] = 'rw' # enable KML support

def polyz_to_poly(gdf):
    return gdf.geometry.map(lambda poly: 
                        shapely.ops.transform(lambda x, y, z: (x, y), poly))
                        
def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)
        
def mask_raster(raster_path, mask_path, crs):
    with fiona.open(mask_path, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    with rasterio.open(raster_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_image[out_image<0] = np.nan
        mask = (out_image!=0)
        out_image = fillnodata(out_image, mask)
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform,
                     "crs": crs})
    return out_image, out_meta
    
def reproject_raster(raster_path, dst_crs, outpul_file):
    with rasterio.open(raster_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(outpul_file, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
                    
def sample_raster(path, gdf):
    with rasterio.open(path) as src:
        return [float(val) for val in src.sample([(x.coords.xy[0][0], 
                                                   x.coords.xy[1][0]) for x in 
                                                   gdf['geometry']])]
                                                   
def merge_rasters(files_path, dst_crs, outpul_file):
    files = glob.glob(files_path)
    src_files_to_mosaic = []
    
    for fp in files:
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)
    
    mosaic, out_trans = merge(src_files_to_mosaic)

    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff",
                 "height": mosaic.shape[1],
                 "width": mosaic.shape[2],
                 "transform": out_trans,
                 "crs": dst_crs
                 }
                )
    
    with rasterio.open(outpul_file, "w", **out_meta) as dest:
        dest.write(mosaic)
        
def get_zonal_stats(vector, path, stats):
    # Run zonal statistics, store result in geopandas dataframe
    # with rasterio.open(path) as src:
    result = zonal_stats(vector, path, stats=stats, geojson_out=True)
    geostats = gpd.GeoDataFrame.from_features(result)
    return geostats

def create_learning_curve(iyear, eyear, cc, rates, method, order=2):
    learning_curve = pd.DataFrame({'Year': range(iyear, eyear+1)})
    if method!='step':
        for year, rate in rates.items():
            learning_curve.loc[learning_curve.Year==year, 'rate'] = rate
        learning_curve['capital_cost'] = (1-learning_curve['rate']) * cc
        learning_curve['capital_cost'] = learning_curve.capital_cost.interpolate(method=method, order=order)
    else:
        for year, rate in rates.items():
            learning_curve.loc[learning_curve.Year>=year, 'rate'] = rate
        learning_curve['capital_cost'] = (1-learning_curve['rate']) * cc
    return learning_curve.set_index('Year').capital_cost