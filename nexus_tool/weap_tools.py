#Standard library imports
import os
import geopandas as gpd
import fiona
import shapely
import rasterio
import rasterio.mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.fill import fillnodata
fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support
fiona.drvsupport.supported_drivers['KML'] = 'rw' # enable KML support

def polyz_to_poly(gdf):
    return gdf.geometry.map(lambda poly: 
                        shapely.ops.transform(lambda x, y, z: (x, y), poly))
                        
def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)
    