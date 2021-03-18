import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import rasterio
from rasterio.plot import show
import nexus_tool.weap_tools as wp
import matplotlib.pyplot as plt
import os


country_border = str(snakemake.input.country_border)
wtd = str(snakemake.input.wtd)
dem = str(snakemake.input.dem)

wtd_masked = str(snakemake.output.wtd_masked )
dem_masked  = str(snakemake.output.dem_masked )
wtd_mask_projected = str(snakemake.output.wtd_mask_projected)
dem_projected = str(snakemake.output.dem_projected)
dem_mask_projected = str(snakemake.output.dem_mask_projected)

## Mask rasters into country boundaries

#Masking water table depth from the Euroasian model v2 to Jordan boundaries
out_image, out_meta = wp.mask_raster(raster_path=wtd, 
                                     mask_path=country_border, 
                                     crs='EPSG:4326')

#Writing masked raster to file
with rasterio.open(wtd_masked, 
                   "w", **out_meta) as dest:
    dest.write(out_image)

#Masking elevation from DEM model to Jordan boundaries
out_image, out_meta = wp.mask_raster(dem,
                                     country_border, 
                                     'EPSG:4326')

#Writing masked raster to file
with rasterio.open(dem_masked, 
                   "w", **out_meta) as dest:
    dest.write(out_image)

## Reprojecting rasters

dst_crs = 'EPSG:28192' #Define crs system (28192=PalestineBelt)

wp.reproject_raster(wtd_masked, 
                    dst_crs, 
                    wtd_mask_projected)
wp.reproject_raster(dem,
                    dst_crs, 
                    dem_projected)
wp.reproject_raster(dem_masked,
                    dst_crs, 
                    dem_mask_projected)