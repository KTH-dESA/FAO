import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
import nexus_tool.weap_tools as wp
import re
from functools import reduce
import os
from shutil import copyfile

#Read processed schematic files
sp_folder = os.path.join('..', 'Jordan dashboard', 'spatial_data')
demand = gpd.read_file(os.path.join(sp_folder, 'Demand_points.gpkg'))
supply = gpd.read_file(os.path.join(sp_folder, 'Supply_points.gpkg'))
pipelines = gpd.read_file(os.path.join(sp_folder, 'Pipelines.gpkg'))