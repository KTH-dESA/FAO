#Standard library imports
import pandas as pd
from pandas import read_csv
import numpy as np
import math

#Related third party imports
import pyeto

#Local application/library specific imports
from agrodem.water_demand import (
    test,
    set_cropland_share,
    get_ky_list,
    get_kc_list,
    get_evap_i,
    get_eto,
    get_effective_rainfall,
    get_kc_i,
    get_kc_values,
    get_harvest_fraction,
)

math.exp = np.exp
math.pow = np.power
math.sqrt = np.sqrt


class Model():
    property_1 = None
    eto = 'ETo_'
    lat = 'lat'
    elevation = 'elevation'
    wind = 'wind'
    srad = 'srad'
    tmin = 'tmin'
    tmax = 'tmax'
    tavg = 'tavg'
    crop_share = 'CropShare'
    pumping_hours_per_day = 10
    deff = 1
    aeff = 0.65
    ky_values = {}
    kc_values = {}
    def __init__(self, df, eto = eto, lat = lat, elevation = elevation,
                 wind = wind, srad = srad, tmin = tmin, tmax = tmax, 
                 tavg = tavg, crop_share = crop_share,
                 pumping_hours_per_day = pumping_hours_per_day,
                 deff = deff, aeff = aeff):
        self.df = df
        self.eto = eto
        self.lat = lat
        self.elevation = elevation
        self.wind = wind
        self.srad = srad
        self.tmin = tmin
        self.tmax = tmax
        self.tavg = tavg
        self.crop_share = crop_share
        self.pumping_hours_per_day = 10
        self.deff = 1
        self.aeff = 0.65
    
    def test(self, var, inplace = False):
        if inplace:
            self.property_1 = var
        else: 
            return test(var)
            
    def print_properties(self):
        print('Properties names:')
        for val, name in zip([self.eto, self.lat, self.elevation, self.wind, 
                              self.srad, self.tmin, self.tmax, self.tavg],
                             ['Reference evapotranspiration (.eto)', 
                              'Latitude (.lat)', 'Elevation (.elevation)', 
                              'Wind speed (.wind)', 'Solar radiation (.srad)', 
                              'Min temperature (.tmin)', 'Max temperature (.tmax)', 
                              'Avegarage temperature (.tavg)']):
            print('    - {}: {}'.format(name, val))
          
    def set_cropland_share(self, crop_var, geo_boundary = 'global', 
                           boundary_name = None, inplace = False):
        if inplace:
            set_cropland_share(self.df, crop_var, geo_boundary = geo_boundary, 
                       boundary_name = boundary_name, crop_share = self.crop_share)
        else:
            return set_cropland_share(self.df, crop_var, 
                                      geo_boundary = geo_boundary, 
                                      boundary_name = boundary_name, 
                                      crop_share = self.crop_share)
                                      
    def get_ky_list(self, inplace = False):
        if inplace:
            get_ky_list(self.df, crop_share = self.crop_share)
        else:
            return get_ky_list(self.df, crop_share = self.crop_share)
           
    def get_kc_list(self, inplace = False):
        if inplace:
            get_kc_list(self.df, crop_share = self.crop_share)
        else:
            return get_kc_list(self.df, crop_share = self.crop_share)
    
    def get_eto(self, inplace = False):
        if inplace:
            get_eto(self.df, self.eto, self.lat, self.elevation, self.wind, 
                    self.srad, self.tmin, self.tmax, self.tavg)
        else:
            return get_eto(self.df, self.eto, self.lat, self.elevation, self.wind, 
                           self.srad, self.tmin, self.tmax, self.tavg)
                
    get_evap_i = staticmethod(get_evap_i)
    get_effective_rainfall = staticmethod(get_effective_rainfall)
    get_kc_i = staticmethod(get_kc_i)
    get_harvest_fraction = staticmethod(get_harvest_fraction)