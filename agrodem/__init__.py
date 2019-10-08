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
    evap_i,
    eto,
    effective_rainfall,
    kc,
    get_harvest_fraction,
)

math.exp = np.exp
math.pow = np.power
math.sqrt = np.sqrt

class DataFrame(pd.DataFrame):
    property_1 = None
    prefix = 'ETo_'
    lat = 'lat'
    elevation = 'elevation'
    wind = 'wind'
    srad = 'srad'
    tmin = 'tmin'
    tmax = 'tmax',
    tavg = 'tavg'
    def test(self, var, inplace = False):
        if inplace:
            self.property_1 = var
        else: 
            return test(var)
          
    def eto(self, inplace = False):
        if inplace:
            eto(self, self.prefix, self.lat, self.elevation, self.wind, 
                self.srad, self.tmin, self.tmax, self.tavg)
        else:
            return eto(self, self.prefix, self.lat, self.elevation, self.wind, 
                       self.srad, self.tmin, self.tmax, self.tavg)
                
    evap_i = staticmethod(evap_i)
    effective_rainfall = staticmethod(effective_rainfall)
    kc = staticmethod(kc)
    get_harvest_fraction = staticmethod(get_harvest_fraction)