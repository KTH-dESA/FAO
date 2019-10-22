#Standard library imports
import pandas as pd
import numpy as np

#### default values:
gw_depth = 'gw_depth'
tdh_gw = 'tdh_gw'

def get_gw_tdh(df, gw_depth = gw_depth, wdd = 0, oap = 0, pld = 0, 
               interp_method = 'nearest', tdh_gw = tdh_gw):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw] = df[tdh_gw].replace(0,np.nan)
    df = df.interpolate(method = interp_method, axis=0)
    return df