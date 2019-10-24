#Standard library imports
import pandas as pd
import numpy as np

def get_gw_tdh(df, gw_depth, wdd, oap, pld, 
               interp_method = 'nearest', tdh_gw = tdh_gw):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw].replace(0, np.nan, inplace=True)
    # df[tdh_gw].interpolate(method = interp_method, axis=0, inplace=True)
    return df
    
def get_pumping_energy(df, trans_eff, pump_eff, pd_e, pwd, sswd, ed_e, tdh_gw, 
                       des_int, des_ener, desalination = False):
    Epump_plant_eff = trans_eff * pump_eff

    for i in range (1,13):
        _pd_e = '{}{}'.format(pd_e, i)
        _pwd = '{}{}'.format(pwd, i)
        _sswd = '{}{}'.format(sswd, i)
        _ed_e = '{}{}'.format(ed_e, i)
        
        df[_pd_e]=(9.81*(df[_pwd]/1000)*df[tdh_gw])/Epump_plant_eff
        df[_ed_e]=(df[_sswd]*df[tdh_gw]*0.00272)/Epump_plant_eff
        
        if desalination:
            df[_pd_e] += (df[_pwd]*df[des_int]*3600/1000)
            df[_ed_e] += (df['{}{}'.format(des_ener, i)]*1000000)
    return df