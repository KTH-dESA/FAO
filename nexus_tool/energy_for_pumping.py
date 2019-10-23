#Standard library imports
import pandas as pd
import numpy as np

#### default values:
gw_depth = 'gw_depth'
tdh_gw = 'tdh_gw'
des_int = 'Einten_KWh/m3'
des_ener = 'Edesal_GWh_'
acwr = 'ACWR_' 
pcwr = 'PCWR_'
pwd = 'PWD_'
sswd = 'SSWD_'
pd_e = 'PD_E_'
ed_e = 'ED_E_'

def get_gw_tdh(df, gw_depth = gw_depth, wdd = 0, oap = 0, pld = 0, 
               interp_method = 'nearest', tdh_gw = tdh_gw):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw] = df[tdh_gw].replace(0,np.nan)
    df = df.interpolate(method = interp_method, axis=0)
    return df
    
def get_pumping_energy(df, trans_eff, pump_eff, pd_e = pd_e, pwd = pwd, 
                       sswd = sswd, ed_e = ed_e, tdh_gw = tdh_gw, 
                       desalination = False, des_int = des_int,
                       des_ener = des_ener):
    Epump_plant_eff = trans_eff * pump_eff

    for i in range (1,13):
        _pd_e = '{}{}'.format(pd_e, i)
        _pwd = '{}{}'.format(pwd, i)
        _sswd = '{}{}'.format(sswd, i)
        _ed_e = '{}{}'.format(ed_e, i)
        
        df[_pd_e]=(9.81*(df[_pwd]/1000)*df[tdh_gw])/Epump_plant_eff
        df[_ed_e]=(df[_sswd]*df[tdh_gw]*0.00272)/Epump_plant_eff
        
        if desalination:
            df[pd_e] += (df[_pwd]*df[des_int]*3600/1000)
            df[_ed_e] += (df['{}{}'.format(des_ener, i)]*1000000)
    return df