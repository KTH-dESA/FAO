#Standard library imports
import pandas as pd
import numpy as np
from math import pi

def wind_cf(df, wind, mu, t, p_rated, z, zr, es, u_arr, p_curve):
    u_zr = df[wind]
    
    # Adjust for the correct hub height
    alpha = (0.37 - 0.088 * np.log(u_zr)) / (1 - 0.088 * np.log(zr / 10))
    u_z = u_zr * (z / zr) ** alpha

    # Rayleigh distribution and sum of series
    rayleigh = [(pi / 2) * (u / u_z ** 2) * np.exp((-pi / 4) * (u / u_z) ** 2) for u in u_arr]
    energy_produced = np.array(sum([mu * es * t * p * r for p, r in zip(p_curve, rayleigh)]))
    
    return energy_produced/(p_rated * t)
    
def get_wind_cf(df, wind, mu, t, p_rated, z, zr, es, u_arr, p_curve):
    cf_df = pd.DataFrame()
    for i in range (1,13):
        wind='wind{}'.format(i)
        cf_df['cf_{}'.format(i)] = wind_cf(df, wind, mu, t, p_rated, z, 
                                             zr, es, u_arr, p_curve)
    return cf_df
    
def get_installed_capacity(df, cf, pd_e):
    ic_df = pd.DataFrame()
    for i in range(1,13):
        _cf = cf if type(cf) == float else cf[f'cf_{i}']
        ic_df[f'ic_{i}'] = df[f'{pd_e}{i}'] / _cf
    return ic_df
    
def get_max_capacity(df):
    return pd.DataFrame({'max_cap': df.filter(like='ic_').max(axis=1)})
