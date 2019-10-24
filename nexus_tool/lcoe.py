#Standard library imports
import pandas as pd
import numpy as np

#### default values:
mu = 0.97  # availability factor
t = 24*30
p_rated = 600
z = 55  # hub height
zr = 80  # velocity measurement height
es = 0.85  # losses in wind electricity
u_arr = range(1, 26)
p_curve = [0, 0, 0, 0, 30, 77, 135, 208, 287, 371, 450, 514, 558,
           582, 594, 598, 600, 600, 600, 600, 600, 600, 600, 600, 600]

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
    for i in range (1,13):
    wind='wind{}'.format(i)
    df['wind_cf_{}'.format(i)] = get_wind_cf(data, wind, mu, t, p_rated, z, 
                                             zr, es, u_arr, p_curve)