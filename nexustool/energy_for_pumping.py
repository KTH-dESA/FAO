#Standard library imports
import numpy as np
from math import pi


#D in m and A in m2
def get_A(D):
    A = (pi*D**2)/4
    return A

#Q in m3/sec , A in m2 and V in m/sec
def get_V(df, avg_Q, A, mV, pump_hours, axis=0):
    if axis:
        for i in range (1,13):
            _avg_Q = '{}{}'.format(avg_Q, i)
            _mV = '{}{}'.format(mV, i)
            
            df[_mV] = df[_avg_Q]/(30*pump_hours*60*60)/A 
            
            # m3/month / (30d/m * 10h/day * 60min/h * 60s/h)
    else:
        df[mV] = df[avg_Q]/(30*pump_hours*60*60)/A  # convert m3/month to m3/s (30day*pump_hours_day*60min/h*60s/min)
    
    return df

#Re=Reynold number (unitless), Ken_visc=Kinematic viscosity (m2 /s) 
def get_Re(df, Re, mV, D, Ken_visc, axis=0):
    if axis:
        for i in range (1,13):
            _mV = '{}{}'.format(mV, i)
            _Re = '{}{}'.format(Re, i) #already defined in the class properties Re='Re_'
            
            df[_Re] = (df[_mV]*D)/Ken_visc
    else:
        df[Re] = df[mV]*D/Ken_visc
        
    return df

#f=friction coefficient (unitless), k =  Roughness factor (m)
def get_f(df, f, k, D, Re, axis=0):
    if axis:
        for i in range (1,13):
            _Re = '{}{}'.format(Re, i)
            _f = '{}{}'.format(f, i) #already defined in the class properties f='f_'
        
            df[_f] =0.25/(np.log((k/(3.7*D))+(5.74/(df[_Re]**0.9)))**2)
        
    else:
        df[f] = 0.25/(np.log((k/(3.7*D))+(5.74/(df[Re]**0.9)))**2)

    return df

#tds in (m)
def get_tdh(df, tdh, elevation, f, L, sswd, D, g, pump_hours, friction=True, axis=0):
    if axis:
        for i in range (1,13):
            _avg_Q = '{}{}'.format(sswd, i)
            _f = '{}{}'.format(f, i)
            _tdh_sw = '{}{}'.format(tdh, i)

            if friction:
                f_losses = ((df[_f]*L*16*((df[_avg_Q]/(30*pump_hours*60*60))**2))/((D**5)*2*g*(pi**2)))
            else:
                f_losses = 0
            df[_tdh_sw] = df[elevation] + f_losses
            df[_tdh_sw].replace(0, np.nan, inplace=True)
    else:
        if friction:
            f_losses = ((df[f] * L * 16 * ((df[sswd] / (30 * pump_hours * 60 * 60)) ** 2)) / ((D ** 5) * 2 * g * (pi ** 2)))
        else:
            f_losses = 0
        df[tdh] = df[elevation] + f_losses

    return df


#P=  Power in (W), dens=Density (Kg/m3), g=gravitational acceleration in (m/sec2)
def get_pumping_energy(df, tdh, pump_eff, pp_e, pa_e, g, pwd,
                       sswd, dens, axis=0):
    if axis:
        pump_eff = pump_eff
        
        for i in range (1,13):
            _swpp_e = '{}{}'.format(pp_e, i) #surface water pumping peak electric demand
            _swpa_e = '{}{}'.format(pa_e, i) #surface water pumping average electric demand
            _peak_Q = '{}{}'.format(pwd, i) #peak water flow in the pipeline. To be updated WEAP output.
            _avg_Q = '{}{}'.format(sswd, i) #average water flow in the pipeline. To be updated with WEAP output
            _tdh_sw = '{}{}'.format(tdh, i)
            
            df[_swpp_e]=((df[_peak_Q]*df[_tdh_sw]*g*dens) / (pump_eff * 1000)) #to convert E from W to KW
            df[_swpa_e]=((df[_avg_Q]*df[_tdh_sw]*g*dens) / (pump_eff * 1000 * 3600)) #to convert E from J to KWh
            
    else:
        df[pp_e] = ((df[pwd] * df[tdh] * g * dens) / (pump_eff * 1000))
        df[pa_e] = ((df[sswd] * df[tdh] * g * dens) / (pump_eff * 1000 * 3600))
        
    return df
