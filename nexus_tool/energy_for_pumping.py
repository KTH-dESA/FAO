#Standard library imports
import pandas as pd
import numpy as np
from math import pi

def get_gw_tdh(df, gw_depth, wdd, oap, pld, tdh_gw, interp_method = 'nearest'):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw].replace(0, np.nan, inplace=True)
    # df[tdh_gw].interpolate(method = interp_method, axis=0, inplace=True)
    return df
#D in m and A in m2
def get_A(D):
    A = (pi*D**2)/4
    return A

#Q in m3/sec , A in m2 and V in m/sec
def get_V(df, avg_Q, A, mV, pump_hours, axis=1):
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
def get_Re(df,Re,mV,D, Ken_visc, axis=1):
    if axis:
        for i in range (1,13):
            _mV = '{}{}'.format(mV, i)
            _Re = '{}{}'.format(Re, i) #already defined in the class properties Re='Re_'
            
            df[_Re] = (df[_mV]*D)/Ken_visc
    else:
        df[Re] = df[mV]*D/Ken_visc
        
    return df

#f=friction coefficient (unitless), k =  Roughness factor (m)
def get_f(df,f, k,D,Re, axis=1):
    if axis:
        for i in range (1,13):
            _Re = '{}{}'.format(Re, i)
            _f = '{}{}'.format(f, i) #already defined in the class properties f='f_'
        
            df[_f] =0.25/(np.log((k/(3.7*D))+(5.74/(df[_Re]**0.9)))**2)
        
    else:
        df[f] = 0.25/(np.log((k/(3.7*D))+(5.74/(df[Re]**0.9)))**2)

    return df

#tds in (m)
def get_sw_tdh(df, tdh_sw, elevation, f, L, avg_Q, D, g, pump_hours, axis=1):
    if axis:
        for i in range (1,13):
            _avg_Q = '{}{}'.format(avg_Q, i)
            _f = '{}{}'.format(f, i)
            _tdh_sw = '{}{}'.format(_tdh_sw, i)
        
            df[_tdh_sw] = (df[elevation] + ((df[_f]*L*16*((df[_avg_Q]/(30*pump_hours*60*60))**2))/((D**5)*2*g*(pi**2))))
            df[_tdh_sw].replace(0, np.nan, inplace=True) #this line is unecesary and check the above one as it is being overwriten every time
    else:
        df[tdh_sw] = (df[elevation] + ((df[f]*L*16*((df[avg_Q]/(30*pump_hours*60*60))**2))/((D**5)*2*g*(pi**2))))

    return df

def get_GWpumping_energy(df, trans_eff, pump_eff, pd_e, pwd, sswd, ed_e, tdh_gw, 
                       des_int, des_ener, desalination = False):
    GWpump_plant_eff = trans_eff * pump_eff
    
    for i in range (1,13):
        _pd_e = '{}{}'.format(pd_e, i) #in kW
        _pwd = '{}{}'.format(pwd, i) #in l/s
        _sswd = '{}{}'.format(sswd, i) #in m3
        _ed_e = '{}{}'.format(ed_e, i) #in kWh
        
        df[_pd_e]=(9.81*(df[_pwd]/1000)*df[tdh_gw])/GWpump_plant_eff # this should be changed to use m3/s instead of l/s
        df[_ed_e]=(df[_sswd]*df[tdh_gw]*0.00272)/GWpump_plant_eff 
        
        if desalination:
            df[_pd_e] += (df[_pwd]*df[des_int]*3600/1000)
            df[_ed_e] += (df['{}{}'.format(des_ener, i)]*1000000)
    return df


#P=  Power in (W), dens=Density (Kg/m3), g=gravitational acceleration in (m/sec2)
def get_SWpumping_energy(df, tdh_sw, SWpump_eff, swpp_e, swpa_e, g, peak_Q, 
                         avg_Q, dens, axis=1):
    if axis:   
        SWpump_eff = SWpump_eff
        
        for i in range (1,13):
            _swpp_e = '{}{}'.format(swpp_e, i) #surface water pumping peak electric demand 
            _swpa_e = '{}{}'.format(swpa_e, i) #surface water pumping average electric demand
            _peak_Q = '{}{}'.format(peak_Q, i) #peak water flow in the pipeline. To be updated WEAP output. 
            _avg_Q = '{}{}'.format(avg_Q, i) #average water flow in the pipeline. To be updated with WEAP output 
            
            
            df[_swpp_e]=((df[_peak_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000)) #to convert E from W to KW
            df[_swpa_e]=((df[_avg_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000*3600)) #to convert E from J to KWh
            
    else:
        df[swpp_e] = ((df[peak_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000))
        df[swpa_e] = ((df[avg_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000*3600))
        
        #m3/month * m * m/s2 * kg/m3 / (1000*3600s/h)
        
    return df
        
def get_total_pumping_energy(df, swpa_e, ed_e):
    for i in range (1,13):
        _ed_e = '{}{}'.format(ed_e, i)
        _swpa_e = '{}{}'.format(swpa_e, i)
        total_pumping_energy ='{}{}'.format('total_pumping_energy', i)
        
        df[total_pumping_energy]=(df[_swpa_e]+df[_ed_e])     
    
    return df

def get_annual_electricity(df, ed_e):
    
    df['annual_el_demand'] = df.filter(like='total_pumping_energy').sum(axis=1)
    return df