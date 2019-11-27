#Standard library imports
import pandas as pd
import numpy as np

def get_gw_tdh(df, gw_depth, wdd, oap, pld, tdh_gw, interp_method = 'nearest'):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw].replace(0, np.nan, inplace=True)
    # df[tdh_gw].interpolate(method = interp_method, axis=0, inplace=True)
    return df
#D in m and A in m2
def get_A(pi,D):
    A = (pi*D**2)/4
    return A

#Q in m3/sec , A in m2 and V in m/sec
def get_V(df,avg_Q,A,mV):
    for i in range (1,13):
        _avg_Q = '{}{}'.format(avg_Q, i)
        _mV = '{}{}'.format(mV, i)
        
        df[_mV] = df[_avg_Q]/A
    
    return df

#Re=Reynold number (unitless), Ken_visc=Kinematic viscosity (m2 /s) 
def get_Re(df,Re,mV,D, Ken_visc):
    for i in range (1,13):
        _mV = '{}{}'.format(mV, i)
        _Re = '{}{}'.format(Re, i) #already defined in the class properties Re='Re_'
        
        df[_Re] = (df[_mV]*D)/Ken_visc
    return df

#f=friction coefficient (unitless), k =  Roughness factor (m)
def get_f(df,f, k,D,Re):
    for i in range (1,13):
        _Re = '{}{}'.format(Re, i)
        _f = '{}{}'.format(f, i) #already defined in the class properties f='f_'
    
        df[_f] =0.25/(np.log((k/(3.7*D))+(5.74/(df[_Re]*0.9)))**2)
    
    return df

#tds in (m)
def get_sw_tdh(df, tdh_sw, elevation, f, L, avg_Q, D, g, pi, interp_method = 'nearest'):
    for i in range (1,13):
        _avg_Q = '{}{}'.format(avg_Q, i)
        _f = '{}{}'.format(f, i)
    
        df[tdh_sw] = (df[elevation] + ((df[_f]*L*16*((df[_avg_Q])**2))/((D**5)*2*g*(pi**2))))  #the elevation should be delta elevation between the two points
        df[tdh_sw].replace(0, np.nan, inplace=True)
    # df[tdh_sw].interpolate(method = interp_method, axis=0, inplace=True)
    return df


def get_GWpumping_energy(df, trans_eff, pump_eff, pd_e, pwd, sswd, ed_e, tdh_gw, 
                       des_int, des_ener, desalination = False):
    GWpump_plant_eff = trans_eff * pump_eff
    
    for i in range (1,13):
        _pd_e = '{}{}'.format(pd_e, i)
        _pwd = '{}{}'.format(pwd, i)
        _sswd = '{}{}'.format(sswd, i)
        _ed_e = '{}{}'.format(ed_e, i)
        
        df[_pd_e]=(9.81*(df[_pwd]/1000)*df[tdh_gw])/GWpump_plant_eff
        df[_ed_e]=(df[_sswd]*df[tdh_gw]*0.00272)/GWpump_plant_eff
        
        if desalination:
            df[_pd_e] += (df[_pwd]*df[des_int]*3600/1000)
            df[_ed_e] += (df['{}{}'.format(des_ener, i)]*1000000)
    return df


#P=  Power in (W), dens=Density (Kg/m3), g=gravitational acceleration in (m/sec2)
def get_SWpumping_energy(df, tdh_sw, SWpump_eff, swpp_e, swpa_e, g, pwd, avg_Q, dens):
    SWpump_eff = SWpump_eff
    
    for i in range (1,13):
        _swpp_e = '{}{}'.format(swpp_e, i) #surface water pumping peak electric demand 
        _swpa_e = '{}{}'.format(swpa_e, i) #surface water pumping average electric demand
        _peak_Q = '{}{}'.format(pwd, i) #peak water flow in the pipeline. To be updated WEAP output. 
        _avg_Q = '{}{}'.format(avg_Q, i) #average water flow in the pipeline. To be updated with WEAP output 
        
        
        df[_swpp_e]=((df[_peak_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000)) #to convert P from W to KW
        df[_swpa_e]=((df[_avg_Q]*df[tdh_sw]*g*dens)/(SWpump_eff*1000)) #to convert P from W to KW
        
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