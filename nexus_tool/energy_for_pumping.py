#Standard library imports
import pandas as pd
import numpy as np

def get_gw_tdh(df, gw_depth, wdd, oap, pld, tdh_gw, interp_method = 'nearest'):
    df[tdh_gw] = df[gw_depth] + wdd + oap + pld
    df[tdh_gw].replace(0, np.nan, inplace=True)
    # df[tdh_gw].interpolate(method = interp_method, axis=0, inplace=True)
    return df

def get_A(pi,D):
    A = (pi*D**2)/4
    return A

def get_V(Q,A):
    V = Q/A
    return V

def get_Re(V,D,Ken_visc):
    Re = (V*D)/Ken_visc
    return Re

def get_f(k,D,Re):
    f =0.25/(np.log((k/(3.7*D))+(5.74/(Re*0.9)))**2)
    return f

def get_sw_tdh(df, tdh_sw, elevation, f, L, Q, D, g, pi, interp_method = 'nearest'):
    df[tdh_sw] = df[elevation] + ((f*L*16*(Q**2))/((D**5)*2*g*(pi**2)))  #the elevation should be delta elevation between the two points
    df[tdh_sw].replace(0, np.nan, inplace=True)
    # df[tdh_sw].interpolate(method = interp_method, axis=0, inplace=True)
    return df

#def get_pumping_energy(df, trans_eff, pump_eff, pd_e, pwd, sswd, ed_e, tdh_gw, tdh_sw,
#                       des_int, des_ener, desalination = False):
#    Epump_plant_eff = trans_eff * pump_eff
#
#    for i in range (1,13):
#        _pd_e = '{}{}'.format(pd_e, i)
#        _pwd = '{}{}'.format(pwd, i)
#        _sswd = '{}{}'.format(sswd, i)
#        _ed_e = '{}{}'.format(ed_e, i)
#        
#        df[_pd_e]=(9.81*(df[_pwd]/1000)*df[tdh_gw])/Epump_plant_eff
#        df[_ed_e]=(df[_sswd]*df[tdh_gw]*0.00272)/Epump_plant_eff
#        
#        
#        if desalination:
#            df[_pd_e] += (df[_pwd]*df[des_int]*3600/1000)
#            df[_ed_e] += (df['{}{}'.format(des_ener, i)]*1000000)
#    return df

#def get_annual_electricity(df, ed_e):
#    df['annual_el_demand'] = df.filter(like=ed_e).sum(axis=1)
#    return df

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


def get_SWpumping_energy(df, tdh_sw, SWpump_eff, swpp_e, swpa_e, g, pwd, sswd, dens):
    SWpump_eff = SWpump_eff
    
    for i in range (1,13):
        _swpp_e = '{}{}'.format(swpp_e, i) #surface water pumping peak electric demand 
        _swpa_e = '{}{}'.format(swpa_e, i) #surface water pumping average electric demand
        _pQ = '{}{}'.format(pwd, i) #peak water flow in the pipeline. To be updated WEAP output. 
        _aQ = '{}{}'.format(sswd, i) #average water flow in the pipeline. To be updated with WEAP output 
        
        
        df[_swpp_e]=(df[_pQ]*df[tdh_sw]*g*dens)/SWpump_eff
        df[_swpa_e]=(df[_aQ]*df[tdh_sw]*g*dens)/SWpump_eff
        
    return df

def get_totalpumping_energy(df, swpa_e, ed_e):
    for i in range (1,13):
        df[total_pumping_energy]=(df[_swpa_e]+df[_ed_e])     
    
    return df

def get_annual_electricity(df, ed_e):
    df['annual_el_demand'] = df.filter(like=total_pumping_energy).sum(axis=1)
    return df