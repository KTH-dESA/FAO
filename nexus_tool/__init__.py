#Standard library imports
import pandas as pd
import numpy as np
from pandas import read_csv
from pandas import read_excel

#Local application/library specific imports
from nexus_tool.water_demand import (
    set_cropland_share,
    get_ky_list,
    get_kc_list,
    get_evap_i,
    get_eto,
    get_eff_rainfall_i,
    get_effective_rainfall,
    get_season_days,
    get_calendar_days,
    get_kc_values,
    get_water_demand,
)

from nexus_tool.energy_for_pumping import (
    get_gw_tdh,
    get_pumping_energy,
    get_annual_electricity,
)

from nexus_tool.lcoe import (
    get_wind_cf,
    get_pv_cf,
    get_installed_capacity,
    get_max_capacity,
    get_lcoe,
    get_least_cost,
)

class Model():
    # water properties:
    eto = 'ETo_'
    lat = 'lat'
    elevation = 'elevation'
    wind = 'wind'
    srad = 'srad'
    tmin = 'tmin'
    tmax = 'tmax'
    tavg = 'tavg'
    eff = 'eff_'
    prec = 'prec'
    kc = 'kc_'
    crop_share = 'crop_share'
    crop_area = 'crop_area'
    seasons = ['init', 'dev', 'mid', 'late']
    start = '_start'
    end = '_end'
    acwr = 'ACWR_' 
    pcwr = 'PCWR_'
    pwd = 'PWD_'
    sswd = 'SSWD_'
    crop_calendar = None
    crop_column = 'crop'
    ky_dict = {}
    kc_dict = {}
    pumping_hours_per_day = 10
    deff = 1
    aeff = 0.45
    ky_values = {}
    kc_values = {}
    # energy properties:
    gw_depth = 'gw_depth'
    tdh_gw = 'tdh_gw'
    des_int = 'Einten_KWh/m3'
    des_ener = 'Edesal_GWh_'
    pd_e = 'PD_E_'
    ed_e = 'ED_E_'
    trans_eff = 0
    pump_eff = 0
    technologies = {}
    discount_rate = 0
    start_year = 0
    end_year = 30
    def __init__(self, df, eto = eto, lat = lat, elevation = elevation,
                 wind = wind, srad = srad, tmin = tmin, tmax = tmax, 
                 tavg = tavg, crop_share = crop_share, crop_area = crop_area,
                 seasons = seasons, start = start, end = end, 
                 crop_calendar = crop_calendar, crop_column = crop_column,
                 pumping_hours_per_day = pumping_hours_per_day,
                 deff = deff, aeff = aeff, gw_depth = gw_depth, 
                 des_int = des_int, des_ener = des_ener, pd_e = pd_e,
                 ed_e = ed_e, trans_eff = trans_eff, pump_eff = trans_eff, ):
        self.df = df
        self.eto = eto
        self.lat = lat
        self.elevation = elevation
        self.wind = wind
        self.srad = srad
        self.tmin = tmin
        self.tmax = tmax
        self.tavg = tavg
        self.crop_share = crop_share
        self.crop_area = crop_area
        self.seasons = seasons
        self.start = start 
        self.end = end
        self.crop_calendar = crop_calendar
        self.crop_column = crop_column
        self.pumping_hours_per_day = pumping_hours_per_day
        self.deff = deff
        self.aeff = aeff
        self.gw_depth = gw_depth
        self.des_int = des_int
        self.des_ener = des_ener
        self.pd_e = pd_e
        self.ed_e = ed_e
        self.trans_eff = trans_eff
        self.pump_eff = trans_eff
            
    def print_properties(self):
        print('Properties names:')
        for val, name in zip([self.eto, self.lat, self.elevation, self.wind, 
                              self.srad, self.tmin, self.tmax, self.tavg,
                              self.crop_share, self.crop_area, self.seasons,
                              self.start, self.end, self.crop_column,
                              self.gw_depth, self.tdh_gw],
                             ['Reference evapotranspiration (.eto)', 
                              'Latitude (.lat)', 'Elevation (.elevation)', 
                              'Wind speed (.wind)', 'Solar radiation (.srad)', 
                              'Min temperature (.tmin)', 'Max temperature (.tmax)', 
                              'Avegarage temperature (.tavg)', 
                              'Cropland share column (.crop_share)', 
                              'Cropland area column (.crop_area)',
                              'Harvest seasons names (.seasons)', 
                              'Seasson start suffix (.start)',
                              'Seasson end suffix (.end)', 
                              'Cropland column (.crop_column)',
                              'Groundwater table depth (.gw_depth)',
                              'Total dynamic head (.tdh_gw)']):
            print('    - {}: {}'.format(name, val))
          
    ####### water related methods ###########
    def set_cropland_share(self, crop_var, geo_boundary = 'global', 
                           boundary_name = None, inplace = False):
        if inplace:
            set_cropland_share(self.df, crop_var, geo_boundary = geo_boundary, 
                       boundary_name = boundary_name, crop_share = self.crop_share)
        else:
            return set_cropland_share(self.df.copy(), crop_var, 
                                      geo_boundary = geo_boundary, 
                                      boundary_name = boundary_name, 
                                      crop_share = self.crop_share)
                                      
    def get_ky_list(self, inplace = False):
        if inplace:
            get_ky_list(self.df, crop_share = self.crop_share)
        else:
            return get_ky_list(self.df.copy(), crop_share = self.crop_share)
           
    def get_kc_list(self, inplace = False):
        if inplace:
            get_kc_list(self.df, crop_share = self.crop_share)
        else:
            return get_kc_list(self.df.copy(), crop_share = self.crop_share)
    
    def get_eto(self, inplace = False):
        if inplace:
            get_eto(self.df, eto = self.eto, lat = self.lat, 
                    elevation = self.elevation, wind = self.wind, 
                    srad = self.srad, tmin = self.tmin, 
                    tmax = self.tmax, tavg = self.tavg)
        else:
            return get_eto(self.df.copy(), eto = self.eto, lat = self.lat, 
                           elevation = self.elevation, wind = self.wind, 
                           srad = self.srad, tmin = self.tmin, 
                           tmax = self.tmax, tavg = self.tavg)
    
    def get_effective_rainfall(self, inplace = False):
        if inplace:
            get_effective_rainfall(self.df, eff = self.eff, prec = self.prec, 
                                   eto = self.eto)
        else:
            return get_effective_rainfall(self.df.copy(), eff = self.eff, 
                                          prec = self.prec, eto = self.eto)
                                          
    def get_calendar_days(self, inplace = False):
        if inplace:
            get_calendar_days(self.crop_calendar, seasons = self.seasons, 
                              start = self.start, end = self.end)
        else:
            return get_calendar_days(self.crop_calendar.copy(), seasons = self.seasons, 
                                     start = self.start, end = self.end)
                                     
    def get_kc_values(self, inplace = False):
        if inplace:
            get_kc_values(crop_calendar = self.crop_calendar, 
                          seasons = self.seasons, kc_dict = self.kc_dict,
                          crop_column = self.crop_column, start = self.start, 
                          end = self.end, kc = self.kc)
        else:
            return get_kc_values(crop_calendar = self.crop_calendar.copy(), 
                                 seasons = self.seasons, kc_dict = self.kc_dict,
                                 crop_column = self.crop_column,
                                 start = self.start, end = self.end, 
                                 kc = self.kc)
                                 
    def get_water_demand(self, inplace = False):
        if inplace:
            get_water_demand(self.df, self.crop_calendar, self.ky_dict, 
                             self.crop_column, self.aeff, self.deff, 
                             self.seasons[0], self.seasons[3], 
                             self.pumping_hours_per_day, 
                             crop_area = self.crop_area, _eto = self.eto, 
                             _kc = self.kc, _eff = self.eff, _acwr = self.acwr, 
                             _pcwr = self.pcwr, _pwd = self.pwd, _sswd = self.sswd, 
                             start = self.start, end = self.end, 
                             crop_share = self.crop_share)
        else:
            return get_water_demand(self.df.copy(), self.crop_calendar, 
                             self.ky_dict, self.crop_column, self.aeff, 
                             self.deff, self.seasons[0], self.seasons[3], 
                             self.pumping_hours_per_day, 
                             crop_area = self.crop_area, _eto = self.eto, 
                             _kc = self.kc, _eff = self.eff, _acwr = self.acwr, 
                             _pcwr = self.pcwr, _pwd = self.pwd, 
                             _sswd = self.sswd, start = self.start, 
                             end = self.end, crop_share = self.crop_share)
                             
    ####### energy related methods ###########
    def get_gw_tdh(self, inplace = False, wdd = 0, oap = 0, pld = 0):
        if inplace:
            get_gw_tdh(self.df, gw_depth = self.gw_depth, wdd = 0, oap = 0, pld = 0, 
                       interp_method = 'nearest', tdh_gw = self.tdh_gw)
        else:
            return get_gw_tdh(self.df.copy(), gw_depth = self.gw_depth, wdd = 0, oap = 0, 
                              pld = 0, interp_method = 'nearest', 
                              tdh_gw = self.tdh_gw)
                              
    def get_pumping_energy(self, inplace = False):
        if inplace:
            get_pumping_energy(self.df, self.trans_eff, self.pump_eff, 
                               pd_e = self.pd_e, pwd = self.pwd, 
                               sswd = self.sswd, ed_e = self.ed_e, 
                               tdh_gw = self.tdh_gw, desalination = False, 
                               des_int = self.des_int, des_ener = self.des_ener)
        else:
            return get_pumping_energy(self.df.copy(), self.trans_eff, self.pump_eff, 
                                      pd_e = self.pd_e, pwd = self.pwd, 
                                      sswd = self.sswd, ed_e = self.ed_e, 
                                      tdh_gw = self.tdh_gw, desalination = False, 
                                      des_int = self.des_int, 
                                      des_ener = self.des_ener)
                                      
    def get_annual_electricity(self, inplace = False):
        if inplace:
            get_annual_electricity(self.df, self.ed_e)
        else:
            return get_annual_electricity(self.df.copy(), self.ed_e)
                                      
    ####### technologies and LCOE related methods #########
    def create_wind_turbine(self, wind_turbine, life, om_cost, 
                            capital_cost, efficiency):
        self.technologies[wind_turbine] = self.WindTurbine(life, om_cost, 
                                                           capital_cost, 
                                                           efficiency)
                                                           
    def create_pv_system(self, pv_system, life, om_cost, 
                         capital_cost, efficiency):
        self.technologies[pv_system] = self.PVSystem(life, om_cost, 
                                                     capital_cost, 
                                                     efficiency, None, 
                                                     0, 1, 0, 0)
                                                     
    def create_standard_tech(self, tech_name, life, om_cost, capital_cost, 
                             efficiency, cf, fuel_cost, fuel_req, emission_factor, 
                             env_cost):
        self.technologies[tech_name] = self.Technology(life, om_cost, 
                                                capital_cost, efficiency, cf, 
                                                fuel_cost, fuel_req, 
                                                emission_factor, env_cost)
        
    def get_cf(self, technologies = 'all'):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            if type(self.technologies[technology]) == self.WindTurbine:
                self.get_wind_cf(technology)
            elif type(self.technologies[technology]) == self.PVSystem:
                self.get_pv_cf(technology)
    
    def get_wind_cf(self, wind_turbine):
        tech = self.technologies[wind_turbine]
        self.technologies[wind_turbine].cf = get_wind_cf(self.df, wind = self.wind, 
                    mu = tech.mu, t = tech.t, p_rated = tech.p_rated, 
                    z = tech.z, zr = tech.zr, es = tech.es, u_arr = tech.u_arr,
                    p_curve = tech.p_curve)
                    
    def get_pv_cf(self, pv_system):
        tech = self.technologies[pv_system]
        self.technologies[pv_system].cf = get_pv_cf(self.df, self.srad)
                    
    def get_installed_capacity(self, technologies = 'all'):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            self.technologies[technology].df = get_installed_capacity(self.df, 
                                                                      tech.cf, 
                                                                      self.pd_e)
                                                
    def get_max_capacity(self, technologies = 'all'):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            self.technologies[technology].df = tech.df.join(get_max_capacity(tech.df))
        
    def get_lcoe(self, technologies = 'all'):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            self.technologies[technology].df['lcoe'] = get_lcoe(
                                        max_capacity = tech.df['max_cap'],
                                        total_demand = self.df['annual_el_demand'],
                                        tech_life=tech.life, om_cost = tech.om_cost,
                                        capital_cost = tech.capital_cost,
                                        discount_rate = self.discount_rate,
                                        project_life = self.end_year - self.start_year,
                                        fuel_cost = tech.fuel_cost, 
                                        fuel_req = tech.fuel_req, 
                                        efficiency = tech.efficiency, 
                                        emission_factor = tech.emission_factor,
                                        env_cost = tech.env_cost)
                                        
    def get_least_cost(self,  technologies = 'all', geo_boundary = None):
        self.df['least_cost_tech'] = np.nan
        self.df['lcoe'] = np.nan
        if (geo_boundary != None) and (type(technologies) == dict):
            for key, value in technologies.items():
                _technologies = self.__check_tech_input(value)
                lcoe_df = pd.DataFrame()
                lcoe_df[geo_boundary] = self.df[geo_boundary]
                for _technology in _technologies:
                    lcoe_df[_technology] = self.technologies[_technology].df['lcoe']
                self.df.loc[self.df[geo_boundary]==key, 'least_cost_tech'], \
                self.df.loc[self.df[geo_boundary]==key, 'lcoe'] = \
                                    get_least_cost(lcoe_df, geo_boundary, key)
        else:
            _technologies = self.__check_tech_input(technologies)
            lcoe_df = pd.DataFrame()
            for _technology in _technologies:
                lcoe_df[_technology] = self.technologies[_technology].df['lcoe']
            self.df['least_cost_tech'], self.df['lcoe'] = \
                                    get_least_cost(lcoe_df)
                                      
    ####### additional methods #############
    def __check_tech_input(self, technologies):
        if type(technologies) == str:
            if technologies.lower() in ['all', 'a', 'everything']:
                technologies = self.technologies.keys()
            else:
                technologies = [technologies]
        return technologies
    
    def print_summary(self, geo_boundary = 'global'):
        if geo_boundary == 'global':
            temp_df = self.df.sum()
            summary = pd.DataFrame()
            
            summary['Irrigated area (ha)'] = [temp_df[self.crop_area]]
            summary['Water intensity (m3/ha)'] = [temp_df.filter(like=self.sswd).sum()/temp_df[self.crop_area]]
            summary['Water demand (Mm3)'] = [temp_df.filter(like=self.sswd).sum()/1000000]
            
            summary.index = ['Global']
        else:
            temp_df = self.df.groupby(geo_boundary).sum()
            summary = pd.DataFrame()
            
            summary['Irrigated area (ha)'] = temp_df[self.crop_area]
            summary['Water intensity (m3/ha)'] = temp_df.filter(like=self.sswd).sum(axis=1)/temp_df[self.crop_area]
            summary['Water demand (Mm3)'] = temp_df.filter(like=self.sswd).sum(axis=1)/1000000
            summary['Energy demand (GWh)'] = temp_df.filter(like=self.ed_e).sum(axis=1)/1000000
        
        summary.round(decimals=3)
        return summary
            
    class Technology():
        df = pd.DataFrame()
        fuel_cost = 0
        fuel_req = 0
        efficiency = 1
        emission_factor = 0
        env_cost = 0
        def __init__(self, life, om_cost, capital_cost, efficiency, cf,
                     fuel_cost, fuel_req, emission_factor, env_cost):
            self.life = life
            self.om_cost = om_cost
            self.capital_cost = capital_cost
            self.efficiency = efficiency
            self.cf = cf
            self.fuel_cost = fuel_cost
            self.fuel_req = fuel_req
            self.emission_factor = emission_factor
            self.env_cost = env_cost

    
    class WindTurbine(Technology):
        # properties:
        mu = 0.97  # availability factor
        t = 24*30
        p_rated = 600
        z = 55  # hub height
        zr = 80  # velocity measurement height
        es = 0.85  # losses in wind electricity
        u_arr = range(1, 26)
        p_curve = [0, 0, 0, 0, 30, 77, 135, 208, 287, 371, 450, 514, 558,
                   582, 594, 598, 600, 600, 600, 600, 600, 600, 600, 600, 600]
        def __init__(self, life, om_cost, capital_cost, efficiency, mu = mu, 
                     t = t, p_rated = p_rated, z = z, zr = zr, es = es, 
                     u_arr = u_arr, p_curve = p_curve):
            super().__init__(life, om_cost, capital_cost, efficiency, None, 
                             0, 1, 0, 0)
            self.mu = mu
            self.t = t
            self.p_rated = p_rated
            self.z = z
            self.zr = zr
            self.es = es
            self.u_arr = u_arr
            self.p_curve = p_curve
           
    class PVSystem(Technology):
        pass
            
            
            
            
            
            
            
            