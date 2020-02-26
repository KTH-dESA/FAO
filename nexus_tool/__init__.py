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
    get_A,
    get_V,
    get_Re,
    get_f,
    get_sw_tdh,
    get_GWpumping_energy,
    get_SWpumping_energy,
    get_total_pumping_energy,
    get_annual_electricity,
)

from nexus_tool.least_cost import (
    get_wind_cf,
    get_pv_cf,
    get_installed_capacity,
    get_max_capacity,
    get_lcoe,
    get_least_cost,
    get_tech_generation,
    get_pumping_cost,
    get_unit_pumping_cost,
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
    D='Pipe_diameter' #new
    L='Pipeline_length' #new
    A='Pipe_area' #new
    mV='Flow_velocity_' #new, monthly Velocity mV in (m/sec)
    tdh_sw = 'tdh_sw' #new
    des_int = 'Einten_KWh/m3'
    des_ener = 'Edesal_GWh_'
    pd_e = 'PD_E_'
    ed_e = 'ED_E_'
    peak_Q = 'PWD_' #new, Peak flow in the pipeline, WEAP
    avg_Q = 'SSWD_'   #new, Average flow in the pipeline, WEAP
    swpp_e= 'SWPP_E_' #new, Surface Water Pumping PEAK Electric demand
    swpa_e= 'SWPA_E_' #new, Surface Water Pumping AVERAGE Electric demand
    Re='Re_' #new
    f='f_' #new
    trans_eff = 0
    pump_eff = 0
    SWpump_eff =0.6 #new
    technologies = {}
    discount_rate = 0
    g= 9.81 #new gravitational acceleration in(m/sec2)
    Ken_visc=1.004e-06 #new
    dens = 1000 #new
    k=0.26  #New - roughness for cast iron
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
                 ed_e = ed_e, swpp_e = swpp_e, swpa_e = swpa_e,  trans_eff = trans_eff, SWpump_eff=SWpump_eff, pump_eff = pump_eff, ):
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
        self.swpp_e = swpp_e #new
        self.swpa_e = swpa_e #new
        self.trans_eff = trans_eff
        self.SWpump_eff = SWpump_eff
        self.pump_eff = pump_eff #changed from self.pump_eff = trans_eff
            
    def print_properties(self):
        print('Properties names:')
        for val, name in zip([self.eto, self.lat, self.elevation, self.wind, 
                              self.srad, self.tmin, self.tmax, self.tavg,
                              self.crop_share, self.crop_area, self.seasons,
                              self.start, self.end, self.crop_column,
                              self.gw_depth, self.tdh_gw,
                              self.tdh_sw,],
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
                              'Total dynamic head ground water (.tdh_gw)',
                              'Total dynamic head surface water (.tdh_sw)',
                              'Pipe area (.A)',
                              'Flow velocity (.V)',
                              'Friction losses (.f)',
                              'Reynolds Number (.Re)',]):
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
                              
    def get_GWpumping_energy(self, inplace = False):
        if inplace:
            self.GWpumping_energy=get_GWpumping_energy(self.df, self.trans_eff, self.pump_eff, 
                               pd_e = self.pd_e, pwd = self.pwd, 
                               sswd = self.sswd, ed_e = self.ed_e, 
                               tdh_gw = self.tdh_gw, desalination = False, 
                               des_int = self.des_int, des_ener = self.des_ener)
        else:
            return get_GWpumping_energy(self.df.copy(), self.trans_eff, self.pump_eff, 
                                      pd_e = self.pd_e, pwd = self.pwd, 
                                      sswd = self.sswd, ed_e = self.ed_e, 
                                      tdh_gw = self.tdh_gw, desalination = False, 
                                      des_int = self.des_int, 
                                      des_ener = self.des_ener)
    
    def get_A(self, inplace=False):
        if inplace:
            self.df[self.A]= get_A(D=self.df[self.D])
            
        else:
            return get_A(D=self.df[self.D])
    
    def get_V(self, inplace=False, axis=1):
        if inplace:
            self.df=get_V(self.df, avg_Q=self.avg_Q, A=self.df[self.A], 
                          mV=self.mV, pump_hours = self.pumping_hours_per_day, 
                          axis=axis)
            
        else:
            return get_V(self.df.copy(), avg_Q=self.avg_Q, A=self.df[self.A], 
                         mV=self.mV, axis=axis)
    
    
    def get_Re(self, inplace=False, axis=1):
        if inplace: 
            self.df=get_Re(self.df, Re=self.Re, mV=self.mV, D=self.df[self.D], 
                           Ken_visc=1000, axis=axis)
            
        else:
            return get_Re(self.df.copy(), Re=self.Re, mV=self.mV, 
                          D=self.df[self.D], Ken_visc=1000, axis=axis)
    
    def get_f(self, inplace=False, axis=1):
        if inplace:
            self.df=get_f(self.df, f=self.f, k=0.26, D=self.df[self.D], 
                          Re=self.Re, axis=axis)
        else:
            return get_f(self.df.copy(), f=self.f, k=0.26, D=self.df[self.D], 
                         Re=self.Re, axis=axis)
    
                                   
    def get_sw_tdh(self, inplace = False, axis=1):
        if inplace:
            self.df=get_sw_tdh(self.df, tdh_sw=self.tdh_sw, 
                               elevation=self.elevation, f =self.f, 
                               L=self.df[self.L], avg_Q=self.avg_Q, 
                               D=self.df[self.D], g= 9.81,
                               pump_hours = self.pumping_hours_per_day, 
                               axis=axis)
        else:
            return get_sw_tdh(self.df.copy(), tdh_sw=self.tdh_sw, 
                              elevation=self.elevation, f =self.f, 
                              L=self.df[self.L], Q=self.avg_Q, 
                              D=self.df[self.D], g= 9.81, 
                              pump_hours = self.pumping_hours_per_day, 
                              axis=axis)   
    
    
    def get_SWpumping_energy(self, inplace = False, axis=1):
        if inplace:
            self.SWpumping_energy=get_SWpumping_energy(self.df, 
                    SWpump_eff = self.SWpump_eff, tdh_sw = self.tdh_sw, 
                    swpp_e = self.swpp_e, peak_Q = self.peak_Q, swpa_e = self.swpa_e,
                    avg_Q = self.avg_Q, g=self.g, dens=self.dens, axis=axis)
        else:
            return get_SWpumping_energy(self.df.copy(), 
                           SWpump_eff=self.SWpump_eff, tdh_sw = self.tdh_sw, 
                           swpp_e = self.swpp_e, peak_Q = self.peak_Q, 
                           swpa_e = self.swpa_e, sswd = self.sswd, g=self.g, 
                           dens=self.dens, axis=axis)
    
    def get_total_pumping_energy(self, inplace =False):
        if inplace:
            get_total_pumping_energy(self.df, swpa_e = self.swpa_e, ed_e = self.ed_e)
        else:
            return get_total_pumping_energy(self.df.copy(), swpa_e=self.swpa_e, ed_e = self.ed_e)
    
    
    
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
        
    def get_cf(self, technologies = 'all', axis=1):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            if type(self.technologies[technology]) == self.WindTurbine:
                self.get_wind_cf(technology, axis)
            elif type(self.technologies[technology]) == self.PVSystem:
                self.get_pv_cf(technology, axis)
    
    def get_wind_cf(self, wind_turbine, axis=1):
        tech = self.technologies[wind_turbine]
        self.technologies[wind_turbine].cf = get_wind_cf(self.df, wind = self.wind, 
                    mu = tech.mu, t = tech.t, p_rated = tech.p_rated, 
                    z = tech.z, zr = tech.zr, es = tech.es, u_arr = tech.u_arr,
                    p_curve = tech.p_curve, axis = axis)
                    
    def get_pv_cf(self, pv_system, axis=1):
        tech = self.technologies[pv_system]
        self.technologies[pv_system].cf = get_pv_cf(self.df, self.srad, axis)
                    
    def get_installed_capacity(self, technologies = 'all', axis=1):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            self.technologies[technology].df = get_installed_capacity(self.df, 
                                                                      tech.cf, 
                                                                      self.pd_e,
                                                                      axis)
                                                
    def get_max_capacity(self, technologies = 'all', axis=1):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            if axis:
                self.technologies[technology].df = tech.df.join(get_max_capacity(tech.df, axis))
            else:
                self.technologies[technology].max_cap = get_max_capacity(tech.df, axis)
        
    def get_lcoe(self, technologies = 'all', years = 'all', axis=1):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            if axis:
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
                                            env_cost = tech.env_cost,
                                            start_year = self.start_year,
                                            end_year = self.end_year,
                                            axis = axis)
            else:
                years = self.__get_years(years)
                self.technologies[technology].lcoe = pd.DataFrame()
                for year in years:
                    self.technologies[technology].lcoe = \
                            self.technologies[technology].lcoe.append(get_lcoe(
                                max_capacity = tech.max_cap.reset_index(),
                                total_demand = self.df,
                                tech_life=tech.life, om_cost = tech.om_cost,
                                capital_cost = tech.capital_cost,
                                discount_rate = self.discount_rate,
                                project_life = self.end_year - self.start_year,
                                fuel_cost = tech.fuel_cost, 
                                fuel_req = tech.fuel_req, 
                                efficiency = tech.efficiency, 
                                emission_factor = tech.emission_factor,
                                env_cost = tech.env_cost,
                                start_year = year,
                                end_year = self.end_year,
                                axis = axis), ignore_index=True)
                                
                # self.technologies[technology].lcoe.reset_index(inplace=True)
                                                            
    def get_least_cost(self,  technologies = 'all', years = 'all',
                       geo_boundary = None, axis=1):
        if axis:
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
                self.df['least_cost_tech'], self.df['lcoe'] = get_least_cost(lcoe_df)
        else:
            _technologies = self.__check_tech_input(technologies)
            lcoe_df = pd.DataFrame()
            for _technology in _technologies:
                dff = self.technologies[_technology].lcoe.set_index(['Demand point', 'year'])
                lcoe_df[_technology] = dff['lcoe']
            
            years = self.__get_years(years)
            self.lcoe = self.df.loc[self.df.Year.isin(years)]
            self.lcoe = self.lcoe.groupby(['Demand point', 'Year']).agg(
                                                      {
                                                       # 'Supply point': 'first',
                                                       # 'links': 'first',
                                                       'province': 'first',
                                                       'sswd': 'sum',
                                                       # 'type': 'first',
                                                       'swpp_e': 'max',
                                                       'swpa_e': 'sum'})
            self.lcoe.rename(columns={'links': 'link', 'sswd': 'water demand',
                                      'swpp_e': 'required capacity',
                                      'swpa_e': 'energy demand'},
                            inplace=True)
            lcoe = get_least_cost(lcoe_df)
            
            self.lcoe['least_cost_technology'] = lcoe['least_cost_technology']
            self.lcoe['lcoe'] = lcoe['lcoe']
    
    def get_tech_generation(self):
        get_tech_generation(self.df, self.technologies.keys())
        
    def get_pumping_cost(self, inplace = False):
        if inplace:
            get_pumping_cost(self.df, 'annual_el_demand', 'lcoe')
        else:
            return get_pumping_cost(self.df.copy(), 'annual_el_demand', 'lcoe')
            
    def get_unit_pumping_cost(self, inplace = False):
        if inplace:
            get_unit_pumping_cost(self.df, 'pumping_cost',
                                  self.df.filter(like=self.sswd).sum(axis=1))
        else:
            return get_unit_pumping_cost(self.df.copy(), 'pumping_cost',
                                self.df.filter(like=self.sswd).sum(axis=1))
                                      
    ####### additional methods #############
    def __check_tech_input(self, technologies):
        if type(technologies) == str:
            if technologies.lower() in ['all', 'a', 'everything']:
                technologies = self.technologies.keys()
            else:
                technologies = [technologies]
        return technologies
    
    def __get_years(self, years):
        if type(years) == str:
            if years.lower() in ['all', 'a', 'everything']:
                years = range(self.start_year, self.end_year + 1)
        elif type(years) == int:
            years = [years]
        return years
    
    def print_summary(self, geo_boundary = 'global'):
        if 'month' in geo_boundary:
            _id_vars = [geo_boundary] if type(geo_boundary) == str else geo_boundary.copy()
            _id_vars.remove('month')
            temp_df = self.df.melt(id_vars=['crop_area']+_id_vars, 
                        value_vars=self.df.columns[self.df.columns.str.contains(self.sswd)])
            temp_df.rename(columns={'crop_area': 'Irrigated area (ha)', 
                                    'variable': 'month', 
                                    'value': 'Water demand (Mm3)'}, 
                                    inplace=True)
            for i in range(1,13):
                temp_df.loc[temp_df['month']==f'{self.sswd}{i}','month'] = i
            
            temp_df['energy demand'] = self.df.melt(value_vars=self.df.columns[self.df.columns.str.contains(self.ed_e)])['value']
            summary = temp_df.groupby(geo_boundary).agg({'Irrigated area (ha)': 'sum',
                                                         'Water demand (Mm3)': lambda row: sum(row)/1000000})
            summary.insert(2, 'Water intensity (m3/ha)', 
                           summary['Water demand (Mm3)'] * 1000000 / \
                           summary['Irrigated area (ha)'])
            try:
                temp_df['Energy demand (GWh)'] = self.df.melt(value_vars=self.df.columns[self.df.columns.str.contains(self.total_pumping_energy)])['value']/1000000
                summary = summary.join(temp_df.groupby(geo_boundary)['Energy demand (GWh)'].sum())
            except:
                pass           
                
        else:
        
            temp_df = pd.DataFrame()
            temp_df['Irrigated area (ha)'] = self.df[self.crop_area]
            
            if geo_boundary == 'global':
                temp_df[geo_boundary] = geo_boundary
            else:
                temp_df[geo_boundary] = self.df[geo_boundary]
            
            summary = temp_df.groupby(geo_boundary).agg({'Irrigated area (ha)': 'sum'})
                
            ### water related data
            try:
                temp_df['Water demand (Mm3)'] = self.df.filter(like=self.sswd).sum(axis=1) / 1000000
                summary = summary.join(temp_df.groupby(geo_boundary)['Water demand (Mm3)'].sum())
                summary.insert(2, 'Water intensity (m3/ha)', 
                               summary['Water demand (Mm3)'] * 1000000 / \
                               summary['Irrigated area (ha)'])
            except:
                pass
                
            ### energy related data
            try:
                temp_df['Energy demand (GWh)'] = self.df['annual_el_demand']/1000000
                summary = summary.join(temp_df.groupby(geo_boundary)['Energy demand (GWh)'].sum())
            except:
                pass
            
            ### lcoe related data
            try:
                temp_df['Average lcoe ($/kWh)'] = self.df['lcoe']
                temp_df['Pumping cost (M$)'] = self.df['pumping_cost']/1000000
                summary = summary.join(temp_df.groupby(geo_boundary)['Average lcoe ($/kWh)'].mean())
                summary = summary.join(temp_df.groupby(geo_boundary)['Pumping cost (M$)'].sum())
                summary['Pumping cost ($/m3)'] = summary['Pumping cost (M$)'] / \
                                                      summary['Water demand (Mm3)']
            except:
                pass
        
        summary.round(decimals=3)
        return summary
            
    class Technology():
        df = pd.DataFrame()
        max_cap = pd.DataFrame()
        lcoe = pd.DataFrame()
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
            
            
            
            
            
            
            
            