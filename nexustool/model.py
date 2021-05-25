# Standard library imports
import pandas as pd
import numpy as np
from math import pi

# Local application/library specific imports
from .water_demand import (
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


from .least_cost import (
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
    elevation_diff = 'elevation'
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
    acwr = 'ACWR'
    pcwr = 'PCWR'
    pwd = 'PWD'
    sswd = 'SSWD'
    crop_calendar = None
    crop_column = 'crop'
    ky_dict = {}
    kc_dict = {}
    pipeline = 'pipeline'
    well = 'point'
    wwtp = 'point'
    desal_plant = 'point'
    pumping_hours_per_day = 'pumpng_hours_per_day'
    deff = 1
    aeff = 0.45
    ky_values = {}
    kc_values = {}
    # energy properties:
    gw_depth = 'gw_depth'
    D = 'Pipe_diameter'
    L = 'Pipeline_length'
    A = 'Pipe_area'
    mV = 'Flow_velocity'  # monthly Velocity mV in (m/sec)
    tdh = 'tdh'
    des_int = 'Einten_KWh/m3'
    des_ener = 'Edesal_GWh'
    pp_e = 'PP_E'  # Pumping PEAK Electric demand
    pa_e = 'PA_E'  # Pumping AVERAGE Electric demand
    Re = 'Re'
    f = 'f'
    trans_eff = 0
    pump_eff = 'pump_eff'
    technologies = {}
    discount_rate = 0
    g = 9.81  # gravitational acceleration in(m/sec2)
    Ken_visc = 1.004e-06
    dens = 1000
    k = 0.061e-3  # roughness for ductile iron pipe
    ww_energy_int = 'ww_energy_int'
    desal_energy_int = 'desal_energy_int'
    start_year = 0
    end_year = 30

    def __init__(self, df):
        self.df = df

    def print_properties(self):
        print('Properties names:')
        for val, name in zip([self.eto, self.lat, self.elevation_diff, self.wind,
                              self.srad, self.tmin, self.tmax, self.tavg,
                              self.crop_share, self.crop_area, self.seasons,
                              self.start, self.end, self.crop_column,
                              self.gw_depth, self.tdh],
                             ['Reference evapotranspiration (.eto)',
                              'Latitude (.lat)', 'Elevation difference (.elevation_diff)',
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
                              'Reynolds Number (.Re)', ]):
            print('    - {}: {}'.format(name, val))

            ####### water related methods ###########

    def set_cropland_share(self, crop_var, geo_boundary='global',
                           boundary_name=None, inplace=False):
        if inplace:
            set_cropland_share(self.df, crop_var, geo_boundary=geo_boundary,
                               boundary_name=boundary_name, crop_share=self.crop_share)
        else:
            return set_cropland_share(self.df.copy(), crop_var,
                                      geo_boundary=geo_boundary,
                                      boundary_name=boundary_name,
                                      crop_share=self.crop_share)

    def get_ky_list(self, inplace=False):
        if inplace:
            get_ky_list(self.df, crop_share=self.crop_share)
        else:
            return get_ky_list(self.df.copy(), crop_share=self.crop_share)

    def get_kc_list(self, inplace=False):
        if inplace:
            get_kc_list(self.df, crop_share=self.crop_share)
        else:
            return get_kc_list(self.df.copy(), crop_share=self.crop_share)

    def get_eto(self, inplace=False):
        if inplace:
            get_eto(self.df, eto=self.eto, lat=self.lat,
                    elevation=self.elevation_diff, wind=self.wind,
                    srad=self.srad, tmin=self.tmin,
                    tmax=self.tmax, tavg=self.tavg)
        else:
            return get_eto(self.df.copy(), eto=self.eto, lat=self.lat,
                           elevation=self.elevation_diff, wind=self.wind,
                           srad=self.srad, tmin=self.tmin,
                           tmax=self.tmax, tavg=self.tavg)

    def get_effective_rainfall(self, inplace=False):
        if inplace:
            get_effective_rainfall(self.df, eff=self.eff, prec=self.prec,
                                   eto=self.eto)
        else:
            return get_effective_rainfall(self.df.copy(), eff=self.eff,
                                          prec=self.prec, eto=self.eto)

    def get_calendar_days(self, inplace=False):
        if inplace:
            get_calendar_days(self.crop_calendar, seasons=self.seasons,
                              start=self.start, end=self.end)
        else:
            return get_calendar_days(self.crop_calendar.copy(), seasons=self.seasons,
                                     start=self.start, end=self.end)

    def get_kc_values(self, inplace=False):
        if inplace:
            get_kc_values(crop_calendar=self.crop_calendar,
                          seasons=self.seasons, kc_dict=self.kc_dict,
                          crop_column=self.crop_column, start=self.start,
                          end=self.end, kc=self.kc)
        else:
            return get_kc_values(crop_calendar=self.crop_calendar.copy(),
                                 seasons=self.seasons, kc_dict=self.kc_dict,
                                 crop_column=self.crop_column,
                                 start=self.start, end=self.end,
                                 kc=self.kc)

    def get_water_demand(self, inplace=False):
        if inplace:
            get_water_demand(self.df, self.crop_calendar, self.ky_dict,
                             self.crop_column, self.aeff, self.deff,
                             self.seasons[0], self.seasons[3],
                             self.pumping_hours_per_day,
                             crop_area=self.crop_area, _eto=self.eto,
                             _kc=self.kc, _eff=self.eff, _acwr=self.acwr,
                             _pcwr=self.pcwr, _pwd=self.pwd, _sswd=self.sswd,
                             start=self.start, end=self.end,
                             crop_share=self.crop_share)
        else:
            return get_water_demand(self.df.copy(), self.crop_calendar,
                                    self.ky_dict, self.crop_column, self.aeff,
                                    self.deff, self.seasons[0], self.seasons[3],
                                    self.pumping_hours_per_day,
                                    crop_area=self.crop_area, _eto=self.eto,
                                    _kc=self.kc, _eff=self.eff, _acwr=self.acwr,
                                    _pcwr=self.pcwr, _pwd=self.pwd,
                                    _sswd=self.sswd, start=self.start,
                                    end=self.end, crop_share=self.crop_share)

    ####### energy related methods ########### 
    def set_efficiency_goal(self, dff, initial_eff, final_eff, 
                            init_year, end_year):
        years = np.arange(year2 - year1 + 1) + 2020
        efficiencies = np.linspace(initial_eff, final_eff, 
                                   end_year - init_year + 1)
        df_eff = pd.DataFrame({'Year': years, self.pump_eff: efficiencies})
    return dff.merge(df_eff, on='Year')[self.pump_eff]
    
    
    def set_specifications(self, category, diameter, pumping_hours_per_day=None, 
                           pump_efficiency=1, efficiency_goal=None, amount=None, 
                           distribute=None, names=None):
        '''Sets the specification for pipelines or wells

        Parameters
        ----------
        category : ['pipeline', 'well']
            The category to configure
        diameter : float, int
            The diameter of the pipeline or well
        pumping_hours_per_day : float, int
            The average amount of hours pumped per day
        pump_efficiency : float
            Pumping energy efficiency
        amount : int
            Amount of pipelines or wells
        distribute : [None, 'proportionally']
            If 'proportionally, distributes the amount of specified wells 
            proportionally based on the water extractions per well field
        names: string, list
            Contains the name or names of the pipelines/wells to configure

        Raises
        ------
        TypeError
            Invalid type for names
        '''
        if category == 'pipeline':
            feature = self.pipeline
        elif category == 'well':
            feature = self.well
        
        if isinstance(names, type(None)):
            self.df[self.D] = diameter #in m
            self.df[self.pumping_hours_per_day] = pumping_hours_per_day
            if efficiency_goal:
                for index, group in self.df.groupby(['Year', feature]):
                    group[self.pump_eff] = set_efficiency_goal(group, 
                                                               pump_efficiency, 
                                                               efficiency_goal['efficiency'], 
                                                               efficiency_goal['init_year'], 
                                                               efficiency_goal['end_year'])
            else:
                self.df[self.pump_eff] = pump_efficiency
            if distribute=='proportionally':
                dff = self.df.groupby(['Year', feature])[['sswd']].sum().reset_index()
                df_amount = dff.loc[dff['Year']==2020].copy()
                df_amount['pipe_amount'] = np.ceil(df_amount['sswd'] / df_amount['sswd'].sum() * amount)
                self.df = self.df.merge(df_amount[[feature, 'pipe_amount']], on=feature)
            else:
                self.df['pipe_amount'] = amount
        elif isinstance(names, list):
            for name in names:
                self.df.loc[self.df[feature]==name, self.D] = diameter #in m
                self.df.loc[self.df[feature]==name, self.pumping_hours_per_day] = pumping_hours_per_day
                self.df.loc[self.df[feature]==name, self.pump_eff] = pump_efficiency
                if amount:
                    self.df.loc[self.df[feature]==name, 'pipe_amount'] = amount
        elif isinstance(names, str):
            self.df.loc[self.df[feature]==names, self.D] = diameter #in m
            self.df.loc[self.df[feature]==names, self.pumping_hours_per_day] = pumping_hours_per_day
            self.df.loc[self.df[feature]==names, self.pump_eff] = pump_efficiency
            if amount:
                self.df.loc[self.df[feature]==names, 'pipe_amount'] = amount
        else:
            raise TypeError(f'Invalid type {type(names)} for names. Names should be a string or a list of strings')
    
    def pipe_area(self):
        self.df[self.A] = (pi*self.df[self.D]**2)/4

    def flow_velocity(self, axis=0):
        if axis:
            for i in range(1, 13):
                _avg_Q = '{}{}'.format(self.sswd, i)
                _mV = '{}{}'.format(self.mV, i)
                self.df[_mV] = self.df[_avg_Q] / (30 * self.df[self.pumping_hours_per_day] * 60 * 60) / self.df[self.A]
        else:
            self.df[self.mV] = self.df[self.sswd] / self.df['pipe_amount'] / (30 * self.df[self.pumping_hours_per_day] * 60 * 60) / self.df[self.A]  # convert m3/month to m3/s (30day*pump_hours_day*60min/h*60s/min)

    def reynolds(self, axis=0):
        if axis:
            for i in range(1, 13):
                _mV = '{}{}'.format(self.mV, i)
                _Re = '{}{}'.format(self.Re, i)
                self.df[_Re] = (self.df[_mV] * self.df[self.D]) / self.Ken_visc
        else:
            self.df[self.Re] = self.df[self.mV] * self.df[self.D] / self.Ken_visc

    def friction_factor(self, axis=0):
        if axis:
            for i in range(1, 13):
                _Re = '{}{}'.format(self.Re, i)
                _f = '{}{}'.format(self.f, i)
                self.df[_f] = 0.25 / (np.log((self.k / (3.7 * self.df[self.D])) + (5.74 / (self.df[_Re] ** 0.9))) ** 2)
        else:
            self.df[self.f] = 0.25 / (np.log((self.k / (3.7 * self.df[self.D])) + (5.74 / (self.df[self.Re] ** 0.9))) ** 2)

    def get_tdh(self, friction=True, axis=0):
        if axis:
            for i in range(1, 13):
                _avg_Q = '{}{}'.format(self.sswd, i)
                _f = '{}{}'.format(self.f, i)
                _tdh_sw = '{}{}'.format(self.tdh, i)

                if friction:
                    f_losses = (self.df[_f] * self.df[self.L] * 16 * ((self.df[_avg_Q] / (30 * self.df[self.pumping_hours_per_day] * 60 * 60)) ** 2)) / ((self.df[self.D] ** 5) * 2 * self.g * (pi ** 2))
                else:
                    f_losses = 0
                self.df[_tdh_sw] = self.df[self.elevation_diff] + f_losses
                self.df[_tdh_sw].replace(0, np.nan, inplace=True)
        else:
            if friction:
                f_losses = (self.df[self.f] * self.df[self.L] * (self.df[self.mV] ** 2)) / (self.df[self.D] * 2 * self.g)
            else:
                f_losses = 0
            self.df[self.tdh] = self.df[self.elevation_diff] + f_losses

    def get_pumping_energy(self, axis=0):
        if axis:
            for i in range(1, 13):
                _swpp_e = '{}{}'.format(self.pp_e, i)  # surface water pumping peak electric demand
                _swpa_e = '{}{}'.format(self.pa_e, i)  # surface water pumping average electric demand
                _peak_Q = '{}{}'.format(self.pwd, i)  # peak water flow in the pipeline. To be updated WEAP output.
                _avg_Q = '{}{}'.format(self.sswd, i)  # average water flow in the pipeline. To be updated with WEAP output
                _tdh_sw = '{}{}'.format(self.tdh, i)

                self.df[_swpp_e] = ((self.df[_peak_Q] * self.df[_tdh_sw] * self.g * self.dens) / (self.df[self.pump_eff] * 1000))  # to convert E from W to KW
                self.df[_swpa_e] = ((self.df[_avg_Q] * self.df[_tdh_sw] * self.g * self.dens) / (
                            self.df[self.pump_eff] * 1000 * 3600))  # to convert E from J to KWh

        else:
            if self.pwd in self.df.columns:
                self.df[self.pp_e] = ((self.df[self.pwd] * self.df[self.tdh] * self.g * self.dens) / (self.df[self.pump_eff] * 1000))
            self.df[self.pa_e] = ((self.df[self.sswd] * self.df[self.tdh] * self.g * self.dens) / (self.df[self.pump_eff] * 1000 * 3600))
            self.df.loc[self.df[self.pa_e]<0, self.pa_e] = 0
            
    def set_treatment_energy(self, category, energy_int, names=None):
        '''Sets the specification for pipelines or wells

        Parameters
        ----------
        category : ['wastewater', 'desalination']
            The category to configure
        energy_int : float
            Energy intensity of treatment
        names: string, list
            Contains the name or names of the treatment plants to configure

        Raises
        ------
        TypeError
            Invalid type for names
        '''
        if category == 'wastewater':
            feature = self.wwtp
            energy_feature = self.ww_energy_int
        elif category == 'desalination':
            feature = self.desal_plant
            energy_feature = self.desal_energy_int
        
        if isinstance(names, type(None)):
            self.df[energy_feature] = energy_int
        elif isinstance(names, list):
            for name in names:
                self.df.loc[self.df[feature]==name, energy_feature] = energy_int #in m
        elif isinstance(names, str):
            self.df.loc[self.df[feature]==names, energy_feature] = energy_int #in m
        else:
            raise TypeError(f'Invalid type {type(names)} for names. Names should be a string or a list of strings')

    def get_treatment_energy(self, category):
        if category == 'wastewater':
            energy_feature = self.ww_energy_int
        elif category == 'desalination':
            energy_feature = self.desal_energy_int
            
        self.df[self.pa_e] = self.df[self.sswd] * self.df[energy_feature]

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

    def get_cf(self, technologies='all', axis=0):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            if type(self.technologies[technology]) == self.WindTurbine:
                self.get_wind_cf(technology, axis)
            elif type(self.technologies[technology]) == self.PVSystem:
                self.get_pv_cf(technology, axis)

    def get_wind_cf(self, wind_turbine, axis=0):
        tech = self.technologies[wind_turbine]
        self.technologies[wind_turbine].cf = get_wind_cf(self.df, wind=self.wind,
                                                         mu=tech.mu, t=tech.t, p_rated=tech.p_rated,
                                                         z=tech.z, zr=tech.zr, es=tech.es, u_arr=tech.u_arr,
                                                         p_curve=tech.p_curve, axis=axis)

    def get_pv_cf(self, pv_system, axis=0):
        tech = self.technologies[pv_system]
        self.technologies[pv_system].cf = get_pv_cf(self.df, self.srad, axis)

    def get_installed_capacity(self, technologies='all', axis=0):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            self.technologies[technology].df = get_installed_capacity(self.df,
                                                                      tech.cf,
                                                                      self.pp_e,
                                                                      axis)

    def get_max_capacity(self, technologies='all', axis=0):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            if axis:
                self.technologies[technology].df = tech.df.join(get_max_capacity(tech.df, axis))
            else:
                self.technologies[technology].max_cap = get_max_capacity(tech.df, axis)

    def get_lcoe(self, technologies='all', years='all', axis=0):
        technologies = self.__check_tech_input(technologies)
        for technology in technologies:
            tech = self.technologies[technology]
            if axis:
                self.technologies[technology].df['lcoe'] = get_lcoe(
                    max_capacity=tech.df['max_cap'],
                    total_demand=self.df['annual_el_demand'],
                    tech_life=tech.life, om_cost=tech.om_cost,
                    capital_cost=tech.capital_cost,
                    discount_rate=self.discount_rate,
                    project_life=self.end_year - self.start_year,
                    fuel_cost=tech.fuel_cost,
                    fuel_req=tech.fuel_req,
                    efficiency=tech.efficiency,
                    emission_factor=tech.emission_factor,
                    env_cost=tech.env_cost,
                    start_year=self.start_year,
                    end_year=self.end_year,
                    axis=axis)
            else:
                years = self.__get_years(years)
                self.technologies[technology].lcoe = pd.DataFrame()
                for year in years:
                    self.technologies[technology].lcoe = \
                        self.technologies[technology].lcoe.append(get_lcoe(
                            max_capacity=tech.max_cap.reset_index(),
                            total_demand=self.df,
                            tech_life=tech.life, om_cost=tech.om_cost,
                            capital_cost=tech.capital_cost,
                            discount_rate=self.discount_rate,
                            project_life=self.end_year - self.start_year,
                            fuel_cost=tech.fuel_cost,
                            fuel_req=tech.fuel_req,
                            efficiency=tech.efficiency,
                            emission_factor=tech.emission_factor,
                            env_cost=tech.env_cost,
                            start_year=year,
                            end_year=self.end_year,
                            axis=axis), ignore_index=True)

                # self.technologies[technology].lcoe.reset_index(inplace=True)

    def get_least_cost(self, technologies='all', years='all',
                       geo_boundary=None, axis=0):
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
                    lcoe = get_least_cost(lcoe_df, geo_boundary, key)
                    self.df.loc[self.df[geo_boundary] == key, 'least_cost_tech'] = \
                        lcoe.loc[self.df[geo_boundary] == key, 'least_cost_technology']
                    self.df.loc[self.df[geo_boundary] == key, 'lcoe'] = \
                        lcoe.loc[self.df[geo_boundary] == key, 'lcoe']
            else:
                _technologies = self.__check_tech_input(technologies)
                lcoe_df = pd.DataFrame()
                for _technology in _technologies:
                    lcoe_df[_technology] = self.technologies[_technology].df['lcoe']
                lcoe = get_least_cost(lcoe_df)
                self.df['least_cost_tech'] = lcoe['least_cost_technology']
                self.df['lcoe'] = lcoe['lcoe']
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

    def get_pumping_cost(self, inplace=False):
        if inplace:
            get_pumping_cost(self.df, 'annual_el_demand', 'lcoe')
        else:
            return get_pumping_cost(self.df.copy(), 'annual_el_demand', 'lcoe')

    def get_unit_pumping_cost(self, inplace=False):
        if inplace:
            get_unit_pumping_cost(self.df, 'pumping_cost',
                                  self.df.filter(like=self.sswd).sum(axis=0))
        else:
            return get_unit_pumping_cost(self.df.copy(), 'pumping_cost',
                                         self.df.filter(like=self.sswd).sum(axis=0))

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
        t = 24 * 30
        p_rated = 600
        z = 55  # hub height
        zr = 80  # velocity measurement height
        es = 0.85  # losses in wind electricity
        u_arr = range(1, 26)
        p_curve = [0, 0, 0, 0, 30, 77, 135, 208, 287, 371, 450, 514, 558,
                   582, 594, 598, 600, 600, 600, 600, 600, 600, 600, 600, 600]

        def __init__(self, life, om_cost, capital_cost, efficiency, mu=mu,
                     t=t, p_rated=p_rated, z=z, zr=zr, es=es,
                     u_arr=u_arr, p_curve=p_curve):
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
