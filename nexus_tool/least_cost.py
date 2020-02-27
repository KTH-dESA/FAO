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
    
def get_wind_cf(df, wind, mu, t, p_rated, z, zr, es, u_arr, p_curve, axis=1):
    cf_df = pd.DataFrame()
    if axis:
        for i in range (1,13):
            _wind = f'{wind}{i}'
            cf_df['cf_{}'.format(i)] = wind_cf(df, _wind, mu, t, p_rated, z, 
                                                 zr, es, u_arr, p_curve)
    else:
        cf_df['Year'] = df.Year
        cf_df['Month'] = df.Month
        cf_df['cf'] = wind_cf(df, wind, mu, t, p_rated, z, 
                                                 zr, es, u_arr, p_curve)
    return cf_df
    
def get_pv_cf(df, srad, axis=1):
    cf_df = pd.DataFrame()
    if axis:
        for i in range (1,13):
            cf_df['cf_{}'.format(i)] = df[f'{srad}{i}'] / (60*60*24) # solar rad: (kJ/(m2.day))*30.day/month*1h/(60*60s) = kWh/(m2.month)*30/(60*60)
                                                                            # pv_cf: energy produced / (p_rated * t)
    else:
        cf_df['Year'] = df.Year
        cf_df['Month'] = df.Month
        cf_df['cf'] = df[srad] / (60*60*24)
    return cf_df
    
def get_installed_capacity(df, cf, pd_e, axis=1):
    ic_df = pd.DataFrame()
    if axis:
        for i in range(1,13):
            _cf = cf if type(cf) == float else cf[f'cf_{i}']
            ic_df[f'ic_{i}'] = df[f'{pd_e}{i}'] / _cf
    else:
        ic_df['Year'] = df.Year
        ic_df['Month'] = df.Month
        ic_df['Demand point'] = df['Demand point']
        cf = cf if type(cf) == float else cf.cf
        ic_df['ic'] = df[pd_e] / cf
    return ic_df
    
def get_max_capacity(df, axis=1):
    if axis:
        return pd.DataFrame({'max_cap': df.filter(like='ic_').max(axis=1)})
    else:
        return df[['Demand point', 'Year', 'ic']].groupby(['Demand point', 'Year']).max()


def get_fuel_cost(fuel_cost, el_gen, efficiency, fuel_req):
    if type(fuel_cost) == int:
        fuel = el_gen * fuel_req * fuel_cost / efficiency
        fuel[0] = 0
    else:
        fuel = np.array([gen * fuel_req * cost / efficiency for gen, cost in zip(el_gen, fuel_cost)])
    return fuel
    
def get_emissions(el_gen, efficiency, fuel_req, emission_factor):
    emissions = el_gen * fuel_req * emission_factor / efficiency
    return emissions

def get_lcoe(max_capacity, total_demand, tech_life, om_cost, capital_cost,
             discount_rate, project_life, fuel_cost, fuel_req, 
             efficiency, emission_factor, env_cost, start_year, end_year, axis=1):
    if axis:
        # Perform the time-value LCOE calculation
        reinvest_year = 0
        capital_cost = capital_cost * max_capacity
        om_cost = om_cost * capital_cost

        # If the technology life is less than the project life, we will have to invest twice to buy it again
        if tech_life < project_life:
            reinvest_year = tech_life

        year = np.arange(project_life)
        el_gen = np.array([demand * np.ones(project_life -1) for demand in total_demand])
        el_gen = np.insert(el_gen,0,0,axis=1)
        discount_factor = np.array([(1 + discount_rate) ** year for i in total_demand])
        investments = np.array([np.zeros(project_life-1) for i in capital_cost])
        investments = np.insert(investments,0,capital_cost,axis=1)
        
        if reinvest_year:
            investments = np.delete(investments,reinvest_year,axis=1)
            investments = np.insert(investments,reinvest_year,capital_cost,axis=1)

        salvage = np.array([np.zeros(project_life-1) for i in capital_cost])
        used_life = project_life
        if reinvest_year:
            # salvage will come from the remaining life after the re-investment
            used_life = project_life - tech_life
        salvage = np.insert(salvage,-1,capital_cost * (1 - used_life / tech_life),axis=1)

        operation_and_maintenance = np.array([i * np.ones(project_life-1) for i in om_cost])
        operation_and_maintenance = np.insert(operation_and_maintenance,0,0,axis=1)
        
        fuel = get_fuel_cost(fuel_cost,el_gen,efficiency,fuel_req)
        emissions = get_emissions(el_gen,efficiency,fuel_req,emission_factor)
        
        discounted_costs = (investments + operation_and_maintenance + fuel + 
                            emissions * env_cost - salvage) / discount_factor
        discounted_generation = el_gen / discount_factor
        
        return discounted_costs.sum(axis=1) / discounted_generation.sum(axis=1)
    else:
        reinvest_year = 0
        df = get_capital_cost(max_capacity, start_year, end_year, 
                              tech_life, capital_cost)
        df = df.loc[(df.Year>=(start_year))&(df.Year<=(end_year))]
        data = total_demand.copy()
        total_demand = data.groupby(['Demand point', 'Year']).swpa_e.sum()
        df['total_demand'] = total_demand.reset_index().swpa_e
        if type(capital_cost)==pd.core.series.Series:
            capital_cost = capital_cost[start_year]
        df['om_cost'] = om_cost * capital_cost
        df['discount_factor'] = (1 + discount_rate) ** (df.Year - start_year)
        df = get_salvage(df, start_year, end_year, tech_life)
        df['emissions'] = get_emissions(df['total_demand'], efficiency, 
                                        fuel_req, emission_factor)
        if type(fuel_cost)==pd.core.series.Series:
            fuel_cost = df.Year.map(fuel_cost)
        df['fuel_cost'] = df['total_demand'] * fuel_req * fuel_cost / efficiency
        df['discounted_costs'] = (df['capital_cost'] + df['om_cost'] + 
                                  df['fuel_cost'] + df['emissions'] * env_cost - 
                                  df['salvage']) / df['discount_factor']
        df['discounted_generation'] = df['total_demand'] / df['discount_factor']
        dff = df.groupby('Demand point')[['discounted_costs', 
                                          'discounted_generation']].sum()
        dff.reset_index(inplace=True)
        dff['lcoe'] = dff['discounted_costs'] / dff['discounted_generation']
        dff.loc[dff.lcoe==np.inf, 'lcoe'] = np.nan
        dff['year'] = start_year
        
        return dff
        
def get_salvage(df, start_year, end_year, tech_life):
    year_used = ((end_year-start_year) % tech_life)
    years_left = tech_life - year_used - 1
    df['salvage'] = 0
    df.loc[df.Year==end_year, 'salvage'] = np.array(df.loc[(df.inv_period==df.inv_period.max())&(df.Year==(end_year-year_used)), 'capital_cost']  * (years_left/tech_life))
    return df
        
def get_capital_cost(max_capacity, start_year, end_year, tech_life, capital_cost):
    df = max_capacity.copy()
    df.loc[(df.Year<=start_year)&(df.Year<=start_year)]
    df['inv_period'] = (df.Year - start_year + tech_life)//tech_life
    df['invest_year'] = (df.Year - start_year + tech_life)%tech_life == 0
    dff = df.groupby(['Demand point', 'inv_period']).ic.max()
    dff = dff.reset_index()
    dff['invest_year'] = True
    dff.set_index(['Demand point', 'inv_period', 'invest_year'], inplace=True)
    df['new_capacity'] = df.set_index(['Demand point', 'inv_period', 'invest_year']).index.map(dff.ic)
    df['new_capacity'] = df['new_capacity'].fillna(0)
    if type(capital_cost)==pd.core.series.Series:
        capital_cost = df.Year.map(capital_cost)
    df['capital_cost'] = capital_cost * df.new_capacity
    return df
   
def get_least_cost(df, geo_boundary_col = None, geo_boundary_name = None):
    if geo_boundary_col == None:
        filter_vec = [True] * df.shape[0]
        i = 0
    else:
        filter_vec = df[geo_boundary_col]==geo_boundary_name
        i = 1
    df.loc[filter_vec, 'least_cost_technology'] = \
                                df.loc[filter_vec, list(df)[i:]].idxmin(axis=1)
    df.loc[filter_vec, 'lcoe'] = df.loc[filter_vec, list(df)[i:]].min(axis=1)
    return df.loc[filter_vec, ['least_cost_technology', 'lcoe']]

def get_tech_generation(df, technologies):
    for key in technologies:
                df.loc[df['least_cost_tech']==key, f'{key} generation'] = \
                        df.loc[df['least_cost_tech']==key, 'annual_el_demand']
    return df

def get_pumping_cost(df, energy_demand, lcoe):
    df['pumping_cost'] = df[energy_demand] * df[lcoe]
    return df

def get_unit_pumping_cost(df, pumping_cost, water_demand):
    df['unit_pumping_cost'] = df[pumping_cost] / water_demand
    return df