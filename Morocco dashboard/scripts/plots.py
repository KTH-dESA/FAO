import pandas as pd
import plotly.express as px
import os

def load_results_data(working_directory, scenarios, climates):
    df_results = pd.DataFrame()
    df_desal = pd.DataFrame()
    df_wwtp = pd.DataFrame()
    for scenario in scenarios:
        for climate in climates:
            folder = os.path.join(working_directory, 'data', scenario, climate)
            # read general results
            dff = pd.read_csv(os.path.join(folder, 'results.gz'))
            dff['Scenario'] = scenario.split('Wastewater Reuse')[0].strip()
            dff['Reuse'] = 'Yes' if 'Reuse' in scenario else 'No'
            df_results = df_results.append(dff, ignore_index=True)
            # read desal results
            dff = pd.read_csv(os.path.join(folder, 'desal_data.gz'))
            dff['Scenario'] = scenario.split('Wastewater Reuse')[0].strip()
            dff['Reuse'] = 'Yes' if 'Reuse' in scenario else 'No'
            df_desal = df_desal.append(dff, ignore_index=True)
            # read wastewater treatment results
            dff = pd.read_csv(os.path.join(folder, 'wwtp_data.gz'))
            dff['Scenario'] = scenario.split('Wastewater Reuse')[0].strip()
            dff['Reuse'] = 'Yes' if 'Reuse' in scenario else 'No'
            df_wwtp = df_wwtp.append(dff, ignore_index=True)
    return df_results, df_desal, df_wwtp

def water_supply_compare_plot(water_delivered, time_frame, title):
    dff_delivered = water_delivered.copy()
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    dff_delivered.loc[dff_delivered['type'].str.contains('WWR'), 'type'] = 'Reused Wastewater'
    dff_delivered = dff_delivered.loc[dff_delivered['type'] != 'Transmission Pipeline']
    dff_delivered = dff_delivered.groupby(['Scenario', time_frame, 'type', 'Reuse'])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()

    df = dff_delivered.reset_index()
    
    for t in df['type'].unique():
        if df.loc[df['type']==t, 'sswd'].sum() == 0:
            df = df.loc[df['type']!=t]
    fig = px.area(df, x='Year', y='sswd', color='type', facet_col='Scenario', facet_row='Reuse',
                  labels={"sswd": "Water (Mm<sup>3</sup>)"}, title=title,
                  facet_col_spacing=0.06, color_discrete_sequence=px.colors.qualitative.Dark2)
    return fig
    
def energy_demand_compare_plot(water_delivered, wwtp_data, desal_data, time_frame, title):
    dff_energy = water_delivered.copy()
    dff_energy.loc[dff_energy['type'].str.contains('GW'), 'type'] = 'Groundwater pumping'
    dff_energy.loc[dff_energy['type'].str.contains('SW|Pipeline'), 'type'] = 'Surface water conveyance'
    dff_energy.loc[dff_energy['type'].str.contains('DS'), 'type'] = 'Desalinated water conveyance'
    dff_energy.loc[dff_energy['type'].str.contains('WWR'), 'type'] = 'Wastewater reuse conveyance'
    dff_energy = dff_energy.groupby(['Scenario', time_frame, 'type', 'Reuse'])['swpa_e'].sum() / 1000000
    dff_energy = dff_energy.reset_index()
    wwtp_data['type'] = 'Wastewater treatment'
    desal_data['type'] = 'Desalination energy'
    for df in [wwtp_data, desal_data]:
        dff = df.groupby(['Scenario', time_frame, 'type', 'Reuse'])['swpa_e'].sum() / 1000000
        dff = dff.reset_index()
        dff_energy = dff_energy.append(dff, sort=False)
    for t in dff_energy['type'].unique():
        if dff_energy.loc[dff_energy['type']==t, 'swpa_e'].sum() == 0:
            dff_energy = dff_energy.loc[dff_energy['type']!=t]
            
    fig = px.area(dff_energy, x='Year', y='swpa_e', color='type', facet_col='Scenario', facet_row='Reuse',
                  labels={"swpa_e": "Energy (GWh)"}, title=title,
                  facet_col_spacing=0.06, color_discrete_sequence=px.colors.qualitative.T10)
    return fig
    
def unmet_demand_compare_plot(water_delivered, time_frame, title):
    dff_unmet = water_delivered.copy()
    dff_unmet.loc[dff_unmet['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_unmet.loc[dff_unmet['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_unmet = dff_unmet.loc[~dff_unmet['type'].str.contains('Aquifer')]
    dff_unmet = dff_unmet.loc[dff_unmet['type'] != 'Transmission Pipeline']
    water_req_year = \
        dff_unmet.groupby(['Scenario','Year', 'Date', 'Demand point', 'category', 'Reuse'])['water_required'].mean().reset_index().groupby(
            ['Scenario', time_frame, 'category', 'Reuse'])['water_required'].sum()
    unment_demand = pd.DataFrame()
    unment_demand['value'] = 1 - (dff_unmet.groupby(['Scenario', time_frame, 'category', 'Reuse'])['sswd'].sum() /
                         water_req_year)

    df = unment_demand.reset_index()
    df = df.loc[df['category'] == 'Agriculture']
    reuse = {'Yes': ' with reuse', 'No': '', }
    df['FullScenario'] = df['Scenario'] + [reuse[i] for i in df['Reuse']]

    fig = px.line(df, x='Year', y='value', color='FullScenario',  # facet_col='Reuse', #facet_row='Reuse',
                        labels={"value": "Unmet demand (%)", 'FullScenario': 'Scenario'},
                        title=title,
                        facet_col_spacing=0.06, color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_layout(yaxis={'tickformat': "%"})
    return fig
    
def load_butane_results(working_directory, scenarios, pv_levels):
    all_dfs = pd.DataFrame()
    for scenario in scenarios:
        phaseout = f'by {scenario}' if scenario else 'None'
        for pv_level in pv_levels:
            file = os.path.join(working_directory, 'Butane_calculations', f'phaseout_{scenario}', f'{pv_level}_PV', 'summary_results.gz')
            dff = pd.read_csv(file)
            dff['butane_phaseout'] = phaseout
            dff['pv_adoption'] = f'{pv_level}%'
            
            all_dfs = all_dfs.append(dff, ignore_index=True)
    
    return all_dfs
    

def total_costs_plot(df, title):    
    dff = df.melt(value_vars=['butane_SUBSIDY(mMAD)', 'grid_cost(mMAD)', 'PV_Capex(mMAD)'], 
             id_vars=['butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable']=="butane_SUBSIDY(mMAD)", 'variable'] = "Butane"
    dff.loc[dff['variable']=="grid_cost(mMAD)", 'variable'] = "Grid"
    dff.loc[dff['variable']=="PV_Capex(mMAD)", 'variable'] = "PV"

    fig = px.bar(dff, x='butane_phaseout', y='value', color='variable',
                 facet_col='pv_adoption', title=title,
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 labels={
                         "butane_phaseout": "Butane phase-out",
                         "pv_adoption": "PV adoption by 2040",
                         "value": "Total costs (mMAD)",
                         "variable": 'Source'
                         }
                )
    fig.update_layout(xaxis={'categoryorder':'total descending'})
    return fig
    
def emissions_compare_plot(df, title):
    df['Total_emissions(MtCO2)'] = df[['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)']].sum(axis=1)
    dff = df.groupby(['Year','butane_phaseout', 'pv_adoption'])['Total_emissions(MtCO2)'].sum().reset_index()

    dff = dff.copy().melt(value_vars=['Total_emissions(MtCO2)'], 
             id_vars=['Year','butane_phaseout', 'pv_adoption'])

    fig = px.line(dff, x='Year', y='value', color='butane_phaseout',
                 facet_col='pv_adoption', title=title,
                 color_discrete_sequence=px.colors.qualitative.T10,
                 labels={
                         "butane_phaseout": "Butane phase-out",
                         "pv_adoption": "PV adoption by 2040",
                         "value": "Emissions (MtCO2)"
                         },
                 facet_col_spacing=0.06,
                 )
    return fig
    
def total_emissions_compare_plot(df, title):
    dff = df.melt(value_vars=['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)'], 
                  id_vars=['butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable']=="butane_emissions(MtCO2)", 'variable'] = "Butane"
    dff.loc[dff['variable']=="grid_emissions(MtCO2)", 'variable'] = "Grid"

    fig = px.bar(dff, x='butane_phaseout', y='value', color='variable',
                 facet_col='pv_adoption', title=title,
                 color_discrete_sequence=px.colors.qualitative.Vivid,
                 labels={
                         "butane_phaseout": "Butane phase-out",
                         "pv_adoption": "PV adoption by 2040",
                         "value": "Emissions (MtCO2)",
                         "variable": 'Source'
                         }
                )

    fig.update_layout(xaxis={'categoryorder':'array', 'categoryarray':['None','by 2040','by 2030']})
    return fig
    
def emisions_vs_costs(df, title):
    df = df.loc[(df.Year >= 2021) & (df.Year < 2050)].copy()

    df['Total_emissions(MtCO2)'] = df[['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)']].sum(axis=1)
    df['Total_costs(mMAD)'] = df[['butane_SUBSIDY(mMAD)', 'grid_cost(mMAD)', 'PV_Capex(mMAD)']].sum(axis=1)

    dff = pd.DataFrame()
    for key, group in df.groupby(['butane_phaseout', 'pv_adoption']):
        _dff = group.copy()
        
        _dff['Cumulative emissions (MtCO2)'] = _dff['Total_emissions(MtCO2)'].cumsum()
        _dff['Cumulative costs (mMAD)'] = _dff['Total_costs(mMAD)'].cumsum()
        dff = dff.append(_dff, ignore_index=True)

    pv_adop_dic = {'10%': 0.1, '20%': 0.2, '50%': 0.5}

    dff = dff.groupby(['butane_phaseout', 'pv_adoption']).sum().reset_index()

    dff['pv_adoption_number'] = [pv_adop_dic[i] for i in dff['pv_adoption']]

    fig = px.scatter(dff, x='Total_costs(mMAD)', y='Total_emissions(MtCO2)',
                     title=title, 
                     color='butane_phaseout', size='pv_adoption_number',
                     color_discrete_sequence=px.colors.qualitative.Vivid,
                     labels={
                             "butane_phaseout": "Butane phase-out",
                             "pv_adoption": "PV adoption share by 2040",
                             "Total_emissions(MtCO2)": "Emissions (MtCO2)",
                             "Total_costs(mMAD)": 'Total costs (mMAD)',
                             "pv_adoption_number": 'PV share by 2040'
                             }
                    )
    return fig
    
def energy_resources_share_plot(df, title): 
    dff = df.groupby(['Year','butane_phaseout', 'pv_adoption'])[['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum().reset_index()

    dff['pv_demand(%)'] = dff['pv_demand(KWh)'] / dff[['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(axis=1)
    dff['butane_demand(%)'] = dff['butane_demand(KWh)'] / dff[['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(axis=1)
    dff['grid_demand(%)'] = dff['grid_demand(KWh)'] / dff[['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(axis=1)

    dff = dff.melt(value_vars=['pv_demand(%)', 'butane_demand(%)', 'grid_demand(%)'], 
             id_vars=['Year','butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable']=="butane_demand(%)", 'variable'] = "Butane"
    dff.loc[dff['variable']=="grid_demand(%)", 'variable'] = "Grid"
    dff.loc[dff['variable']=="pv_demand(%)", 'variable'] = "PV"

    fig = px.bar(dff, x='Year', y='value', color='variable',
                 facet_col='pv_adoption', facet_row='butane_phaseout', title=title,
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 labels={
                         "butane_phaseout": "phaseout",
                         "pv_adoption": "PV adoption by 2040",
                         "value": "Energy share",
                         "variable": 'Source'
                         },
                 category_orders={'butane_phaseout': ['None', 'by 2040', 'by 2030']},
                 facet_row_spacing=0.04,
                 facet_col_spacing=0.04,
                )
    fig.update_layout(height=600, bargap=0)
    return fig
    
def energy_demand_ag(df):
    df = df.melt(value_vars=['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)'], 
                   id_vars=['Year'])
    df['value'] /= 1000000

    df.loc[df['variable']=='butane_demand(KWh)', 'variable'] = "Butane"
    df.loc[df['variable']=='grid_demand(KWh)', 'variable'] = "Grid"
    df.loc[df['variable']=='pv_demand(KWh)', 'variable'] = "PV"

    fig = px.bar(df, x='Year',  y='value', color='variable',
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 labels={'value': 'Energy demand (GWh)',
                         'variable': 'Source'}
                )
    return fig
    
def pv_installed_capacity(df):
    df = df
    df.loc[df.Year==2020,'PV_new_cap(MW)'] = 0
    df['previous_capacity(MW)'] = df['cap_a(MW)'] - df['PV_new_cap(MW)']
    df = df.melt(value_vars=['previous_capacity(MW)', 'PV_new_cap(MW)'], id_vars=['Year'])
    
    df.loc[df['variable']=='previous_capacity(MW)', 'variable'] = 'Previous capacity'
    df.loc[df['variable']=='PV_new_cap(MW)', 'variable'] = 'New capacity'

    fig = px.bar(df, x='Year', y='value', color='variable',
                 color_discrete_sequence=px.colors.qualitative.Pastel,
                )
    return fig

def emissions_ag(df):
    df = df.melt(value_vars=['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)'], 
                       id_vars=['Year'])

    df.loc[df['variable']=='butane_emissions(MtCO2)', 'variable'] = "Butane"
    df.loc[df['variable']=='grid_emissions(MtCO2)', 'variable'] = "Grid"

    fig = px.bar(df, x='Year', y='value', color='variable',
                 title='Yearly emissions',
                 color_discrete_sequence=px.colors.qualitative.Vivid,
                 labels={'value': 'Emissions (MtCO2)',
                         'variable': 'Source'}
                )
    return fig