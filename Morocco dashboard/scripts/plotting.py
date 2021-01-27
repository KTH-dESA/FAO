import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


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
        if df.loc[df['type'] == t, 'sswd'].sum() == 0:
            df = df.loc[df['type'] != t]
    fig = px.area(df, x='Year', y='sswd', color='type', facet_col='Scenario', facet_row='Reuse',
                  labels={"sswd": "Water (Mm<sup>3</sup>)"}, title=title,
                  facet_col_spacing=0.06, color_discrete_sequence=px.colors.qualitative.Dark2)
    fig.update_layout(height=600)
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
        if dff_energy.loc[dff_energy['type'] == t, 'swpa_e'].sum() == 0:
            dff_energy = dff_energy.loc[dff_energy['type'] != t]

    fig = px.area(dff_energy, x='Year', y='swpa_e', color='type', facet_col='Scenario', facet_row='Reuse',
                  labels={"swpa_e": "Energy (GWh)"}, title=title,
                  facet_col_spacing=0.06, color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_layout(height=600)
    return fig


def unmet_demand_compare_plot(water_delivered, time_frame, title):
    dff_unmet = water_delivered.copy()
    dff_unmet.loc[dff_unmet['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_unmet.loc[dff_unmet['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_unmet = dff_unmet.loc[~dff_unmet['type'].str.contains('Aquifer')]
    dff_unmet = dff_unmet.loc[dff_unmet['type'] != 'Transmission Pipeline']
    water_req_year = \
        dff_unmet.groupby(['Scenario', 'Year', 'Date', 'Demand point', 'category', 'Reuse'])[
            'water_required'].mean().reset_index().groupby(
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
    fig.update_layout(height=600)
    return fig


def total_costs_plot(df, title):
    dff = df.melt(value_vars=['butane_SUBSIDY(mMAD)', 'grid_cost(mMAD)', 'PV_Capex(mMAD)'],
                  id_vars=['butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable'] == "butane_SUBSIDY(mMAD)", 'variable'] = "Butane"
    dff.loc[dff['variable'] == "grid_cost(mMAD)", 'variable'] = "Grid"
    dff.loc[dff['variable'] == "PV_Capex(mMAD)", 'variable'] = "PV"

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
    fig.update_layout(xaxis={'categoryorder': 'total descending'})
    # fig.update_layout(height=600)
    return fig


def emissions_compare_plot(df, title):
    df['Total_emissions(MtCO2)'] = df[['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)']].sum(axis=1)
    dff = df.groupby(['Year', 'butane_phaseout', 'pv_adoption'])['Total_emissions(MtCO2)'].sum().reset_index()

    dff = dff.copy().melt(value_vars=['Total_emissions(MtCO2)'],
                          id_vars=['Year', 'butane_phaseout', 'pv_adoption'])

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
    # fig.update_layout(height=600)
    return fig


def total_emissions_compare_plot(df, title):
    dff = df.melt(value_vars=['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)'],
                  id_vars=['butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable'] == "butane_emissions(MtCO2)", 'variable'] = "Butane"
    dff.loc[dff['variable'] == "grid_emissions(MtCO2)", 'variable'] = "Grid"

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

    fig.update_layout(xaxis={'categoryorder': 'array', 'categoryarray': ['None', 'by 2040', 'by 2030']})
    # fig.update_layout(height=600)
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
    # fig.update_layout(height=600)
    return fig


def energy_resources_share_plot(df, title):
    dff = df.groupby(['Year', 'butane_phaseout', 'pv_adoption'])[
        ['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum().reset_index()

    dff['pv_demand(%)'] = dff['pv_demand(KWh)'] / dff[['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(
        axis=1)
    dff['butane_demand(%)'] = dff['butane_demand(KWh)'] / dff[
        ['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(axis=1)
    dff['grid_demand(%)'] = dff['grid_demand(KWh)'] / dff[
        ['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)']].sum(axis=1)

    dff = dff.melt(value_vars=['pv_demand(%)', 'butane_demand(%)', 'grid_demand(%)'],
                   id_vars=['Year', 'butane_phaseout', 'pv_adoption'])

    dff.loc[dff['variable'] == "butane_demand(%)", 'variable'] = "Butane"
    dff.loc[dff['variable'] == "grid_demand(%)", 'variable'] = "Grid"
    dff.loc[dff['variable'] == "pv_demand(%)", 'variable'] = "PV"

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


def energy_demand_ag(df, layout):
    df = df.melt(value_vars=['pv_demand(KWh)', 'butane_demand(KWh)', 'grid_demand(KWh)'],
                 id_vars=['Year'])
    df['value'] /= 1000000

    df.loc[df['variable'] == 'butane_demand(KWh)', 'variable'] = "Butane"
    df.loc[df['variable'] == 'grid_demand(KWh)', 'variable'] = "Grid"
    df.loc[df['variable'] == 'pv_demand(KWh)', 'variable'] = "PV"

    fig = px.bar(df, x='Year', y='value', color='variable',
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 labels={'value': 'Energy demand (GWh)',
                         'variable': 'Source'}
                 )
    fig.update_layout(layout)
    return fig


def pv_installed_capacity(df, layout):
    df = df
    df.loc[df.Year == 2020, 'PV_new_cap(MW)'] = 0
    df['previous_capacity(MW)'] = df['cap_a(MW)'] - df['PV_new_cap(MW)']
    df = df.melt(value_vars=['previous_capacity(MW)', 'PV_new_cap(MW)'], id_vars=['Year'])

    df.loc[df['variable'] == 'previous_capacity(MW)', 'variable'] = 'Previous capacity'
    df.loc[df['variable'] == 'PV_new_cap(MW)', 'variable'] = 'New capacity'

    fig = px.bar(df, x='Year', y='value', color='variable',
                 color_discrete_sequence=px.colors.qualitative.Pastel,
                 )
    fig.update_layout(layout)
    return fig


def emissions_ag(df, layout):
    df = df.melt(value_vars=['butane_emissions(MtCO2)', 'grid_emissions(MtCO2)'],
                 id_vars=['Year'])

    df.loc[df['variable'] == 'butane_emissions(MtCO2)', 'variable'] = "Butane"
    df.loc[df['variable'] == 'grid_emissions(MtCO2)', 'variable'] = "Grid"

    fig = px.bar(df, x='Year', y='value', color='variable',
                 title='Yearly emissions',
                 color_discrete_sequence=px.colors.qualitative.Vivid,
                 labels={'value': 'Emissions (MtCO2)',
                         'variable': 'Source'}
                 )
    fig.update_layout(layout)
    return fig

def water_delivered_plot(water_delivered, time_frame, layout):
    dff_delivered = water_delivered.copy()
    dff_delivered.loc[dff_delivered['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_delivered.loc[dff_delivered['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_delivered.loc[dff_delivered['type'].str.contains('Aquifer'), 'category'] = 'Aquifer recharge'
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    dff_delivered.loc[dff_delivered['type'].str.contains('WWR'), 'type'] = 'Wastewater reuse'
    dff_delivered = dff_delivered.loc[dff_delivered['type'] != 'Transmission Pipeline']
    dff_delivered = dff_delivered.groupby([time_frame, 'category', 'type'])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()

    dff_delivered_cat = dff_delivered.groupby([time_frame, 'category']).sswd.sum()
    dff_delivered['share'] = dff_delivered.sswd / dff_delivered.set_index([time_frame, 'category']).index.map(
        dff_delivered_cat)

    df = dff_delivered_cat.reset_index()
    fig = px.area(df, x=time_frame, y='sswd', color='category',
                  color_discrete_sequence=px.colors.qualitative.Dark2)
    fig.update_traces(fill='tonexty')
    fig.update_layout(layout)
    return fig


def unmet_demand_plot(water_delivered, time_frame, layout):
    dff_unmet = water_delivered.copy()
    dff_unmet.loc[dff_unmet['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_unmet.loc[dff_unmet['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_unmet = dff_unmet.loc[~dff_unmet['type'].str.contains('Aquifer')]
    dff_unmet = dff_unmet.loc[dff_unmet['type'] != 'Transmission Pipeline']
    water_req_year = \
        dff_unmet.groupby(['Year', 'Date', 'Demand point', 'category'])['water_required'].mean().reset_index().groupby(
            [time_frame, 'category'])['water_required'].sum()
    unmet_demand = 1 - (dff_unmet.groupby([time_frame, 'category'])['sswd'].sum() /
                         water_req_year)

    df = unmet_demand.reset_index().rename(columns={0: 'unmet_demand'})
    # df.loc[df['unmet_demand']<0.001, 'unmet_demand'] = 0

    fig = px.line(df, x=time_frame, y='unmet_demand', color='category',
                  color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:,.2%}' + '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, yaxis={'tickformat': '%'})

    return fig


def water_supply_plot(water_delivered, time_frame, layout, by='type'):
    dff_delivered = water_delivered.copy()
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    dff_delivered.loc[dff_delivered['type'].str.contains('WWR'), 'type'] = 'Reused Wastewater'
    dff_delivered = dff_delivered.loc[dff_delivered['type'] != 'Transmission Pipeline']
    dff_delivered = dff_delivered.groupby([time_frame, by])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    # dff_delivered = dff_delivered.loc[dff_delivered['sswd'] != 0]

    df = dff_delivered.reset_index()

    fig = px.area(df, x=time_frame, y='sswd', color=by,
                  color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_traces(fill='tonexty',
                      hovertemplate='<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}')
    fig.update_layout(layout)

    return fig

def energy_demand_plot(water_delivered, wwtp_data, desal_data,
                       time_frame, layout, group_by='type'):
    # emission_factor = 1.76
    dff_energy = water_delivered.copy()
    dff_energy.loc[dff_energy['type'].str.contains('GW'), 'type'] = 'Groundwater pumping'
    dff_energy.loc[dff_energy['type'].str.contains('Pipeline'), 'type'] = 'Surface water conveyance'
    dff_energy.loc[dff_energy['type'].str.contains('DS'), 'type'] = 'Desalinated water conveyance'
    dff_energy = dff_energy.loc[~dff_energy['Supply point'].str.contains('Complexe Aoulouz Mokhtar Soussi')]
    dff_energy = dff_energy[[time_frame, group_by, 'swpa_e']]
    wwtp_data['type'] = 'Wastewater treatment'
    desal_data['type'] = 'Desalination energy'
    for df in [wwtp_data, desal_data]:
        dff = df[[time_frame, group_by, 'swpa_e']]
        dff_energy = dff_energy.append(dff, sort=False)

    df = dff_energy.groupby([time_frame, group_by])[['swpa_e']].sum() / 1000000
    df.reset_index(inplace=True)
    df = df.loc[df['swpa_e'] != 0]

    fig = px.area(df.reset_index(), x=time_frame, y='swpa_e', color=group_by,
                  color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout)

    return fig

def wtd_plot(water_delivered, time_frame, layout):
    dff_wtd = water_delivered.copy()
    dff_wtd = dff_wtd.loc[dff_wtd['type'].str.contains('GW')]
    wtd = dff_wtd.groupby([time_frame, 'Supply point'])['wtd'].mean().reset_index()
    df = wtd

    fig = px.line(df, x=time_frame, y='wtd', color='Supply point',
                  color_discrete_sequence=px.colors.qualitative.Vivid)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}')
    fig.update_layout(layout,
                      yaxis=dict(range=[df.wtd.max()*1.05, df.wtd.min()*0.95]))

    return fig


# TODO: this function would need to be updated if the lcoe results will be shown again
def lcoe_plot(data, ag_lcoe, cropland):
    dff_lcoe = ag_lcoe.dropna().copy()
    dff_lcoe['area_m2'] = dff_lcoe['Demand point'].map(cropland.area_m2)

    dff_lcoe.loc[dff_lcoe['least_cost_technology'].str.contains('Grid'), 'least_cost_technology'] = 'Grid electricity'
    count = len(dff_lcoe['Demand point'].unique())
    # total = dff_lcoe.loc[dff_lcoe['Year']==2020, 'area_m2'].sum()
    name = 'least-cost'
    data[name] = [{'x': dff_lcoe.loc[dff_lcoe['least_cost_technology'] == tech].groupby('Year')[
        'least_cost_technology'].count().index,
                   'y': dff_lcoe.loc[dff_lcoe['least_cost_technology'] == tech].groupby('Year')[
                            'least_cost_technology'].count() / count,
                   'name': tech,
                   'type': 'bar',
                   'hovertemplate': '<b>Value</b>: %{y:,.2%}' + '<br><b>Year</b>: %{x}'
                   } for tech in sorted(dff_lcoe['least_cost_technology'].unique())]
    return data

def plot_borders(geojson):
    line = dict(width=1, color='rgb(80,100,80)')
    colorscale = ((0, 'rgba(80,100,80,0.1)'), (1, 'rgba(80,100,80,0.1)'))
    df = pd.DataFrame(geojson['features'])
    df['color'] = 0

    trace = go.Choroplethmapbox(geojson=geojson, locations=df.id, z=df['color'],
                                colorscale=colorscale,
                                marker=dict(line=line),
                                hoverinfo='skip', showscale=False,
                                showlegend=False)
    return trace


def plot_points(df):
    trace = px.scatter_mapbox(df, lat="lat", lon="lon", color="type",
                              color_discrete_sequence=px.colors.qualitative.Dark2,
                              custom_data=['type', 'point'], hover_name='point').data
    return trace


def plot_pipelines(df):
    color = 'rgb(100,100,100)'
    trace = px.line_mapbox(df, lat="lat", lon="lon", color='type',
                           color_discrete_map={'pipeline': color},
                           custom_data=['type', 'name'], hover_name='name'
                           )
    trace.update_traces(mode='markers+lines', marker=dict(opacity=0))
    return trace.data


def choroplet_map(geojson, df):
    fig = px.choropleth_mapbox(df, geojson=geojson, locations='id',
                               color='color',
                               color_continuous_scale=px.colors.sequential.Viridis,
                               custom_data=['id'])

    fig.update_layout(coloraxis_colorbar=dict(
        len=0.5,
        xanchor="right", x=1,
        yanchor='bottom', y=0.1,
        thickness=10,
    ))
    return fig


def wastewater_treated(df, layout, colors):
    fig = px.area(df, x='Year', y='water',
                  color_discrete_sequence=colors,
                  labels={'water': 'Water (Mm<sup>3</sup>)',
                          'energy': 'Energy (GWh)'})
    fig.update_layout(layout)
    return fig


def water_delivered(df, layout, group_by='type'):
    df.loc[df['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    df.loc[df['type'].str.contains('GW'), 'type'] = 'Groundwater'
    df.loc[df['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    df.loc[df['type'].str.contains('WWR'), 'type'] = 'Reused Wastewater'
    fig = px.area(df, x='Year', y='water',
                  color=group_by,
                  labels={'water': 'Water (Mm<sup>3</sup>)',
                          'energy': 'Energy (GWh)'})
    fig.update_layout(layout)
    return fig


def wastewater_supply(df, layout):
    fig = px.area(df, x='Year', y='water', color='Demand point',
                  labels={'water': 'Water (Mm<sup>3</sup>)',
                          'energy': 'Energy (GWh)'}
                  )
    fig.update_layout(layout)
    return fig


def wwt_energy(df, layout, colors):
    fig = px.area(df, x='Year', y='energy',
                  color_discrete_sequence=colors,
                  labels={'water': 'Water (Mm<sup>3</sup>)',
                          'energy': 'Energy (GWh)'}
                  )
    fig.update_layout(layout)
    return fig


def desal_energy(df, layout):
    # emission_factor = 1.76
    names = {'swpa_e': 'Desalination energy'}
    dff = df.rename(columns=names).melt(id_vars=['Year', 'Supply point'],
                                        value_vars=names.values())
    dff = dff.groupby(['Year', 'variable'])['value'].sum() / 1000000
    dff = dff.reset_index()

    fig = px.area(dff, x='Year', y='value', color='variable',
                  color_discrete_sequence=["#F9ADA0"])
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout)

    return fig


def pumping_energy(df, layout):
    # emission_factor = 1.76
    dff = df.groupby(['Year', 'Supply point'])['swpa_e'].sum() / 1000000
    dff = dff.reset_index()

    fig = px.area(dff, x='Year', y='swpa_e',
                  color_discrete_sequence=["#F9ADA0"])
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout)

    return fig


def crop_production(df, group_by, layout):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']
    dff = df.groupby(['Year', group_by])[['production_kg']].sum() / 1000
    dff = dff.reset_index()
    fig = px.area(dff, x='Year', y='production_kg', color=group_by,
                  color_discrete_sequence=colors,
                  labels={'production_kg': 'Production (ton)'})

    fig.update_layout(layout)
    return fig


def crop_production_per_crop(df, group_by, layout):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']
    dff = df.groupby(['Province', group_by])[['production_kg']].mean() / 1000
    dff.reset_index(inplace=True)
    dff = dff.loc[dff['production_kg'] > 0.0001]
    fig = px.bar(dff, x='production_kg', y=group_by, color='Province',
                 orientation="h",
                 color_discrete_sequence=colors)
    fig.update_layout(layout, height=500,
                      yaxis={'categoryorder': 'total ascending'},
                      legend=dict(
                          orientation="h",
                          # yanchor="top",
                          # y=0,
                          # xanchor="left",
                          x=-0.2
                      ))
    return fig