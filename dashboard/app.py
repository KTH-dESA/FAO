import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import json
import pandas as pd
import numpy as np
import geopandas as gpd
import os.path
import itertools

my_path = os.path.abspath(os.path.dirname(__file__))
# my_path = ''
spatial_data = os.path.join(my_path, 'spatial_data')

# governorates = gpd.read_file(os.path.join(spatial_data,'Admin','JOR_adm1.shp'))
demand_points = gpd.read_file(os.path.join(spatial_data, 'Demand_points.geojson'))
supply_points = gpd.read_file(os.path.join(spatial_data, 'Supply_points.geojson'))
pipelines = gpd.read_file(os.path.join(spatial_data, 'Pipelines.geojson'))
WebMercator = 4326
# governorates.to_crs(epsg=WebMercator, inplace=True)
for gdf in [demand_points, supply_points, pipelines]:
    gdf.to_crs(epsg=WebMercator, inplace=True)

def load_data(scenario, eto, level):
    data_folder = os.path.join(my_path, 'data')
    if not eto:
        eto = ['Without Eto trend']
    data = os.path.join(data_folder, scenario, eto[0], level)

    water_delivered = pd.read_csv(os.path.join(data, 'Water_delivered.csv'))
    water_required = pd.read_csv(os.path.join(data, 'Water_requirements.csv'))
    gw_pumped = pd.read_csv(os.path.join(data, 'Groundwater_pumping.csv'))
    pl_flow = pd.read_csv(os.path.join(data, 'Pipelines_data.csv'))
    wwtp_data = pd.read_csv(os.path.join(data, 'wwtp_data.csv'))
    desal_data = pd.read_csv(os.path.join(data, 'desal_data.csv'))

    return water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data


# water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data = load_data('New Resources', ['Eto trend'],
#                                                                                       'level_1')

token = "pk.eyJ1IjoiY2FtaWxvcmciLCJhIjoiY2p1bTl0MGpkMjgyYjQ0b2E0anRibWJ1MSJ9.GhUUGD6gok1d36lvP17CQQ"

# with open('governorates.json') as response:
#     counties = json.load(response)
#     for feature in counties['features']:
#         feature['id'] = feature['properties']['id']
#
# df = governorates

layout = dict(
    autosize=True,
    automargin=True,
    margin=dict(l=30, r=30, b=100, t=100),
    hovermode="closest",
    plot_bgcolor="#f8f9fa",
    paper_bgcolor="#f8f9fa",
    legend=dict(font=dict(size=10), orientation="h"),
    xaxis={'tickformat': 'd'},
    showlegend=True,
    hoverlabel= {'align': 'left'},
)

app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.BOOTSTRAP],
                # these meta_tags ensure content is scaled correctly on different devices
                # see: https://www.w3schools.com/css/css_rwd_viewport.asp for more
                meta_tags=[
                    {"name": "viewport", "content": "width=device-width, initial-scale=1"}
                ],
                )

server = app.server

# we use the Row and Col components to construct the sidebar header
# it consists of a title, and a toggle, the latter is hidden on large screens
sidebar_header = dbc.Row(
    [
        dbc.Col([html.H3("Jordan"), html.H6('NEXUS model', style={'color': 'gray'})]),
        dbc.Col(
            [
                html.Button(
                    # use the Bootstrap navbar-toggler classes to style
                    html.Span(className="navbar-toggler-icon"),
                    className="navbar-toggler",
                    # the navbar-toggler classes don't set color
                    style={
                        "color": "rgba(0,0,0,.5)",
                        "border-color": "rgba(0,0,0,.1)",
                    },
                    id="navbar-toggle",
                ),
                html.Button(
                    # use the Bootstrap navbar-toggler classes to style
                    html.Span(className="navbar-toggler-icon"),
                    className="navbar-toggler",
                    # the navbar-toggler classes don't set color
                    style={
                        "color": "rgba(0,0,0,.5)",
                        "border-color": "rgba(0,0,0,.1)",
                    },
                    id="sidebar-toggle",
                ),
            ],
            # the column containing the toggle will be only as wide as the
            # toggle, resulting in the toggle being right aligned
            width="auto",
            # vertically align the toggle in the center
            align="center",
        ),
    ],
    id='sidebar-header'
)

scenario_options = html.Div(
    [
        html.H6('Select scenario'),
        html.Div(
            dbc.RadioItems(
                id="rb-scenario",
                options=[
                    {"label": "Reference", "value": 'Reference'},
                    {"label": "Improved agricultural efficiency", "value": 'Improve AG eff'},
                    {"label": "New water resources", "value": 'New Resources'},
                    {"label": "Reduce non-revenue water", "value": 'Reduce NRW'},
                ],
                value='Reference',
            ),
        ),
        html.Hr(),
    ],
)

eto_options = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(
                    dbc.Checklist(
                        options=[
                            {"label": "", "value": 'Eto trend'},
                        ],
                        value=['Eto trend'],
                        id="eto-input",
                        switch=True,
                    ),
                    width=2,
                ),
                dbc.Col(html.H6('Increase Evapotranspiration (climate change)'), width=8),
            ],
            no_gutters=True
        ),
        html.Hr(),
    ]
)

level_options = html.Div(
    [
        html.H6('Select level of variable'),
        html.Div(
            dcc.Dropdown(
                id="drop-level",
                options=[
                    {"label": "Level 1", "value": 'level_1'},
                    # {"label": "Level 2", "value": 2},
                ],
                value='level_1',
                clearable=False
            ),
        ),
        html.Hr(),
    ]
)

energy_options = html.Div(
    [
        html.H6('Select energy pumping efficiency'),
        dbc.Row([
            # dcc.Slider(
            #     id="pump-eff-level",
            #     min=0.4,
            #     max=0.9,
            #     value=0.6,
            #     step=0.05,
            #     marks={i/100: f'{i}%' for i in range(40, 100, 10)}
            #     #marks={y: y for y in range(int(years.min()), 2031, 2)},
            # ),
            dbc.Col(
                dbc.InputGroup(
                    [
                        dbc.InputGroupAddon("Current", addon_type="prepend"),
                        dbc.Input(placeholder="value", type="number", value=0.45,
                                  step=0.05, min=0.3, max=0.85, id='pump-eff-init'),
                    ],
                    size="sm"
                ),
                width='6'
            ),  # html.Br(),
            dbc.Col(
                dbc.InputGroup(
                    [
                        dbc.InputGroupAddon("Goal", addon_type="prepend"),
                        dbc.Input(placeholder="value", type="number", value=0.45,
                                  step=0.05, min=0.3, max=0.85, id='pump-eff-end'),
                    ],
                    size="sm"
                ),
                width={"size": 5, "offset": 1}

            )
        ],
            no_gutters=True,
        ),
    ]
)

scenario_tools = html.Div(
    [
        scenario_options,
        eto_options,
        level_options,
        energy_options
    ],
    id='tools',
    className='container',
)

footer = dbc.Row(
            [
                dbc.Button("Apply", color="info", className="mr-1",
                           style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-apply'),
            ],
            align='center',
            justify="end",
            className='footer',
        )

sidebar = html.Div(
    [
        sidebar_header,
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("Scenario", active=True, href="#")),
                dbc.NavItem(dbc.NavLink("Visualisation", href="#")),
            ],
            id="tabs",
        ),
        # we wrap the horizontal rule and short blurb in a div that can be
        # hidden on a small screen
        scenario_tools,
        footer
        # id="blurb",
        # use the Collapse component to animate hiding / revealing links
    ],
    id="sidebar",
)

results_header = dbc.Row(
    [
        dbc.Col([html.H5("Summary of results:"), html.Div([html.H6(id='resultsTitle', style={'color': 'gray',
                                                                                             'display': 'inline'}),
                                                           html.H6(id='scenarioTitle', style={'color': 'gray',
                                                                                              'display': 'inline'})],
                                                          className='flex')]),
    ],
    id='results-header',
)

footer_results = dbc.Row(
            [
                dbc.Button("Download", color="info", className="mr-1",
                           style={'font-size': '0.85rem', 'font-weight': '600'}, id='button-download'),
            ],
            align='center',
            justify="center",
            className='footer',
        )

map = dcc.Graph(id="map",
                className='col-lg-7 col-md-12 col-sm-12 col-xs-12')

graphs = html.Div([results_header,
                   html.Div(dbc.Col(id='graphs'), id='graphs-container'),
                   footer_results
                   ],
                  id='results-container',
                  className='col-lg-5 col-md-12 col-sm-12 col-xs-12')

content = html.Div([map, graphs], id="page-content", className='row')

app.layout = html.Div([dcc.Store(id='water'), sidebar, content])


# Helper funtions
# def choroplethmap():
#     data = [dict(
#         type="choroplethmapbox",
#         geojson=counties,
#         locations=df.NAME_1,
#         name=df.NAME_1,
#         z=df.GLCV2sum,
#         customdata=[{'Water demand': np.random.randint(1, 101, 12),
#                      'Energy demand': np.random.randint(1, 101, 12),
#                      'Agricultural yield': np.random.randint(1, 101, 12)} for i in range(0, df.shape[0])],
#         colorbar={'len': 0.3, 'thicknessmode': 'fraction',
#                   'thickness': 0.01, 'x': 0.01, 'y': 0.85,
#                   'title': {'text': 'Cropland density'}},
#         showscale=True
#     )]
#     return data


def scatterpointmap(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data):
    df_pipelines = pipelines.groupby('index').agg({'pipeline': 'first',
                                                   'segment_length_m': 'first',
                                                   'geometry': 'first'}).reset_index()
    dff_pipeflow = pl_flow.groupby(['Year', 'pipeline']).agg({'water_use': lambda value: sum(value) / 1000000,
                                                              'SWPA_E_': lambda value: sum(value) / 1000000,
                                                              'pipeline_length': 'first'}).reset_index()
    data = [dict(
        type="scattermapbox",
        mode='lines+markers',
        marker=dict(opacity=0),
        # line=dict(color='gray'),
        lon=list(itertools.chain.from_iterable([list(line.coords.xy[0]) + [None] for line in df_pipelines.geometry])),
        lat=list(itertools.chain.from_iterable([list(line.coords.xy[1]) + [None] for line in df_pipelines.geometry])),
        text=list(itertools.chain.from_iterable([len(list(line.coords.xy[0])) * [pipe] + [None] for line, pipe in
                                                 zip(df_pipelines.geometry, df_pipelines.pipeline)])),
        hoverinfo='name+text',
        name='Pipeline',
        showlegend=True,
        customdata=list(itertools.chain.from_iterable([[{
            'water': [dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].water_use,
                      dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].Year],
            'energy': [dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].SWPA_E_,
                       dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].Year],
            'type': "pipeline"} for i in range(len(list(line.coords.xy[0])))] + [{}] for line, pipe in
                                                       zip(df_pipelines.geometry, df_pipelines.pipeline)])),
    )]

    df_demand = demand_points.groupby('point').agg({'type': 'first',
                                                    'elevation_m': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_delivered = water_delivered.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_required = water_required.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_required = dff_required.reset_index()
    # dff_unmet = dff_required.copy()
    # dff_unmet['value'] = dff_unmet.value - dff_unmet.set_index(['Year','point']).index.map(dff_delivered.set_index(['Year','point']).value)
    # dff_unmet.loc[dff_unmet['value']<0.001, 'value'] = 0
    names = ['Water delivered (Mm3)', 'Water required (Mm3)']
    # dff['Date'] = [pd.Timestamp(year=i, month=j, day=1) for i, j in zip(dff.Year, dff.Month)]
    data += [dict(
        type="scattermapbox",
        lon=[point.x for point in df_demand.loc[df_demand['type'] == type].geometry],
        lat=[point.y for point in df_demand.loc[df_demand['type'] == type].geometry],
        text=df_demand.loc[df_demand['type'] == type].point,
        ids=df_demand.loc[df_demand['type'] == type].point,
        hoverinfo='name+text',
        name=type,
        showlegend=True,
        customdata=[{names[0]: [dff_delivered.loc[(dff_delivered.point == point)].value,
                                dff_delivered.loc[(dff_delivered.point == point)].Year],
                     names[1]: [dff_required.loc[(dff_required.point == point)].value,
                                dff_required.loc[(dff_required.point == point)].Year],
                     'type': 'demand',
                     # names[2]: [dff_unmet.loc[(dff_unmet.point == point)].value,
                     #            dff_unmet.loc[(dff_unmet.point == point)].Year],
                     } for point in df_demand.loc[df_demand.type == type].point],
    ) for type in sorted(df_demand['type'].unique())]

    df_supply = supply_points.groupby('point').agg({'type': 'first',
                                                    'wtd_m': 'first',
                                                    'elevation_m': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_gw = gw_pumped.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    dff_wwtp = wwtp_data.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    dff_desal = desal_data.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    dff_surface = pl_flow.groupby(['Year', 'type', 'point'])[['water_use']].sum() / 1000000
    dff_surface.rename(columns={'water_use': 'value'}, inplace=True)
    dff_surface = dff_surface.reset_index()
    dff_gw = dff_gw.reset_index()
    dff_wwtp = dff_wwtp.reset_index()
    dff_desal = dff_desal.reset_index()
    dff = dff_gw.append([dff_wwtp, dff_surface, dff_desal], sort=False, ignore_index=True)
    data += [dict(
        type="scattermapbox",
        lon=[point.x for point in df_supply.loc[df_supply['type'] == type].geometry],
        lat=[point.y for point in df_supply.loc[df_supply['type'] == type].geometry],
        text=df_supply.loc[df_supply['type'] == type].point,
        ids=df_supply.loc[df_supply['type'] == type].point,
        hoverinfo='name+text',
        name=type,
        showlegend=True,
        customdata=[{'water': [dff.loc[(dff.point == point)].value,
                               dff.loc[(dff.point == point)].Year],
                     'energy': [dff.loc[(dff.point == point)].SWPA_E_,
                                dff.loc[(dff.point == point)].Year],
                     'type': type,
                     } for point in df_supply.loc[df_supply.type == type].point],
    ) for type in sorted(df_supply['type'].unique())]

    return data


def plot_map(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data):
    layout_map = layout.copy()

    data = scatterpointmap(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data)

    layout_map["mapbox"] = {"center": {"lon": 36.5, 'lat': 31.2}, 'zoom': 6,
                            'style': "light", 'accesstoken': token}
    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'event+select'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)
    map = dict(data=data, layout=layout_map)
    return map


# @app.callback(
#     Output("water", "data"),
#     [Input("button-apply", "n_clicks"), Input("sidebar-toggle", "n_clicks"), Input('water', 'modified_timestamp')],
#     [State('rb-scenario', 'value'), State('eto-input', 'value'), State('drop-level', 'value')]
# )
# def update_production_text(n_1, n_2, ts, scenario, eto, level):
#     if ts is None:
#         raise PreventUpdate
#     water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data = load_data(scenario,
#                                                                                            eto,
#                                                                                            level)
#     return [water_delivered.to_dict(), water_required.to_dict()]

@app.callback(
    Output("sidebar", "className"),
    [Input("sidebar-toggle", "n_clicks")],
    [State("sidebar", "className")],
)
def toggle_classname(n, classname):
    if n and classname == "":
        return "collapsed"
    return ""


@app.callback(
    [Output("graphs", "children"), Output('resultsTitle', 'children')],
    [Input('map', 'selectedData'),
     Input("sidebar-toggle", "n_clicks"),
     Input('scenarioTitle', 'children'),
     Input('pump-eff-init', 'value'),
     Input('pump-eff-end', 'value')],
    [State('rb-scenario', 'value'), State('eto-input', 'value'), State('drop-level', 'value')]
)
def update_results(selection, n_1, title, eff_init, eff_end, scenario, eto, level):
    # time.sleep(1)

    names = ['Water delivered (Mm3)', 'Water required (Mm3)']
    colors = {'water': "#59C3C3", 'energy': "#F9ADA0", 'food': "#849E68"}
    data = {}
    if type(selection) == type(None):
        emission_factor = 0.643924449
        water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data = load_data(scenario,
                                                                                               eto,
                                                                                               level)
        # water_delivered = pd.DataFrame(water_data[0])
        # water_required = pd.DataFrame(water_data[0])
        name = 'Jordan country'
        dff_delivered = water_delivered.groupby(['Year', 'type'])['value'].sum() / 1000000
        dff_delivered = dff_delivered.reset_index()
        dff_required = water_required.groupby(['Year', 'type'])['value'].sum() / 1000000
        dff_required = dff_required.reset_index()

        dff_unmet = dff_required.copy()
        dff_unmet['value'] = (dff_unmet.value - \
                             dff_unmet.set_index(['Year', 'type']).index.map(
                                 dff_delivered.set_index(['Year', 'type']).value)) / dff_unmet.value * 100
        dff_unmet.loc[dff_unmet['value'] < 0.001, 'value'] = 0

        dff_energy = pd.DataFrame()
        pl_flow['type'] = 'Water conveyance'
        l = len(pl_flow.Year.unique())
        eff = np.array([(eff_end - eff_init) / (l - 1) * i + eff_init for i in range(l)])
        for df in [gw_pumped, pl_flow]:
            dff = df.groupby(['Year', 'type'])['SWPA_E_'].sum() / (eff * 1000000)
            dff = dff.reset_index()
            dff_energy = dff_energy.append(dff, sort=False)
        for df in [wwtp_data, desal_data]:
            dff = df.groupby(['Year', 'type'])['SWPA_E_'].sum() / (1000000)
            dff = dff.reset_index()
            dff_energy = dff_energy.append(dff, sort=False)
        dff_energy.rename(columns={'SWPA_E_': 'value'}, inplace=True)

        for df, name in zip([dff_delivered, dff_required], [names[0], names[1]]):
            data[name] = [{'x': df.loc[df['type'] == type].Year,
                           'y': df.loc[df['type'] == type].value,
                           'name': type,
                           'stackgroup': 'one',
                           'mode': 'lines',
                           'text': dff_unmet.loc[dff_unmet['type'] == type].value,
                           'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}' +
                                            '<br><b>Unmet demand</b>: %{text:.2f}%'
                           } for type in sorted(df['type'].unique())]

        for df, name in zip([dff_energy], ['Energy demand (GWh)']):
            data[name] = [{'x': df.loc[df['type'] == type].Year,
                           'y': df.loc[df['type'] == type].value,
                           'name': type,
                           'stackgroup': 'one',
                           'mode': 'lines',
                           'text': df.loc[df['type'] == type].value * emission_factor/1000,
                           'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}' +
                                            '<br><b>Emissions</b>: %{text: 0.2f} MtCO2'
                           } for type in sorted(df['type'].unique())]


    elif selection['points'][0]['customdata']['type'] in ['demand']:
        name = selection['points'][0]['text']
        name_key = selection['points'][0]['customdata']['type']
        name_dict = {'demand': {'water': 'Water demand (Mm3)'}}
        data[name_dict[name_key]['water']] = [{'x': selection['points'][0]['customdata'][name][1],
                                               'y': selection['points'][0]['customdata'][name][0],
                                               'fill': 'tonexty',
                                               'mode': 'lines',
                                               'name': name,
                                               } for name in names]
    elif selection['points'][0]['customdata']['type'] in ['Groundwater supply', 'WWTP', 'pipeline', 'Desalination']:
        name = selection['points'][0]['text']
        name_key = selection['points'][0]['customdata']['type']
        name_dict = {'Groundwater supply': {'water': 'Water supplied (Mm3)', 'energy': 'Energy for pumping (GWh)'},
                     'WWTP': {'water': 'Treated wastewater (Mm3)', 'energy': 'Energy for treatment (GWh)'},
                     'pipeline': {'water': 'Water conveyed (Mm3)', 'energy': 'Energy for pumping (GWh)'},
                     'Desalination': {'water': 'Water desalinated (Mm3)', 'energy': 'Energy for desalination (GWh)'}}

        if name_key in ['Groundwater supply', 'pipeline']:
            l = len(selection['points'][0]['customdata']['water'][1])
            eff = np.array([(eff_end - eff_init) / (l - 1) * i + eff_init for i in range(l)])
        else:
            eff = 1
        data[name_dict[name_key]['water']] = [{'x': selection['points'][0]['customdata']['water'][1],
                                               'y': selection['points'][0]['customdata']['water'][0],
                                               'fill': 'tozeroy', 'mode': 'lines', 'showlegend': False,
                                               'line': {'color': colors['water']}}]
        data[name_dict[name_key]['energy']] = [{'x': selection['points'][0]['customdata']['energy'][1],
                                                'y': np.array(selection['points'][0]['customdata']['energy'][0])/eff,
                                                'fill': 'tozeroy', 'mode': 'lines', 'showlegend': False,
                                                'line': {'color': colors['energy']}}]
    elif selection['points'][0]['customdata']['type'] in ['River/pipeline supply']:
        name = selection['points'][0]['text']
        name_key = selection['points'][0]['customdata']['type']
        name_dict = {'River/pipeline supply': {'water': 'Water supplied (Mm3)'}}
        data[name_dict[name_key]['water']] = [{'x': selection['points'][0]['customdata']['water'][1],
                                               'y': selection['points'][0]['customdata']['water'][0],
                                               'fill': 'tozeroy', 'mode': 'lines', 'showlegend': False,
                                               'line': {'color': colors['water']}}]

    plots = []
    for key, value in data.items():
        layout_plot = layout.copy()
        layout_plot['title'] = dict(text=key)
        layout_plot['height'] = 400
        layout_plot['font'] = dict(size=11, color="#7f7f7f")
        # layout_plot['showlegend'] = False
        plots.append(dcc.Graph(figure=dict(data=value, layout=layout_plot)))
        # data[value][0]['line'] = dict(shape='lines+markers', color=colors[i])
    return plots, name


@app.callback(
    [Output('drop-level', 'options'), Output('drop-level', 'value'), Output('drop-level', 'disabled')],
    [Input('rb-scenario', 'value')]
)
def update_level_dropdown(scenario):
    level_dict = {'Reference': {'Select...': 'level_1'},
                  'Improve AG eff': {'by 10 percent': 'level_1',
                                     'by 20 percent': 'level_2'},
                  'New Resources': {'Select...': 'level_1'},
                  'Reduce NRW': {'to 40 percent': 'level_1',
                                 'to 20 percent': 'level_2'}}
    options = [{"label": key, 'value': value} for key, value in level_dict[scenario].items()]
    disable = False
    if len(level_dict[scenario].keys()) == 1:
        disable = True
    return options, 'level_1', disable


@app.callback(
    [Output("map", "figure"), Output('scenarioTitle', 'children')],
    [Input("button-apply", "n_clicks"), Input("sidebar-toggle", "n_clicks")],
    [State('rb-scenario', 'value'), State('eto-input', 'value'), State('drop-level', 'value')]
)
def update_level_dropdown(n_1, n_2, scenario, eto, level):
    level_dict = {'Reference': {'level_1': ''},
                  'Improve AG eff': {'level_1': 'by 10 percent',
                                     'level_2': 'by 20 percent'},
                  'New Resources': {'level_1': ''},
                  'Reduce NRW': {'level_1': 'to 40 percent',
                                 'level_2': 'to 20 percent'}}

    name = f' - {scenario} {level_dict[scenario][level]} scenario'

    # global water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data, map

    water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data = load_data(scenario,
                                                                                           eto,
                                                                                           level)
    map = plot_map(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data)

    return map, name


if __name__ == "__main__":
    app.run_server(debug=True)
