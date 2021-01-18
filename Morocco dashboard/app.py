# -*- coding: utf-8 -*-

import itertools
import os.path

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import geopandas as gpd
import numpy as np
import pandas as pd
import yaml
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame
from dash_extensions.snippets import send_file
import plotly.graph_objects as go
import plotly.io as pio
import json
import boto3

from scripts import plotting

pio.templates.default = "plotly_white"

client = boto3.client('s3')
# my_path = os.path.abspath(os.path.dirname(__file__))
my_path = os.path.join('..', 'Morocco model', 'test dash results')
spatial_data = os.path.join(my_path, 'spatial_data')

# provinces = gpd.read_file(os.path.join(spatial_data, 'Admin', 'provinces.geojson'))
cropland = pd.read_csv(os.path.join(spatial_data, 'cropland.csv'))
demand_points = gpd.read_file(os.path.join(spatial_data, 'Demand_points.gpkg'))
demand_points.loc[demand_points['type'] == 'Catchment', 'type'] = 'Agriculture'
demand_points.loc[demand_points['type'] == 'Demand site', 'type'] = 'Municipality' #TODO: move this to the schematic processing
supply_points = gpd.read_file(os.path.join(spatial_data, 'Supply_points.gpkg'))
supply_points.loc[supply_points['type'] == 'Other supply', 'type'] = 'Desalination plant'
pipelines = gpd.read_file(os.path.join(spatial_data, 'Pipelines.gpkg'))
WebMercator = 4326

# centroids = provinces.centroid

with open(os.path.join(spatial_data, 'Admin', 'provinces.geojson')) as response:
    provinces = json.load(response)
    for value in provinces['features']:
        value['id'] = value['properties']['id']

for gdf in [demand_points, supply_points, pipelines]:
    gdf.to_crs(epsg=WebMercator, inplace=True)

points_coords, pipe_coords = plotting.data_merging(demand_points, supply_points, pipelines)


def load_data(path, scenario, climate, phaseout_year, pv_level, files='all'):
    init_year = 2020
    butane_scenario = f'phaseout_{phaseout_year}' if phaseout_year != 2050 else 'phaseout_None'
    data_folder = os.path.join(path, 'data')
    if not climate:
        climate = ['Trend']
    data = os.path.join(data_folder, scenario, climate[0])
    # lcoe = os.path.join(data_folder, scenario, climate[0], level)

    if files == 'all':
        files = ['results.gz', 'wwtp_data.gz', 'desal_data.gz', 'summary_results.gz']

    if isinstance(files, str):
        files = [files]

    if len(files) == 1:
        dff = pd.read_csv(os.path.join(data, files[0]))
        output = dff
    else:
        output = []
        for file in files:
            if file == 'summary_results.gz':
                dff = pd.read_csv(os.path.join(path,
                                               'Butane_calculations',
                                               butane_scenario,
                                               f'{pv_level}_PV',
                                               file))
            else:
                dff = pd.read_csv(os.path.join(data, file))
            dff = dff.loc[dff.Year >= init_year]
            output.append(dff)
    return output


# def load_data(scenario, climate, level, phaseout_year, pv_level):
#     init_year = 2020
#     butane_scenario = f'phaseout_{phaseout_year}' if phaseout_year != 2050 else 'phaseout_None'
#     data_folder = 's3://souss-massa-project/data'
#     butane_folder = 's3://souss-massa-project/Butane_calculations'
#     if not climate:
#         climate = ['Trend']
#     data = f'{data_folder}/{scenario}/{climate[0]}'
#     # lcoe = f'{data_folder}/{scenario}/{climate[0]}/{level}'
#
#     water_delivered = pd.read_csv(f'{data}/results.gz')
#     # ag_lcoe = pd.read_csv(f'{lcoe}/lcoe.gz')
#     ag_lcoe = pd.DataFrame({'Year': []})
#     wwtp_data = pd.read_csv(f'{data}/wwtp_data.gz')
#     desal_data = pd.read_csv(f'{data}/desal_data.gz')
#     desal_data['swpa_e'] = desal_data['sswd'] * 3.31
#     butane_data = pd.read_csv(f'{butane_folder}/{butane_scenario}/{pv_level}_PV/summary_results.gz')
#
#     return water_delivered.loc[water_delivered.Year>=init_year], ag_lcoe.loc[ag_lcoe.Year>=init_year], \
#            wwtp_data.loc[wwtp_data.Year>=init_year], desal_data.loc[desal_data.Year>=init_year], \
#            butane_data.loc[butane_data.Year >= init_year]


button_color = 'primary'
info_ids = []
title_ids = []


def title_info(title_id, info_id, modal_size):
    info_ids.append(info_id)
    title_ids.append(title_id)
    return dbc.Row([dbc.Col(html.H6(id=title_id), width=10),
                    dbc.Col(html.Button(html.Span(className="fa fa-info-circle"), className='info', id=info_id),
                            width=2),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(id=f'{info_id}-title'),
                            dbc.ModalBody(id=f'{info_id}-body'),
                            dbc.ModalFooter(
                                dbc.Button("Close", id=f"{info_id}-close", className="ml-auto")
                            ),
                        ],
                        size=modal_size,
                        id=f"{info_id}-modal",
                    )])


token = open('.mapbox_token').read()

layout = dict(
    autosize=True,
    height=350,
    margin=dict(l=30, r=20, b=50, t=100),
    hovermode="closest",
    plot_bgcolor="#fff",
    paper_bgcolor="#fff",
    legend=dict(font=dict(size=10), orientation="h"),
    showlegend=True,
    hoverlabel={'align': 'left'},
    xaxis={'title': ''},
    font=dict(size=10, color="#7f7f7f"),
    yaxis={'title': ''},
)

app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.BOOTSTRAP,
                                      'https://use.fontawesome.com/releases/v5.12.0/css/all.css'],
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
        dbc.Col([html.H3("Souss-Massa"), html.H6('NEXUS model', id='title', style={'color': 'gray'})]),
        dbc.Col(
            [
                dbc.Button(
                    html.Span(className="fa fa-bars"),
                    color="secondary",
                    outline=True,
                    id="navbar-toggle",
                ),
                dbc.Button(
                    html.Span(className="fa fa-bars"),
                    color="secondary",
                    outline=True,
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
        title_info(title_id='scenarios-title', info_id='scenario-info',
                   modal_size='lg'),
        html.Div(
            dbc.RadioItems(
                id="rb-scenario",
                value='Reference',
                className='checklist-selected-style',
            ),
        ),
    ],
    className='options'
)

climate_options = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(
                    dbc.Checklist(
                        options=[
                            {"label": "", "value": 'Climate Change'},
                        ],
                        value=[],
                        id="climate-input",
                        switch=True,
                    ),
                    width=2,
                ),
                dbc.Col(title_info(title_id='climate', info_id='climate-info', modal_size='lg')),
            ],
            no_gutters=True,
        )
    ],
    className='options'
)

technology_options = html.Div(
    [
        title_info(title_id='cost-reduction', info_id='technology-info', modal_size='lg'),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-wind"), html.Label('Wind power', id='wind-title',
                                                                    style={'padding': '0 0 0 10px'})],
                        style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=0.7,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'},
                                          0.7: {'label': '70%'}},
                                   included=False, id='rate-wind'),
                        style={'margin-top': '5px', 'margin-left': '0.5em'}
                        ),
            ],
            no_gutters=True,
            style={'margin-top': '1em', 'margin-bottom': '1em'}
        ),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-solar-panel"), html.Label('Solar panels', id='pv-title',
                                                                           style={'padding': '0 0 0 8px'})],
                        style={'font-size': '14px', 'color': 'gray', 'margin-right': '0em'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=None,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'},
                                          0.7: {'label': '70%'}},
                                   included=False, id='rate-pv'),
                        style={'margin-top': '5px', 'margin-left': '0.5em'}
                        ),
            ],
            no_gutters=True,
        ),
    ],
    className='options',
    hidden=True,
)

grid_options = html.Div(
    [
        title_info(title_id='price-increase', info_id='grid-info', modal_size='lg'),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-plug"), html.Label('Grid electricity', id='grid-title',
                                                                    style={'padding': '0 0 0 10px'})],
                        style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=None,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'},
                                          0.7: {'label': '70%'}},
                                   included=False, id='rate-grid'),
                        style={'margin-top': '5px', 'margin-left': '0.5em'}
                        ),
            ],
            no_gutters=True, style={'margin-top': '1em'}
        ),
    ],
    className='options',
    hidden=True,
)

butane_options = html.Div(
    [
        title_info(title_id='butane-title', info_id='butane-info', modal_size='lg'),
        dbc.Row(
            [
                dbc.Col(
                    [html.I(className="fa fa-burn"), html.Label(id='butane-phaseout', style={'padding': '0 0 0 10px'})],
                    style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=2030, max=2050, value=2050, step=None,
                                   marks={2030: {'label': '2030'}, 2040: {'label': '2040'}, 2050: {'label': 'none'}},
                                   included=False, id='butane-year'),
                        style={'margin-top': '5px', 'margin-left': '0.5em'}
                        ),
            ],
            no_gutters=True, style={'margin-top': '1em'}
        ),
    ],
    className='options',
    # hidden=True,
)

pv_options = html.Div(
    [
        title_info(title_id='pv-adoption-title', info_id='pv-info', modal_size='lg'),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-solar-panel"),
                         html.Label(id='pv-adoption', style={'padding': '0 0 0 10px'})],
                        style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=10, max=50, value=10, step=None,
                                   marks={10: {'label': '10%'}, 20: {'label': '20%'}, 50: {'label': '50%'}},
                                   included=False, id='pv-share'),
                        style={'margin-top': '5px', 'margin-left': '0.5em'}
                        ),
            ],
            no_gutters=True, style={'margin-top': '1em'}
        ),
    ],
    className='options',
    # hidden=True,
)


# map_options = html.Div(
#     [
#         title_info(title='Map visualization options', info_id='map-info',
#                    info_title='Map visualization info', modal_content='This is the content of the modal'),
#         html.Div(
#             dbc.RadioItems(
#                 id="map-options",
#                 options=[
#                     {"label": "System schematic", "value": 'sch-map'},
#                     {"label": "By province", "value": 'cho-map'},
#                 ],
#                 value='sch-map',
#                 className='checklist-selected-style',
#             ),
#         ),
#         dbc.Collapse(dbc.Card(dbc.CardBody(
#             [
#                 dcc.Dropdown(id='cho-map-drop',
#                              options=[
#                                  {"label": "Cropland density", "value": 'Cropland density'},
#                                  {"label": "Water delivered", "value": 'Water delivered'},
#                                  {"label": "Energy demand", "value": 'Energy demand'},
#                              ],
#                              value='Cropland density',
#                              placeholder='Select variable...',
#                              clearable=False,
#                              style={'marginBottom': '0em'}
#                              ),
#                 # dcc.Dropdown(id='cho-map-filter',
#                 #              options=[
#                 #                  {"label": "Province", "value": 'prov'},
#                 #                  {"label": "Irrigation district", "value": 'irr'},
#                 #              ],
#                 #              value='crop',
#                 #              placeholder='Filter by...',
#                 #              clearable=False
#                 #              ),
#             ])),
#             id="cho-map-collapse",
#         ),
#     ],
#     className='options'
# )

scenario_tools = html.Div(
    [
        scenario_options,
        html.Hr(),
        climate_options,
        html.Hr(),
        technology_options,
        # html.Hr(),
        grid_options,
        # html.Hr(),
        butane_options,
        html.Hr(),
        pv_options
        # unit_options,
        # html.Hr(),
        # map_options,
        # html.Hr(),
        # compare_scenarios
    ],
    id='tools',
)

visual_tools = html.Div(
    # [
    #     unit_options,
    #     html.Hr(),
    #     map_options,
    #     html.Hr(),
    #     compare_scenarios
    # ],
    id='visual-tools',
)

footer = dbc.Row(
    [
        dbc.Button([html.I(className='fa fa-redo-alt'), " Reset"], color=button_color, outline=True, className="mr-1",
                   style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-reset'),
        dbc.Button([html.I(className='fa fa-check-double'), " Apply"], color=button_color, className="mr-1",
                   style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-apply'),
    ],
    align='center',
    justify="around",
    className='footer',
)

sidebar = html.Div(
    [
        sidebar_header,
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("Scenario options", active=True, href="#", id='page-1', className='tabs')),
                # dbc.NavItem(dbc.NavLink("Visualisation", href="#", id='page-2', className='tabs')),
            ],
            id="tabs",
        ),
        dbc.Collapse([
            scenario_tools,
            visual_tools,
            # dbc.Input(type='range', className='custom-range')
            # use the Collapse component to animate hiding / revealing links
        ],
            id="collapse",
        ),
        footer,
    ],
    id="sidebar",
)

results_header = dbc.Row(
    [
        dbc.Col([html.H5("Summary of results:", id='resultsTitle'),
                 html.Div([html.H6(id='resultsSubTitle', style={'color': 'gray',
                                                                'display': 'inline'}),
                           html.H6(id='scenarioTitle', style={'color': 'gray',
                                                              'display': 'inline'})],
                          className='flex')]),
    ],
    id='results-header',
    style={'margin-top': '-1em'}
)

footer_results = dbc.Row(
    [
        dbc.Button([html.I(className='fa fa-chart-pie'), " Compare"], color=button_color, outline=True,
                   className="mr-1", style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='compare'),
        dbc.Modal(
            [
                dbc.ModalBody([dbc.Tabs([dbc.Tab(tab_id="tab-1", id="tab-1"),
                                         dbc.Tab(tab_id="tab-2", id="tab-2"),
                                         dbc.Tab(tab_id="tab-3", id="tab-3"),
                                         dbc.Tab(tab_id="tab-4", id="tab-4")],
                                        id="compare-tabs",
                                        active_tab="tab-1",
                                        style={'fontSize': '1.5rem', 'fontWeight': '600'}
                                        ),
                               dcc.Store(id='compare-data'),
                               dcc.Loading(dbc.ModalBody([html.Div(html.Div(style={'height': '300px'}),
                                                                   id="compare-content")]))]),
                dbc.ModalFooter(
                    dbc.Button("Close", id="compare-close", className="ml-auto")
                ), ],
            size="xl",
            id="compare-modal",
        ),
        Download(id="download"),
        dbc.Button([html.I(className='fa fa-download'), " Download"], color=button_color,
                   className="mr-1", style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-download'),
    ],
    align='center',
    justify="around",
    className='footer',
)

map = html.Div(dcc.Loading(
    id="loading-1",
    type="default",
    children=dcc.Graph(id="map",
                       # className='col-lg-7 col-md-7 col-sm-12 col-xs-12',
                       config=dict(showSendToCloud=True,
                                   toImageButtonOptions=dict(format='png', filename='map', height=700,
                                                             width=700, scale=2)
                                   )
                       ),
),
    id='map-div',
    className='col-xl-7 col-lg-12',
)
# map = dcc.Graph(id="map",
#                className='col-xl-7 col-lg-12',
#                config=dict(showSendToCloud=True,
#                            toImageButtonOptions=dict(format='png', filename='map', height=700,
#                                                      width=700, scale=2)
#                            )
#                )


graphs = html.Div([dcc.Store(id='current-language'), dbc.Nav([dbc.DropdownMenu(
    [dbc.DropdownMenuItem(["English ", html.I(className='fa fa-language')], className="drop-items", id="english"),
     dbc.DropdownMenuItem(["Spanish ", html.I(className='fa fa-language')], className="drop-items", id="spanish"),
     dbc.DropdownMenuItem(["French ", html.I(className='fa fa-language')], className="drop-items", id="french")],
    label="Language", id='language', nav=True),
    dbc.NavLink("About", id='about', href="#"),
    dbc.Modal(
        [
            dbc.ModalHeader(id='about-title'),
            dbc.ModalBody(id='about-body'),
            dbc.ModalFooter(
                dbc.Button("Close", id="about-close", className="ml-auto")
            ),
        ],
        size="lg",
        id="about-modal",
    )
],
    horizontal='end', style={'margin-bottom': '0em'}, id='about-nav'),
                   results_header,
                   html.Div(dbc.Col(dcc.Loading(id='graphs', type="default")),
                            id='graphs-container'),
                   footer_results
                   ],
                  id='results-container',
                  className='col-xl-5 col-lg-12')

content = html.Div([map, graphs], id="page-content")

app.layout = html.Div([dcc.Store(id='current'), sidebar, content])


# Helper funtions
# with open(os.path.join(spatial_data, 'Admin', 'Provinces.geojson')) as response:
#     provinces_chmap = json.load(response)
#     for feature in provinces_chmap['features']:
#         # print(feature)
#         feature['id'] = feature['properties']['Province']

def get_language(language):
    file = f"assets/{language}.yaml"
    with open(file, 'rt', encoding='utf8') as yml:
        language_dic = yaml.load(yml, Loader=yaml.FullLoader)
    return language_dic


def choroplethmap(geojson, locations, title):
    layout_map = layout.copy()
    data = [dict(
        type="choroplethmapbox",
        geojson=geojson,
        locations=locations,
        name=locations,
        z=np.random.randint(1, 101, locations.shape[0]),
        # customdata=[{'Water demand': np.random.randint(1, 101, 12),
        #              'Energy demand': np.random.randint(1, 101, 12),
        #              'Agricultural yield': np.random.randint(1, 101, 12)} for i in range(0, locations.shape[0])],
        colorbar={'len': 0.3, 'thicknessmode': 'fraction',
                  'thickness': 0.01, 'x': 0.01, 'y': 0.85,
                  'title': {'text': title}},
        showscale=True,
    )]

    layout_map["mapbox"] = {"center": {"lon": -8.9, 'lat': 30.2}, 'zoom': 7,
                            'style': "outdoors", 'accesstoken': token}
    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'event+select'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)
    map = dict(data=data, layout=layout_map)

    return map


def plot_map(background):
    layout_map = {}

    fig = go.Figure()
    fig.add_traces(plotting.plot_borders(provinces))
    fig.add_traces(plotting.plot_pipelines(pipe_coords))
    fig.add_traces(plotting.plot_points(points_coords))

    layout_map["mapbox"] = {"center": {"lon": -8.9, 'lat': 30.2}, 'zoom': 7,
                            'style': background, 'accesstoken': token}

    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'select+event'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)

    fig.update_layout(layout_map)
    return fig


def get_graphs(data, water_delivered, wwtp_data, desal_data, ag_lcoe, butane_data):
    data['water delivered'] = plotting.water_delivered_plot(water_delivered, 'Year', layout)
    data['unmet water'] = plotting.unmet_demand_plot(water_delivered, 'Year', layout)
    data['water supplied'] = plotting.water_supply_plot(water_delivered, 'Year', layout)
    data['energy demand'] = plotting.energy_demand_plot(water_delivered, wwtp_data, desal_data, 'Year', layout)
    data['depth to groundwater'] = plotting.wtd_plot(water_delivered, 'Date', layout)
    data['energy ag'] = plotting.energy_demand_ag(butane_data, layout)
    data['pv capacity'] = plotting.pv_installed_capacity(butane_data, layout)
    data['emissions ag'] = plotting.emissions_ag(butane_data, layout)
    # data = lcoe_plot(data, ag_lcoe)
    return data


@app.callback(
    [Output("current", "data"), Output('map', 'selectedData')],
    [
        Input("button-apply", "n_clicks"),
    ],
    [State('rate-wind', 'value'), State('rate-pv', 'value'), State('rate-grid', 'value'), State('rb-scenario', 'value'),
     State('climate-input', 'value'), State('butane-year', 'value'), State('pv-share', 'value')]
    # State("map-options", 'value'), State('cho-map-drop', 'value')]
)
def update_current_data(n_1, rate_wind, rate_pv, rate_grid, scenario, climate, butane_year, pv_share):
    water_delivered, wwtp_data, desal_data, butane_data = load_data(my_path, scenario, climate,
                                                                    butane_year, pv_share, 'all')
    # if map_op == 'cho-map':
    #     #     map = choroplethmap(provinces_chmap, provinces['Province'], label)
    #     # elif map_op == 'sch-map':
    #     #
    #     # else:
    #     #     map = {}

    # map = plot_map(water_delivered)
    map = {}
    data = {}
    graphs = get_graphs(data, water_delivered, wwtp_data, desal_data, None, butane_data)

    data_dict = dict(map=map, graphs=graphs, scenario=scenario, level=climate, pv_share=pv_share,
                     butane_year=butane_year)
    return data_dict, None


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
    Output("collapse", "is_open"),
    [Input("navbar-toggle", "n_clicks")],
    [State("collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    [Output("graphs", "children"), Output('resultsSubTitle', 'children')],
    [Input('map', 'selectedData'),
     Input('current-language', "modified_timestamp")],
    [State('current-language', 'data'), State('current', 'data')]
)
def update_results(selection, ts2, language, data_current):
    if data_current is None:
        raise PreventUpdate

    language_dic = get_language(language)

    colors = {'water': "#59C3C3", 'energy': "#F9ADA0", 'food': "#849E68"}
    data = {}

    if selection is None:
        name = 'Souss-Massa'
        data = data_current['graphs']

    elif selection['points'][0]['customdata'][0] in ['WWTP']:
        name = selection['points'][0]['customdata'][1]
        water_delivered, wwtp_data = load_data(my_path, data_current['scenario'],
                                               data_current['level'],
                                               data_current['butane_year'],
                                               data_current['pv_share'],
                                               ['results.gz', 'wwtp_data.gz'])
        df_tww = wwtp_data.loc[wwtp_data.point == name].groupby('Year')[['value', 'swpa_e']].sum().reset_index()
        df_tww.rename(columns={'value': 'water', 'swpa_e': 'energy'}, inplace=True)
        df_tww['water'] /= 1000000
        df_tww['energy'] /= 1000000
        data['water treated'] = plotting.wastewater_treated(df_tww, layout, [colors['water']])

        df_sww = water_delivered.loc[water_delivered['Supply point'] == name].groupby(['Year',
                                                                                     'Demand point'])['sswd'].sum().reset_index()
        df_sww.rename(columns={'sswd': 'water'}, inplace=True)
        df_tww['water'] /= 1000000
        data['water supplied'] = plotting.wastewater_supply(df_sww, layout)

        data['energy demand'] = plotting.wwt_energy(df_tww, layout, [colors['energy']])

    elif selection['points'][0]['customdata'][0] in ['Agriculture', 'Municipality']:
        name = selection['points'][0]['customdata'][1]

        water_delivered = load_data(my_path, data_current['scenario'],
                                    data_current['level'],
                                    data_current['butane_year'],
                                    data_current['pv_share'],
                                    ['results.gz'])

        df = water_delivered.loc[water_delivered['Demand point'] == name].groupby(['Year',
                                                                                   'Demand point',
                                                                                   'type'])['sswd'].sum().reset_index()
        df.rename(columns={'sswd': 'water'}, inplace=True)
        df['water'] /= 1000000
        data['water delivered'] = plotting.water_delivered(df, layout)

        data['unmet water'] = plotting.unmet_demand_plot(water_delivered.loc[water_delivered['Demand point'] == name], 'Year', layout)

    elif selection['points'][0]['customdata'][0] in ['Groundwater supply',
                                                     'Surfacewater withdrawal',
                                                     'Reservoir supply']:
        name = selection['points'][0]['customdata'][1]
        water_delivered = load_data(my_path, data_current['scenario'],
                                    data_current['level'],
                                    data_current['butane_year'],
                                    data_current['pv_share'],
                                    ['results.gz'])
        df = water_delivered.loc[water_delivered['Supply point'] == name]
        print(list(df))

        data['water supplied'] = plotting.water_delivered_plot(df, 'Year', layout)
        if selection['points'][0]['customdata'][0] == 'Groundwater supply':
            data['depth to groundwater'] = plotting.wtd_plot(df, 'Date', layout)
        data['energy demand'] = plotting.pumping_energy(df, layout)

    elif selection['points'][0]['customdata'][0] in ['Desalination plant']:
        name = selection['points'][0]['customdata'][1]
        desal_data = load_data(my_path, data_current['scenario'],
                                    data_current['level'],
                                    data_current['butane_year'],
                                    data_current['pv_share'],
                                    ['desal_data.gz'])
        df = desal_data.loc[desal_data['Supply point'] == name]

        data['water supplied'] = plotting.water_supply_plot(df, 'Year', layout, 'Demand point')
        data['energy demand'] = plotting.desal_energy(df, layout)

    else:
        raise PreventUpdate


    plots = []
    for key, value in data.items():
        value['layout']['title'] = language_dic['graphs'][key]

        plots.append(dcc.Graph(figure=value, config=dict(toImageButtonOptions=dict(
            format='png', filename=key, height=400,
            width=400, scale=2))))

    return plots, name


@app.callback(
    [Output("map", "figure"), Output("loading-1", "children")],
    [
        Input('current', 'modified_timestamp'),
    ]
)
def update_level_dropdown(ts):
    map = plot_map("light")
    return {}, dcc.Graph(figure=map, id="map",
                         config=dict(showSendToCloud=True,
                                     toImageButtonOptions=dict(format='png',
                                                               filename='map',
                                                               height=700,
                                                               width=700, scale=2)))


for info_id in info_ids:
    @app.callback(
        Output(f"{info_id}-modal", "is_open"),
        [Input(info_id, "n_clicks"), Input(f"{info_id}-close", "n_clicks")],
        [State(f"{info_id}-modal", "is_open")],
    )
    def toggle_popover(n_1, n_2, is_open):
        if n_1 or n_2:
            return not is_open
        return is_open

for modal_id in ['about', 'compare']:
    @app.callback(
        Output(f"{modal_id}-modal", "is_open"),
        [Input(modal_id, "n_clicks"), Input(f"{modal_id}-close", "n_clicks")],
        [State(f"{modal_id}-modal", "is_open")],
    )
    def toggle_popover(n_1, n_2, is_open):
        if n_1 or n_2:
            return not is_open
        return is_open


# @app.callback(
#     [Output('cho-map-drop', 'disabled'),
#      Output('cho-map-drop', 'value'),
#      # Output('cho-map-filter', 'disabled'),
#      # Output('cho-map-filter', 'value'),
#      Output("cho-map-collapse", "is_open"),
#      ],
#     [Input('map-options', 'value')],
# )
# def enable_choromap_options(value):
#     if value == 'cho-map':
#         return False, None, True
#     return True, None, False

@app.callback(
    [
        Output('rb-scenario', 'value'),
        Output('climate-input', 'value'),
        Output('rate-wind', 'value'),
        Output('rate-pv', 'value'),
        Output('rate-grid', 'value'),
        Output('butane-year', 'value'),
        Output('pv-share', 'value'),
        # Output('map-options', 'value')
        # Output('compare-options', 'value'),
    ],
    [Input('button-reset', 'n_clicks')],
)
def reset_output(n):
    return 'Reference', [], 0, 0, 0, 2050, 10


@app.callback(Output("compare-data", "data"),
              [Input("current-language", "modified_timestamp")],
              [State('current-language', 'data')])
def read_compare_data(ts, language):
    if language is None:
        raise PreventUpdate

    language_dic = get_language(language)

    # df_results = pd.read_csv('s3://souss-massa-project/compare_results.gz')
    # df_desal = pd.read_csv('s3://souss-massa-project/compare_desal.gz')
    # df_wwtp = pd.read_csv('s3://souss-massa-project/compare_wwtp.gz')

    df_results = pd.read_csv(os.path.join(my_path, 'compare_results.gz'))
    df_desal = pd.read_csv(os.path.join(my_path, 'compare_desal.gz'))
    df_wwtp = pd.read_csv(os.path.join(my_path, 'compare_wwtp.gz'))

    df_results = df_results.loc[df_results.Year >= 2020]
    df_desal = df_desal.loc[df_desal.Year >= 2020]
    df_wwtp = df_wwtp.loc[df_wwtp.Year >= 2020]

    fig_supply = plotting.water_supply_compare_plot(df_results, 'Year', language_dic['graphs']['water supply compare'])
    fig_energy = plotting.energy_demand_compare_plot(df_results, df_wwtp, df_desal, 'Year',
                                                     language_dic['graphs']['energy demand compare'])
    fig_unmet = plotting.unmet_demand_compare_plot(df_results, 'Year', language_dic['graphs']['unmet demand compare'])

    df_butane = pd.read_csv(os.path.join(my_path, 'compare_butane.gz'))

    df = df_butane.loc[df_butane.Year >= 2021].groupby(['butane_phaseout', 'pv_adoption'])[
        ['water_demand(m3)', 'energy_demand(KWh)', 'pv_elec(KWh)', 'grid_elec(KWh)', 'butane_cons(tonnes)',
         'butane_FARcost(mMAD)', 'PV_new_cap(MW)', 'reinv_cap(MW)', 'butane_emissions(MtCO2)', 'grid_emissions(MtCO2)',
         'Total_emissions(MtCO2)', 'butane_SUBSIDY(mMAD)', 'grid_cost(mMAD)', 'PV_Capex(mMAD)', 'pv_demand(KWh)',
         'butane_demand(KWh)', 'grid_demand(KWh)']].sum().reset_index()

    fig_butane_costs = plotting.total_costs_plot(df, language_dic['graphs']['total costs compare'])
    fig_butane_emissions = plotting.emissions_compare_plot(df_butane, language_dic['graphs']['emissions compare'])
    fig_butane_emissions_total = plotting.total_emissions_compare_plot(df,
                                                                       language_dic['graphs']['total emissions compare'])
    fig_emissions_costs = plotting.emisions_vs_costs(df_butane, language_dic['graphs']['emissions vs costs compare'])
    fig_energy_share = plotting.energy_resources_share_plot(df_butane, language_dic['graphs']['energy share compare'])

    return {'water supply': fig_supply, 'energy demand': fig_energy, 'unmet demand': fig_unmet,
            'butane costs': fig_butane_costs, 'butane emissions': fig_butane_emissions,
            'butane total emissions': fig_butane_emissions_total, 'emissions vs costs': fig_emissions_costs,
            'resources share': fig_energy_share}


@app.callback(Output("compare-content", "children"),
              [Input("compare-tabs", "active_tab"),
               Input('compare', 'n_clicks')],
              [State("compare-data", "data")])
def switch_tab(at, n, data):
    if data is None:
        raise PreventUpdate

    if at == "tab-1":
        return dcc.Graph(figure=data['water supply'],
                         config=dict(toImageButtonOptions=dict(format='png', filename='water_supply_all',
                                                               height=500, width=900, scale=2)))
    elif at == "tab-2":
        return dcc.Graph(figure=data['energy demand'],
                         config=dict(toImageButtonOptions=dict(format='png', filename='energy_demand_all',
                                                               height=500, width=900, scale=2)))
    elif at == "tab-3":
        return dcc.Graph(figure=data['unmet demand'],
                         config=dict(toImageButtonOptions=dict(format='png', filename='unmet_demand_ag_all',
                                                               height=500, width=900, scale=2)))
    elif at == "tab-4":
        return [dcc.Graph(figure=data['resources share'],
                          config=dict(toImageButtonOptions=dict(format='png', filename='energy_resources_share_ag',
                                                                height=600, width=900, scale=2))),
                dcc.Graph(figure=data['butane costs'],
                          config=dict(toImageButtonOptions=dict(format='png', filename='total_butane_system_cost',
                                                                height=500, width=900, scale=2))),
                dcc.Graph(figure=data['butane emissions'],
                          config=dict(toImageButtonOptions=dict(format='png', filename='annual_emissions_butane',
                                                                height=500, width=900, scale=2))),
                dcc.Graph(figure=data['butane total emissions'],
                          config=dict(toImageButtonOptions=dict(format='png', filename='total_emissions_butane',
                                                                height=500, width=900, scale=2))),

                dcc.Graph(figure=data['emissions vs costs'],
                          config=dict(toImageButtonOptions=dict(format='png', filename='emissions_vs_costs',
                                                                height=500, width=900, scale=2)))]


@app.callback(Output("download", "data"), [Input("button-download", "n_clicks_timestamp")])
def func(ts):
    if ts is None:
        raise PreventUpdate
    df_butane = pd.read_csv('s3://souss-massa-project/compare_butane.gz')
    return send_data_frame(df_butane.to_csv, "butan_phaseout.csv")


@app.callback(
    Output('current-language', 'data'),
    [Input(language, 'n_clicks_timestamp') for language in ['english', 'spanish', 'french']]
)
def current_language(n1, n2, n3):
    language_list = ['english', 'spanish', 'french']
    n_list = []
    for n in [n1, n2, n3]:
        if n is None:
            n_list.append(0)
        else:
            n_list.append(n)

    language_index = np.array(n_list).argmax()
    language = language_list[language_index]
    return language


@app.callback(
    [Output('title', "children"), Output('page-1', 'children'), Output('scenarios-title', "children"),
     Output('rb-scenario', "options"), Output('climate', "children"), Output('cost-reduction', "children"),
     Output('wind-title', 'children'), Output('pv-title', 'children'), Output('price-increase', 'children'),
     Output('grid-title', 'children'), Output('butane-title', 'children'), Output('butane-phaseout', 'children'),
     Output('pv-adoption-title', 'children'), Output('pv-adoption', 'children'),
     Output('button-reset', 'children'), Output('button-apply', 'children'),
     Output('language', 'label'), Output('about', 'children'),
     Output('resultsTitle', 'children'), Output('button-download', 'children'), Output('scenarioTitle', 'children')] +
    [Output(f'{i}-title', 'children') for i in info_ids] + [Output(f'{i}-body', 'children') for i in info_ids] +
    [Output(f'{i}-close', 'children') for i in info_ids] + [Output(f'about-{i}', 'children') for i in
                                                            ['title', 'body', 'close']] +
    [Output('compare', 'children'), Output('compare-close', 'children')] +
    [Output(f'tab-{i + 1}', 'label') for i in range(4)],
    # [Input(language, "n_clicks_timestamp") for language in ['english', 'spanish', 'french']] + [
    [Input('current-language', 'modified_timestamp'), Input('current', 'modified_timestamp')],
    [State('current-language', 'data'), State('current', 'data')]
)
def update_language(ts, ts2, language, data_current):
    if data_current is None:
        raise PreventUpdate
    language_dic = get_language(language)

    options = [
        {"label": language_dic["sidebar"]["scenarios"]["options"][0], "value": 'Reference'},
        {"label": language_dic['sidebar']['scenarios']['options'][1], "value": 'Desalination'},
        {"label": language_dic['sidebar']['scenarios']['options'][2], "value": 'Irrigation intensification'},
        {"label": language_dic['sidebar']['scenarios']['options'][3], "value": "Desalination Wastewater Reuse"},
        {"label": language_dic['sidebar']['scenarios']['options'][4], "value": "Reference Wastewater Reuse"}
    ]
    climate = data_current['level']
    if not climate:
        climate = ['historical trend']

    results_string = {
        'english': f" - {language_dic['results']['scenarios'][data_current['scenario']]} & {language_dic['results']['climate'][climate[0]]} scenario",
        'spanish': f" - Escenario de {language_dic['results']['scenarios'][data_current['scenario']]} & {language_dic['results']['climate'][climate[0]]}",
        'french': f" - {language_dic['results']['scenarios'][data_current['scenario']]} & {language_dic['results']['climate'][climate[0]]} scénario"
    }

    info_scenarios = [html.H6(language_dic['information']['scenarios']['body'][0]),
                      html.P(language_dic['information']['scenarios']['body'][1]),
                      html.H6(language_dic['information']['scenarios']['body'][2]),
                      html.P(language_dic['information']['scenarios']['body'][3]),
                      html.H6(language_dic['information']['scenarios']['body'][4]),
                      html.P(language_dic['information']['scenarios']['body'][5])]

    info_climate = [html.H6(language_dic['information']['climate']['body'][0]),
                    html.P(dcc.Markdown(f'''
                    * {language_dic['information']['climate']['body'][1]} 
                    * {language_dic['information']['climate']['body'][2]}
                    '''))]

    info_butane = [html.H6(language_dic['information']['butane']['body'][0]),
                   html.P(dcc.Markdown(f'''
                        * {language_dic['information']['butane']['body'][1]} 
                        * {language_dic['information']['butane']['body'][2]}
                        * {language_dic['information']['butane']['body'][3]}
                        '''))]

    info_pv_share = [html.H6(language_dic['information']['pv adoption']['body'][0]),
                     html.P(dcc.Markdown(f'''
                            * {language_dic['information']['pv adoption']['body'][1]} 
                            * {language_dic['information']['pv adoption']['body'][2]}
                            * {language_dic['information']['pv adoption']['body'][3]}
                            '''))]

    about_content = [html.P(language_dic['about']['body'][0]), html.P(language_dic['about']['body'][1]),
                     html.P(language_dic['about']['body'][2]), html.P(language_dic['about']['body'][3]),
                     dbc.Row([dbc.Col(html.A(html.Img(src='assets/kth.png', style={'height': '130px'}),
                                             href='https://www.energy.kth.se/energy-systems/about-the-division-of-energy-systems-1.937036'),
                                      width=3),
                              dbc.Col(html.A(html.Img(src='assets/sei.png', style={'height': '130px'}),
                                             href='https://www.sei.org/'), width=4),
                              dbc.Col(html.A(html.Img(src='assets/fao.png', style={'height': '130px'}),
                                             href='http://www.fao.org/home/en/'), width=2)], justify="center")
                     ]

    return language_dic['sidebar']['title'], language_dic['sidebar']['options'], \
           language_dic['sidebar']['scenarios']['title'], \
           options, language_dic['sidebar']['climate']['title'], language_dic['sidebar']['cost reduction']['title'], \
           language_dic['sidebar']['cost reduction']['wind'], language_dic['sidebar']['cost reduction']['pv'], \
           language_dic['sidebar']['price increase']['title'], language_dic['sidebar']['price increase']['grid'], \
           language_dic['sidebar']['butane']['title'], language_dic['sidebar']['butane']['year'], \
           language_dic['sidebar']['pv adoption']['title'], language_dic['sidebar']['pv adoption']['share'], \
           [html.I(className='fa fa-redo-alt'), f" {language_dic['sidebar']['reset']}"], \
           [html.I(className='fa fa-check-double'), f" {language_dic['sidebar']['apply']}"], \
           language_dic['language'], language_dic['about']['header'], \
           language_dic['results']['title'], \
           [html.I(className='fa fa-download'), f" {language_dic['results']['download']}"], \
           results_string[language], language_dic['information']['scenarios']['title'], \
           language_dic['information']['climate']['title'], language_dic['information']['cost']['title'], \
           language_dic['information']['price']['title'], language_dic['information']['butane']['title'], \
           language_dic['information']['pv adoption']['title'], info_scenarios, \
           info_climate, language_dic['information']['cost']['body'], \
           language_dic['information']['price']['body'], info_butane, info_pv_share, \
           language_dic['information']['close'], language_dic['information']['close'], \
           language_dic['information']['close'], language_dic['information']['close'], \
           language_dic['information']['close'], language_dic['information']['close'], \
           language_dic['about']['title'], about_content, language_dic['about']['close'], \
           [html.I(className='fa fa-chart-pie'), f" {language_dic['results']['compare']}"], \
           language_dic['compare']['close'], language_dic['compare']['tabs'][0], \
           language_dic['compare']['tabs'][1], language_dic['compare']['tabs'][2], \
           language_dic['compare']['tabs'][3]


if __name__ == "__main__":
    app.run_server(debug=True)
