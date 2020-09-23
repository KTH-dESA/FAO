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

from scripts import plots

# import boto3

# client = boto3.client('s3')
# my_path = os.path.abspath(os.path.dirname(__file__))
my_path = os.path.join('..', 'Morocco model', 'test dash results')
spatial_data = os.path.join(my_path, 'spatial_data')

# language = open("assets/english.yaml")
# language = yaml.load(language, Loader=yaml.FullLoader)

provinces = gpd.read_file(os.path.join(spatial_data, 'Admin', 'Provinces.gpkg'))
cropland = pd.read_csv(os.path.join(spatial_data, 'cropland.csv'))
demand_points = gpd.read_file(os.path.join(spatial_data, 'Demand_points.gpkg'))
demand_points.loc[demand_points['type'] == 'Catchment', 'type'] = 'Agriculture'
demand_points.loc[demand_points['type'] == 'Demand site', 'type'] = 'Municipality'
supply_points = gpd.read_file(os.path.join(spatial_data, 'Supply_points.gpkg'))
supply_points.loc[supply_points['type'] == 'Other supply', 'type'] = 'Desalination plant'
pipelines = gpd.read_file(os.path.join(spatial_data, 'Pipelines.gpkg'))
WebMercator = 4326

centroids = provinces.centroid

for gdf in [demand_points, supply_points, pipelines]:
    gdf.to_crs(epsg=WebMercator, inplace=True)


def load_data(scenario, climate, level, phaseout_year, pv_level):
    init_year = 2020
    butane_scenario = f'phaseout_{phaseout_year}' if phaseout_year!=2050 else 'phaseout_None'
    data_folder = os.path.join(my_path, 'data')
    if not climate:
        climate = ['Trend']
    data = os.path.join(data_folder, scenario, climate[0])
    lcoe = os.path.join(data_folder, scenario, climate[0], level)

    water_delivered = pd.read_csv(os.path.join(data, 'results.gz'))
    # water_required = pd.read_csv(os.path.join(data, 'Water_requirements.csv'))
    # ag_lcoe = pd.read_csv(os.path.join(lcoe, 'lcoe.gz'))
    ag_lcoe = pd.DataFrame({'Year': []})
    wwtp_data = pd.read_csv(os.path.join(data, 'wwtp_data.gz'))
    desal_data = pd.read_csv(os.path.join(data, 'desal_data.gz'))
    desal_data['swpa_e'] = desal_data['sswd'] * 3.31
    butane_data = pd.read_csv(os.path.join(my_path, 'Butane_calculations', butane_scenario,
                                           f'{pv_level}_PV', 'summary_results.gz'))

    return water_delivered.loc[water_delivered.Year >= init_year], ag_lcoe.loc[ag_lcoe.Year >= init_year], \
           wwtp_data.loc[wwtp_data.Year >= init_year], desal_data.loc[desal_data.Year >= init_year], \
           butane_data.loc[butane_data.Year >= init_year]


# def load_data(scenario, climate, level):
#     init_year = 2020
#     data_folder = 's3://souss-massa-project/data'
#     if not climate:
#         climate = ['Trend']
#     data = f'{data_folder}/{scenario}/{climate[0]}'
#     lcoe = f'{data_folder}/{scenario}/{climate[0]}/{level}'
#
#     water_delivered = pd.read_csv(f'{data}/results.csv')
#     ag_lcoe = pd.read_csv(f'{lcoe}/lcoe.csv', compression='zip')
#     wwtp_data = pd.read_csv(f'{data}/wwtp_data.csv')
#     desal_data = pd.read_csv(f'{data}/desal_data.csv')
#     desal_data['swpa_e'] = desal_data['sswd'] * 3.31
#
#     return water_delivered.loc[water_delivered.Year>=init_year], ag_lcoe.loc[ag_lcoe.Year>=init_year], \
#            wwtp_data.loc[wwtp_data.Year>=init_year], desal_data.loc[desal_data.Year>=init_year]


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


token = "pk.eyJ1IjoiY2FtaWxvcmciLCJhIjoiY2p1bTl0MGpkMjgyYjQ0b2E0anRibWJ1MSJ9.GhUUGD6gok1d36lvP17CQQ"

layout = dict(
    autosize=True,
    automargin=True,
    margin=dict(l=30, r=20, b=50, t=100),
    hovermode="closest",
    plot_bgcolor="#fff",
    paper_bgcolor="#fff",
    legend=dict(font=dict(size=10), orientation="h"),
    showlegend=True,
    hoverlabel={'align': 'left'},
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

# unit_options = html.Div(
#     [
#         title_info(title='Select display units', info_id='units-info',
#                    info_title='Units info', modal_content='This is the content of the modal'),
#         dbc.InputGroup(
#             [
#                 dbc.InputGroupAddon("Water demand", addon_type="prepend"),
#                 dbc.Select(id='water-units',
#                            options=[
#                                {"label": "Mm3", "value": "Mm3"},
#                                {"label": "m3", "value": "m3"},
#                            ],
#                            value='Mm3',
#                            className='form-control'
#                            ),
#             ],
#             size="sm",
#             style={'marginBottom': '1em'}
#         ),
#         dbc.InputGroup(
#             [
#                 dbc.InputGroupAddon("Energy demand", addon_type="prepend"),
#                 dbc.Select(id='energy-units',
#                            options=[
#                                {"label": "GWh", "value": "GWh"},
#                                {"label": "MWh", "value": "MWh"},
#                                {"label": "Wh", "value": "Wh"},
#                                {"label": "PJ", "value": "PJ"},
#                            ],
#                            value='GWh'
#                            ),
#             ],
#             size="sm"
#         ),
#     ],
#     className='options'
# )

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
#
# compare_scenarios = html.Div(
#     [
#         title_info(title='Select scenarios to display', info_id='compare-info',
#                    info_title='Compare scenarios info', modal_content='This is the content of the modal'),
#         html.Div(
#             dcc.Dropdown(
#                 id="compare-options",
#                 options=[
#                     {"label": "Current", "value": 'current'},
#                     {"label": "Scenario 1", "value": 's1'},
#                     {"label": "Scenario 2", "value": 's2'},
#                 ],
#                 value='current',
#                 clearable=False,
#                 searchable=False,
#                 multi=True,
#                 # className='checklist-selected-style',
#             ),
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
                               dcc.Loading(dbc.ModalBody(id="compare-content"))]),
                dbc.ModalFooter(
                    dbc.Button("Close", id="compare-close", className="ml-auto")
                ), ],
            size="xl",
            id="compare-modal",
        ),
        dbc.Button([html.I(className='fa fa-download'), " Download"], color=button_color, className="mr-1",
                   style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-download'),
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


def scatterpointmap(water_delivered):
    data = [dict(
        type="scattermapbox",
        mode='lines',
        line=dict(color='rgb(80,100,80)', width=1),
        fill='toself',
        fillcolor='rgba(80,100,80,0.1)',
        layer="below",
        lon=list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(
            [[list(l.coords.xy[0]) + [None] for l in line.boundary] for line in provinces.geometry])))),
        lat=list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(
            [[list(l.coords.xy[1]) + [None] for l in line.boundary] for line in provinces.geometry])))),
        hoverinfo='skip',
        showlegend=False,
    )]

    data += [dict(
        type='scattermapbox',
        mode='text',
        name='Province',
        textfont=dict(color='rgba(80,80,80,0.8)'),
        # lon=[point.coords[0][0] for point in provinces.centroid],
        # lat=[point.coords[0][1] for point in provinces.centroid],
        lon=[point.coords[0][0] for point in centroids],
        lat=[point.coords[0][1] for point in centroids],
        text=provinces.Province,
        hoverinfo='skip',
    )]

    df_pipelines = pipelines

    data += [dict(
        type="scattermapbox",
        mode='lines+markers',
        marker=dict(opacity=0),
        lon=list(itertools.chain.from_iterable([list(line.coords.xy[0]) + [None] for line in df_pipelines.geometry])),
        lat=list(itertools.chain.from_iterable([list(line.coords.xy[1]) + [None] for line in df_pipelines.geometry])),
        text=list(itertools.chain.from_iterable([len(list(line.coords.xy[0])) * [pipe] + [None] for line, pipe in
                                                 zip(df_pipelines.geometry, df_pipelines.diversion)])),
        hoverinfo='skip',
        name='Pipeline',
        showlegend=True,
    )]

    df_demand = demand_points.groupby('point').agg({'type': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_delivered = water_delivered.groupby(['Year', 'type', 'Demand point']).agg({'sswd': lambda x: sum(x) / 1000000,
                                                                                   'swpa_e': lambda x: sum(x) / 1000000,
                                                                                   'unmet_demand_year': 'mean'})
    dff_delivered = dff_delivered.reset_index()
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'

    data += [dict(
        type="scattermapbox",
        lon=[point.x for point in df_demand.loc[df_demand['type'] == type].geometry],
        lat=[point.y for point in df_demand.loc[df_demand['type'] == type].geometry],
        text=df_demand.loc[df_demand['type'] == type].point,
        ids=df_demand.loc[df_demand['type'] == type].point,
        hoverinfo='name+text',
        name=type,
        showlegend=True,
        customdata=[{'water': [dff_delivered.loc[(dff_delivered['Demand point'] == point)].to_dict()],
                     'unmet': [dff_delivered.loc[dff_delivered['Demand point'] == point].to_dict()],
                     'type': 'demand',
                     } for point in df_demand.loc[df_demand.type == type].point],
    ) for type in sorted(df_demand['type'].unique())]

    df_supply = supply_points.groupby('point').agg({'type': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_delivered = water_delivered.groupby(['Year', 'type', 'Supply point'])[['sswd', 'swpa_e']].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    data += [dict(
        type="scattermapbox",
        lon=[point.x for point in df_supply.loc[df_supply['type'] == type].geometry],
        lat=[point.y for point in df_supply.loc[df_supply['type'] == type].geometry],
        text=df_supply.loc[df_supply['type'] == type].point,
        ids=df_supply.loc[df_supply['type'] == type].point,
        hoverinfo='name+text',
        name=type,
        showlegend=True,
        customdata=[{'water': [dff_delivered.loc[(dff_delivered['Supply point'] == point)].to_dict()],
                     'energy': [dff_delivered.loc[(dff_delivered['Supply point'] == point)].to_dict()],
                     'type': 'supply',
                     } for point in df_supply.loc[df_supply.type == type].point],
    ) for type in sorted(df_supply['type'].unique())]

    return data


def plot_map(water_delivered):
    layout_map = layout.copy()

    data = scatterpointmap(water_delivered)

    layout_map["mapbox"] = {"center": {"lon": -8.9, 'lat': 30.2}, 'zoom': 7,
                            'style': "light", 'accesstoken': token}
    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'event+select'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)
    map = dict(data=data, layout=layout_map)
    return map


def water_delivered_plot(data, water_delivered, time_frame):
    name = 'water delivered'
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
    data[name] = [{'x': df.loc[df['category'] == category][time_frame],
                   'y': df.loc[df['category'] == category].sswd,
                   'name': category,
                   'stackgroup': 'one',
                   'mode': 'lines',
                   'text': ["<br>".join(
                       [f'{row[1]["type"]}: {round(row[1]["share"] * 100, 2)}%' for row in group.iterrows()]) \
                       for year, group in
                       dff_delivered.loc[dff_delivered['category'] == category].groupby(time_frame)],
                   'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}' +
                                    '<br><b>Supply</b>: %{text}'
                   } for category in sorted(df['category'].unique())]
    return data


def energy_demand_plot(data, water_delivered, wwtp_data, desal_data, time_frame):
    emission_factor = 1.76
    dff_energy = water_delivered.copy()
    dff_energy.loc[dff_energy['type'].str.contains('GW'), 'type'] = 'Groundwater pumping'
    dff_energy.loc[dff_energy['type'].str.contains('Pipeline'), 'type'] = 'Surface water conveyance'
    dff_energy.loc[dff_energy['type'].str.contains('DS'), 'type'] = 'Desalinated water conveyance'
    # dff_energy.loc[dff_energy['type'].str.contains('WWR'), 'type'] = 'Wastewater reuse conveyance'
    # dff_energy = dff_energy.loc[~dff_energy['type'].str.contains('WWR')]
    # dff_energy = dff_energy.loc[~dff_energy['type'].str.contains('SW')]
    dff_energy = dff_energy.loc[~dff_energy['Supply point'].str.contains('Complexe Aoulouz Mokhtar Soussi')]
    dff_energy = dff_energy.groupby([time_frame, 'type'])['swpa_e'].sum() / 1000000
    dff_energy = dff_energy.reset_index()
    wwtp_data['type'] = 'Wastewater treatment'
    desal_data['type'] = 'Desalination energy'
    for df in [wwtp_data, desal_data]:
        dff = df.groupby([time_frame, 'type'])['swpa_e'].sum() / 1000000
        dff = dff.reset_index()
        dff_energy = dff_energy.append(dff, sort=False)

    for t in dff_energy['type'].unique():
        if dff_energy.loc[dff_energy['type'] == t, 'swpa_e'].sum() == 0:
            dff_energy = dff_energy.loc[dff_energy['type'] != t]

    name = 'energy demand'
    df = dff_energy
    data[name] = [{'x': df.loc[df['type'] == type][time_frame],
                   'y': df.loc[df['type'] == type].swpa_e,
                   'name': type,
                   'stackgroup': 'one',
                   'mode': 'lines',
                   'text': df.loc[df['type'] == type].swpa_e * emission_factor / 1000,
                   'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}'
                   # '<br><b>Emissions</b>: %{text: 0.2f} MtCO2'
                   } for type in sorted(df['type'].unique())]
    return data


def lcoe_plot(data, ag_lcoe):
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


def unmet_demand_plot(data, water_delivered, time_frame):
    name = 'unmet water'
    dff_unmet = water_delivered.copy()
    dff_unmet.loc[dff_unmet['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_unmet.loc[dff_unmet['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_unmet = dff_unmet.loc[~dff_unmet['type'].str.contains('Aquifer')]
    dff_unmet = dff_unmet.loc[dff_unmet['type'] != 'Transmission Pipeline']
    water_req_year = \
        dff_unmet.groupby(['Year', 'Date', 'Demand point', 'category'])['water_required'].mean().reset_index().groupby(
            [time_frame, 'category'])['water_required'].sum()
    unment_demand = 1 - (dff_unmet.groupby([time_frame, 'category'])['sswd'].sum() /
                         water_req_year)

    df = unment_demand.reset_index()
    data[name] = [{'x': df.loc[df['category'] == category][time_frame],
                   'y': df.loc[df['category'] == category].iloc[:, 2],
                   'name': category,
                   'mode': 'lines',
                   'hovertemplate': '<b>Value</b>: %{y:,.2%}' + '<br><b>Year</b>: %{x}'
                   } for category in sorted(df['category'].unique())]
    return data


def wtd_plot(data, water_delivered, time_frame):
    name = 'depth to groundwater'
    dff_wtd = water_delivered.copy()
    dff_wtd = dff_wtd.loc[dff_wtd['type'].str.contains('GW')]
    wtd = dff_wtd.groupby([time_frame, 'Supply point'])['wtd'].mean().reset_index()
    df = wtd
    data[name] = [{'x': df.loc[df['Supply point'] == point][time_frame],
                   'y': df.loc[df['Supply point'] == point].wtd,
                   'name': point,
                   'mode': 'lines',
                   'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}'
                   } for point in sorted(df['Supply point'].unique())]
    return data


def water_supply_plot(data, water_delivered, time_frame):
    name = 'water supplied'
    dff_delivered = water_delivered.copy()
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    dff_delivered.loc[dff_delivered['type'].str.contains('WWR'), 'type'] = 'Reused Wastewater'
    dff_delivered = dff_delivered.loc[dff_delivered['type'] != 'Transmission Pipeline']
    dff_delivered = dff_delivered.groupby([time_frame, 'type'])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()

    df = dff_delivered.reset_index()
    data[name] = [{'x': df.loc[df['type'] == type][time_frame],
                   'y': df.loc[df['type'] == type].sswd,
                   'name': type,
                   'stackgroup': 'one',
                   'mode': 'lines',
                   'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}'
                   } for type in sorted(df['type'].unique())]
    return data


def get_graphs(data, water_delivered, wwtp_data, desal_data, ag_lcoe, butane_data):
    data = water_delivered_plot(data, water_delivered, 'Year')
    data = unmet_demand_plot(data, water_delivered, 'Year')
    data = water_supply_plot(data, water_delivered, 'Year')
    data = energy_demand_plot(data, water_delivered, wwtp_data, desal_data, 'Year')
    data = wtd_plot(data, water_delivered, 'Date')
    data['energy ag'] = plots.energy_demand_ag(butane_data)
    data['pv capacity'] = plots.pv_installed_capacity(butane_data)
    data['emissions ag'] = plots.emissions_ag(butane_data)
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
    water_delivered, ag_lcoe, wwtp_data, desal_data, butane_data = load_data(scenario, climate,
                                                                             f'W{rate_wind}_PV{rate_pv}_Grid{-rate_grid}',
                                                                             butane_year, pv_share)
    # if map_op == 'cho-map':
    #     #     map = choroplethmap(provinces_chmap, provinces['Province'], label)
    #     # elif map_op == 'sch-map':
    #     #
    #     # else:
    #     #     map = {}

    map = plot_map(water_delivered)
    data = {}
    graphs = get_graphs(data, water_delivered, wwtp_data, desal_data, ag_lcoe, butane_data)

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
     # Input("sidebar-toggle", "n_clicks"),
     Input('current', 'modified_timestamp'), Input('current-language', "modified_timestamp")],
    [State('current-language', 'data'), State('current', 'data')]
)
def update_results(selection, ts, ts2, language, data_current):
    if data_current is None:
        raise PreventUpdate

    language_dic = get_language(language)

    colors = {'water': "#59C3C3", 'energy': "#F9ADA0", 'food': "#849E68"}
    data = {}
    if selection is None:
        name = 'Souss-Massa'
        data = data_current['graphs']

    elif 'WWTP' in selection['points'][0]['text']:
        name = selection['points'][0]['text']
        water_delivered, ag_lcoe, wwtp_data, desal_data, butane_data = load_data(data_current['scenario'], data_current['level'],
                                                                    'W{0}_PV{0}_Grid{-0}', data_current['butane_year'],
                                                                    data_current['pv_share'])
        df = wwtp_data.loc[wwtp_data.point == name].groupby('Year')['value', 'swpa_e'].sum()
        data['water treated'] = [{'x': df.index,
                                  'y': df.value,
                                  'stackgroup': 'one',
                                  'mode': 'lines',
                                  'name': selection['points'][0]['text'],
                                  'line': {'color': colors['water']}
                                  }]

        dff = water_delivered.loc[water_delivered['Supply point'] == name].groupby(['Year', 'Demand point'])[
            'sswd'].sum().reset_index()
        data['water supplied'] = [{'x': dff.loc[dff['Demand point'] == point].Year,
                                   'y': dff.loc[dff['Demand point'] == point].sswd,
                                   'stackgroup': 'one',
                                   'mode': 'lines',
                                   'name': point,
                                   } for point in dff['Demand point'].unique()]

        data['energy demand'] = [{'x': df.index,
                                  'y': df.swpa_e,
                                  'stackgroup': 'one',
                                  'mode': 'lines',
                                  'name': name,
                                  'line': {'color': colors['energy']}
                                  }]


    elif selection['points'][0]['customdata']['type'] in ['demand', 'supply']:
        name = selection['points'][0]['text']
        name_key = selection['points'][0]['customdata']['type']
        name_dict = {'demand': {'water': 'water delivered'}, 'supply': {'water': 'water supplied'}}
        df = pd.DataFrame(selection['points'][0]['customdata']['water'][0])
        data[name_dict[name_key]['water']] = [{'x': group.Year,
                                               'y': group.sswd,
                                               'stackgroup': 'one',
                                               'mode': 'lines',
                                               'name': type,
                                               } for type, group in df.groupby('type')]

        if selection['points'][0]['customdata']['type'] == 'supply':
            name = selection['points'][0]['text']
            name_key = selection['points'][0]['customdata']['type']
            name_dict = {'supply': {'energy': 'energy demand'}}
            df = pd.DataFrame(selection['points'][0]['customdata']['water'][0])
            data[name_dict[name_key]['energy']] = [{'x': df.groupby('Year').swpa_e.sum().index,
                                                    'y': df.groupby('Year').swpa_e.sum(),
                                                    'fill': 'tozeroy',
                                                    'showlegend': False,
                                                    'mode': 'lines',
                                                    'line': {'color': colors['energy']}
                                                    }]

        if selection['points'][0]['customdata']['type'] == 'demand':
            name = selection['points'][0]['text']
            name_key = selection['points'][0]['customdata']['type']
            name_dict = {'demand': {'unmet': 'unmet water'}}
            df = pd.DataFrame(selection['points'][0]['customdata']['unmet'][0])
            data[name_dict[name_key]['unmet']] = [{'x': df.groupby('Year').unmet_demand_year.mean().index,
                                                   'y': df.groupby('Year').unmet_demand_year.mean(),
                                                   'fill': 'tozeroy',
                                                   'showlegend': False,
                                                   'mode': 'lines',
                                                   'line': {'color': colors['water']}
                                                   }]

    plots = []
    for key, value in data.items():
        layout_plot = layout.copy()
        layout_plot['title'] = dict(text=language_dic['graphs'][key])
        layout_plot['height'] = 400
        layout_plot['barmode'] = 'stack'
        layout_plot['font'] = dict(size=10, color="#7f7f7f")

        if selection is None:
            if 'least-cost' in key:
                layout_plot['yaxis'] = {'tickformat': ',.0%', 'range': [0, 1]}
                layout_plot['annotations'] = dict(xref='paper',
                                                   yref='paper',
                                                   x=0, y=1.12,
                                                   showarrow=False,
                                                   text='Percentage of Agricultural sites per technology')
                plots.append(
                    dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                        format='png', filename=key, height=300,
                        width=400, scale=2))))
            elif 'unmet' in key:
                layout_plot['yaxis'] = {'tickformat': '%'}
                plots.append(
                    dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                        format='png', filename=key, height=300,
                        width=400, scale=2))))
            elif ('delivered' in key) or ('supplied' in key):
                # layout_plot['yaxis'] = {'range': [0, 1100]}
                plots.append(
                    dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                        format='png', filename=key, height=300,
                        width=400, scale=2))))
            elif 'energy demand' in key:
                # layout_plot['yaxis'] = {'range': [0, 1600]}
                plots.append(
                    dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                        format='png', filename=key, height=400,
                        width=400, scale=2))))
            elif 'depth' in key:
                layout_plot['yaxis'] = {'range': [200, 0]}
                plots.append(
                    dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                        format='png', filename=key, height=400,
                        width=400, scale=2))))
            else:
                plots.append(
                    dcc.Graph(figure=dict(data=value['data'], layout=layout_plot),
                              config=dict(toImageButtonOptions=dict(
                                  format='png', filename=key, height=300,
                                  width=400, scale=2))))
        else:
            if 'unmet' in key:
                layout_plot['yaxis'] = {'tickformat': '%', 'range': [0, 1]}
            plots.append(dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
                format='png', filename=key, height=400,
                width=400, scale=2))))

    return plots, name


# @app.callback(
#     [Output('drop-level', 'options'), Output('drop-level', 'value'), Output('drop-level', 'disabled')],
#     [Input('rb-scenario', 'value')]
# )
# def update_level_dropdown(scenario):
#     level_dict = {'Reference': {'Select...': 'level_1'},
#                   'Improve AG eff': {'by 10 percent': 'level_1',
#                                      'by 20 percent': 'level_2'},
#                   'New Resources': {'Select...': 'level_1'},
#                   'Reduce NRW': {'to 40 percent': 'level_1',
#                                  'to 20 percent': 'level_2'}}
#     options = [{"label": key, 'value': value} for key, value in level_dict[scenario].items()]
#     disable = False
#     if len(level_dict[scenario].keys()) == 1:
#         disable = True
#     return options, 'level_1', disable


# @app.callback(
#     Output('scenarioTitle', 'children'),
#     [
#         # Input("button-apply", "n_clicks"),
#         # Input("sidebar-toggle", "n_clicks"),
#         Input('current', 'modified_timestamp')
#     ],
#     [State("current", "data")]
# )
# def update_level_dropdown(ts, data):
#     if data is None:
#         raise PreventUpdate
#     if ts is None:
#         raise PreventUpdate
#
#     scenario = data['scenario']
#     level = data['level']
#     if not level:
#         level = ['historical trend']
#     name = f' - {scenario} {level[0]} scenario'
#
#     return name


@app.callback(
    Output("map", "figure"),
    [
        # Input("sidebar-toggle", "n_clicks"),
        Input('current', 'modified_timestamp'),
    ],
    [State("current", "data")]
)
def update_level_dropdown(ts, data):
    if data is None:
        raise PreventUpdate
    map = data['map']
    return map


# @app.callback(
#     [
#         Output('page-1', 'active'),
#          Output('page-2', 'active'),
#          Output('tools', 'hidden'),
#          Output('visual-tools', 'hidden')
#     ],
#     [
#         Input('page-1', 'n_clicks_timestamp'),
#         Input('page-2', 'n_clicks_timestamp')
#      ],
# )
# def selected_tab(n_1, n_2):
#     if n_2 is None:
#         n_1 = 1
#         n_2 = 0
#         tools = False
#         visual = True
#     elif n_1 is None:
#         n_1 = 0
#         n_2 = 1
#         tools = True
#         visual = False
#     if n_1 > n_2:
#         state_1 = True
#         state_2 = False
#         tools = False
#         visual = True
#     else:
#         state_1 = False
#         state_2 = True
#         tools = True
#         visual = False
#     return state_1, state_2, tools, visual


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
              [Input("compare", "n_clicks"), Input("current-language", "modified_timestamp")],
              [State('current-language', 'data')])
def read_compare_data(n, ts, language):
    language_dic = get_language(language)
    # General results
    scenarios = ['Reference', 'Desalination', 'Reference Wastewater Reuse', 'Desalination Wastewater Reuse']
    climates = ['Climate Change']
    df_results, df_desal, df_wwtp = plots.load_results_data(my_path, scenarios, climates)

    df_results = df_results.loc[df_results.Year >= 2020]
    df_desal = df_desal.loc[df_desal.Year >= 2020]
    df_wwtp = df_wwtp.loc[df_wwtp.Year >= 2020]

    fig_supply = plots.water_supply_compare_plot(df_results, 'Year', language_dic['graphs']['water supply compare'])
    fig_energy = plots.energy_demand_compare_plot(df_results, df_wwtp, df_desal, 'Year',
                                                  language_dic['graphs']['energy demand compare'])
    fig_unmet = plots.unmet_demand_compare_plot(df_results, 'Year', language_dic['graphs']['unmet demand compare'])

    # Butane results
    scenarios = [None, 2040, 2030]
    pv_levels = [10, 20, 50]
    df_butane = plots.load_butane_results(my_path, scenarios, pv_levels)

    df = df_butane.loc[df_butane.Year >= 2021].groupby(['butane_phaseout', 'pv_adoption'])[
        ['water_demand(m3)', 'energy_demand(KWh)', 'pv_elec(KWh)', 'grid_elec(KWh)', 'butane_cons(tonnes)',
         'butane_FARcost(mMAD)', 'PV_new_cap(MW)', 'reinv_cap(MW)', 'butane_emissions(MtCO2)', 'grid_emissions(MtCO2)',
         'Total_emissions(MtCO2)', 'butane_SUBSIDY(mMAD)', 'grid_cost(mMAD)', 'PV_Capex(mMAD)', 'pv_demand(KWh)',
         'butane_demand(KWh)', 'grid_demand(KWh)']].sum().reset_index()

    fig_butane_costs = plots.total_costs_plot(df, language_dic['graphs']['total costs compare'])
    fig_butane_emissions = plots.emissions_compare_plot(df_butane, language_dic['graphs']['emissions compare'])
    fig_butane_emissions_total = plots.total_emissions_compare_plot(df, language_dic['graphs']['total emissions compare'])
    fig_emissions_costs = plots.emisions_vs_costs(df_butane, language_dic['graphs']['emissions vs costs compare'])
    fig_energy_share = plots.energy_resources_share_plot(df_butane, language_dic['graphs']['energy share compare'])

    return {'water supply': fig_supply, 'energy demand': fig_energy, 'unmet demand': fig_unmet,
            'butane costs': fig_butane_costs, 'butane emissions': fig_butane_emissions,
            'butane total emissions': fig_butane_emissions_total, 'emissions vs costs': fig_emissions_costs,
            'resources share': fig_energy_share}


@app.callback(Output("compare-content", "children"),
              [Input("compare-tabs", "active_tab"), Input('compare-data', 'modified_timestamp')],
              [State("compare-data", "data")])
def switch_tab(at, ts, data):
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
                                                                height=500, width=900, scale=2)))
                ]

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
        {"label": language_dic['sidebar']['scenarios']['options'][3], "value": "Desalination wastewater reuse"},
        {"label": language_dic['sidebar']['scenarios']['options'][4], "value": "Reference wastewater reuse"}
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
           language_dic['results']['compare'], language_dic['compare']['close'], language_dic['compare']['tabs'][0], \
           language_dic['compare']['tabs'][1], language_dic['compare']['tabs'][2], language_dic['compare']['tabs'][3]


if __name__ == "__main__":
    app.run_server(debug=True)
