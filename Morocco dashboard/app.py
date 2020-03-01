import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import json
import pandas as pd
import numpy as np
import geopandas as gpd
import os.path
import itertools
import boto3

client = boto3.client('s3')           
my_path = os.path.abspath(os.path.dirname(__file__))
# my_path = ''
spatial_data = os.path.join(my_path, 'spatial_data')

provinces = gpd.read_file(os.path.join(spatial_data, 'Admin', 'Provinces.gpkg'))
cropland = pd.read_csv(os.path.join(spatial_data, 'cropland.csv'))
demand_points = gpd.read_file(os.path.join(spatial_data, 'Demand_points.gpkg'))
demand_points.loc[demand_points['type']=='Catchment', 'type'] = 'Agriculture'
demand_points.loc[demand_points['type']=='Demand site', 'type'] = 'Municipality'
supply_points = gpd.read_file(os.path.join(spatial_data, 'Supply_points.gpkg'))
supply_points.loc[supply_points['type']=='Other supply', 'type'] = 'Desalination plant'
pipelines = gpd.read_file(os.path.join(spatial_data, 'Pipelines.gpkg'))
WebMercator = 4326

for gdf in [demand_points, supply_points, pipelines, provinces]:
    gdf.to_crs(epsg=WebMercator, inplace=True)


def load_data(scenario, climate, level):
    init_year = 2020
    data_folder = os.path.join(my_path, 'data')
    if not climate:
        climate = ['Trend']
    data = os.path.join(data_folder, scenario, climate[0])
    lcoe = os.path.join(data_folder, scenario, climate[0], level)

    water_delivered = pd.read_csv(os.path.join(data, 'results.csv'))
    # water_required = pd.read_csv(os.path.join(data, 'Water_requirements.csv'))
    ag_lcoe = pd.read_csv(os.path.join(lcoe, 'lcoe.csv'))
    wwtp_data = pd.read_csv(os.path.join(data, 'wwtp_data.csv'))
    desal_data = pd.read_csv(os.path.join(data, 'desal_data.csv'))
    desal_data['swpa_e'] = desal_data['sswd'] * 3.31

    return water_delivered.loc[water_delivered.Year>=init_year], ag_lcoe.loc[ag_lcoe.Year>=init_year], \
           wwtp_data.loc[wwtp_data.Year>=init_year], desal_data.loc[desal_data.Year>=init_year]

# def load_data(scenario, climate, level):
    # init_year = 2020
    # data_folder = 's3://souss-massa-project/data'
    # if not climate:
        # climate = ['Trend']
    # data = f'{data_folder}/{scenario}/{climate[0]}'
    # lcoe = f'{data_folder}/{scenario}/{climate[0]}/{level}'

    # water_delivered = pd.read_csv(f'{data}/results.csv')
    # ag_lcoe = pd.read_csv(f'{lcoe}/lcoe.csv', compression='zip')
    # wwtp_data = pd.read_csv(f'{data}/wwtp_data.csv')
    # desal_data = pd.read_csv(f'{data}/desal_data.csv')
    # desal_data['swpa_e'] = desal_data['sswd'] * 3.31

    # return water_delivered.loc[water_delivered.Year>=init_year], ag_lcoe.loc[ag_lcoe.Year>=init_year], \
           # wwtp_data.loc[wwtp_data.Year>=init_year], desal_data.loc[desal_data.Year>=init_year]


button_color = 'primary'
info_ids = []


def title_info(title, info_id, info_title, modal_content):
    info_ids.append(info_id)
    return dbc.Row([dbc.Col(html.H6(title), width=10),
                    dbc.Col(html.Button(html.Span(className="fa fa-info-circle"), className='info', id=info_id),
                            width=2),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(info_title),
                            dbc.ModalBody(modal_content),
                            dbc.ModalFooter(
                                dbc.Button("Close", id=f"{info_id}-close", className="ml-auto")
                            ),
                        ],
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
    xaxis={'tickformat': 'd'},
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
        dbc.Col([html.H3("Souss-Massa"), html.H6('NEXUS model', style={'color': 'gray'})]),
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
        title_info(title='Select scenario', info_id='scenario-info', info_title='Scenario information',
                   modal_content='This is the content of the modal'),
        html.Div(
            dbc.RadioItems(
                id="rb-scenario",
                options=[
                    {"label": "Reference", "value": 'Reference'},
                    {"label": "Desalinated water", "value": 'Desalination'},
                    {"label": "Irrigation Intensification", "value": 'Irrigation intensification'},
                ],
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
                        value=['Climate Change'],
                        id="climate-input",
                        switch=True,
                    ),
                    width=2,
                ),
                dbc.Col(title_info(title='Increased droughts events', info_id='climate-info',
                                   info_title='Climate change variables information',
                                   modal_content='This is the content of the modal')),
            ],
            no_gutters=True,
        )
    ],
    className='options'
)

technology_options = html.Div(
    [
        title_info(title='Select the cost reduction for year 2040', info_id='technology-info',
                   info_title='Technology cost reduction information', modal_content='This is the content of the modal'),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-wind"), html.Label('Wind power', style={'padding': '0 0 0 10px'})], style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=0.7,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'}, 0.7: {'label': '70%'}},
                                   included=False, id='rate-wind'),
                        style={'margin-top': '5px'}
                        ),
            ],
            no_gutters=True,
            style={'margin-top': '12px', 'margin-bottom': '10px'}
        ),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-solar-panel"), html.Label('Solar panels', style={'padding': '0 0 0 10px'})], style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=0.7,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'}, 0.7: {'label': '70%'}},
                                   included=False, id='rate-pv'),
                        style={'margin-top': '5px'}
                        ),
            ],
            no_gutters=True,
        ),
    ],
    className='options'
)

grid_options = html.Div(
    [
        title_info(title='Select the price increase for year 2040', info_id='grid-info',
                   info_title='Electricity price reduction information', modal_content='This is the content of the modal'),
        dbc.Row(
            [
                dbc.Col([html.I(className="fa fa-plug"), html.Label('Grid electricity', style={'padding': '0 0 0 10px'})], style={'font-size': '14px', 'color': 'gray'}, width=5),
                dbc.Col(dcc.Slider(min=0, max=0.7, value=0, step=0.7,
                                   marks={0: {'label': '0%'}, 0.3: {'label': '30%'}, 0.5: {'label': '50%'},
                                          0.7: {'label': '70%'}},
                                   included=False, id='rate-grid'),
                        style={'margin-top': '5px'}
                        ),
            ],
            no_gutters=True, style={'margin-top': '12px'}
        ),
    ],
    className='options'
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

map_options = html.Div(
    [
        title_info(title='Map visualization options', info_id='map-info',
                   info_title='Map visualization info', modal_content='This is the content of the modal'),
        html.Div(
            dbc.RadioItems(
                id="map-options",
                options=[
                    {"label": "System schematic", "value": 'sch-map'},
                    {"label": "Choropleth map", "value": 'cho-map'},
                ],
                value='sch-map',
                className='checklist-selected-style',
            ),
        ),
        dbc.Collapse(dbc.Card(dbc.CardBody(
            [
                dcc.Dropdown(id='cho-map-drop',
                             options=[
                                 {"label": "Cropland density", "value": 'Cropland density'},
                                 {"label": "Water delivered", "value": 'Water delivered'},
                                 {"label": "Energy demand", "value": 'Energy demand'},
                             ],
                             value='Cropland density',
                             placeholder='Select variable...',
                             clearable=False,
                             style={'marginBottom': '1em'}
                             ),
                dcc.Dropdown(id='cho-map-filter',
                             options=[
                                 {"label": "Province", "value": 'prov'},
                                 {"label": "Irrigation district", "value": 'irr'},
                             ],
                             value='crop',
                             placeholder='Filter by...',
                             clearable=False
                             ),
            ])),
            id="cho-map-collapse",
        ),
    ],
    className='options'
)

compare_scenarios = html.Div(
    [
        title_info(title='Select scenarios to display', info_id='compare-info',
                   info_title='Compare scenarios info', modal_content='This is the content of the modal'),
        html.Div(
            dcc.Dropdown(
                id="compare-options",
                options=[
                    {"label": "Current", "value": 'current'},
                    {"label": "Scenario 1", "value": 's1'},
                    {"label": "Scenario 2", "value": 's2'},
                ],
                value='current',
                clearable=False,
                searchable=False,
                multi=True,
                # className='checklist-selected-style',
            ),
        ),
    ],
    className='options'
)

scenario_tools = html.Div(
    [
        scenario_options,
        html.Hr(),
        climate_options,
        html.Hr(),
        technology_options,
        html.Hr(),
        grid_options,
        # html.Hr(),
        # unit_options,
        html.Hr(),
        map_options,
        html.Hr(),
        compare_scenarios
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
                dbc.NavItem(dbc.NavLink("Scenario", active=True, href="#", id='page-1', className='tabs')),
                dbc.NavItem(dbc.NavLink("Visualisation", href="#", id='page-2', className='tabs')),
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
        dbc.Button([html.I(className='fa fa-download'), " Download"], color=button_color, className="mr-1",
                   style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-download'),
    ],
    align='center',
    justify="center",
    className='footer',
)

map = dcc.Graph(id="map",
                className='col-lg-7 col-md-7 col-sm-12 col-xs-12',
                config=dict(showSendToCloud=True, toImageButtonOptions=dict(format='png', filename='map', height=700,
                                                                            width=700, scale=2)))

graphs = html.Div([results_header,
                   html.Div(dbc.Col(id='graphs'), id='graphs-container'),
                   footer_results
                   ],
                  id='results-container',
                  className='col-lg-5 col-md-5 col-sm-12 col-xs-12')

content = html.Div([map, graphs], id="page-content")

app.layout = html.Div([dcc.Store(id='current'), sidebar, content])

# Helper funtions
with open(os.path.join(spatial_data, 'Admin', 'Provinces.geojson')) as response:
    provinces_chmap = json.load(response)
    for feature in provinces_chmap['features']:
        # print(feature)
        feature['id'] = feature['properties']['Province']


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


def scatterpointmap(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data):
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
        lon=[point.coords[0][0] for point in provinces.boundary.centroid],
        lat=[point.coords[0][1] for point in provinces.boundary.centroid],
        text=provinces.Province,
        hoverinfo='skip',
    )]

    df_pipelines = pipelines

    # df_pipelines = pipelines.groupby('index').agg({'pipeline': 'first',
    #                                                'segment_length_m': 'first',
    #                                                'geometry': 'first'}).reset_index()
    # dff_pipeflow = pl_flow.groupby(['Year', 'pipeline']).agg({'water_use': lambda value: sum(value) / 1000000,
    #                                                           'SWPA_E_': lambda value: sum(value) / 1000000,
    #                                                           'pipeline_length': 'first'}).reset_index()
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
        # customdata=list(itertools.chain.from_iterable([[{
        #     'water': [dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].water_use,
        #               dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].Year],
        #     'energy': [dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].SWPA_E_,
        #                dff_pipeflow.loc[dff_pipeflow.pipeline == pipe].Year],
        #     'type': "pipeline"} for i in range(len(list(line.coords.xy[0])))] + [{}] for line, pipe in
        #                                                zip(df_pipelines.geometry, df_pipelines.pipeline)])),
    )]

    df_demand = demand_points.groupby('point').agg({'type': 'first',
                                                    # 'elevation_m': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_delivered = water_delivered.groupby(['Year', 'type', 'Demand point'])[['sswd','swpa_e']].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    # dff_required = water_required.groupby(['Year', 'type', 'Demand point'])['value'].sum() / 1000000
    # dff_required = dff_required.reset_index()

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
                     # names[1]: [dff_required.loc[(dff_required.point == point)].value,
                     #            dff_required.loc[(dff_required.point == point)].Year],
                     'type': 'demand',
                     } for point in df_demand.loc[df_demand.type == type].point],
    ) for type in sorted(df_demand['type'].unique())]

    df_supply = supply_points.groupby('point').agg({'type': 'first',
                                                    # 'wtd_m': 'first',
                                                    # 'elevation_m': 'first',
                                                    'geometry': 'first'}).reset_index()
    dff_delivered = water_delivered.groupby(['Year', 'type', 'Supply point'])[['sswd','swpa_e']].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    # dff_gw = gw_pumped.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    # dff_wwtp = wwtp_data.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    # dff_desal = desal_data.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
    # dff_surface = pl_flow.groupby(['Year', 'type', 'point'])[['water_use']].sum() / 1000000
    # dff_surface.rename(columns={'water_use': 'value'}, inplace=True)
    # dff_surface = dff_surface.reset_index()
    # dff_gw = dff_gw.reset_index()
    # dff_wwtp = dff_wwtp.reset_index()
    # dff_desal = dff_desal.reset_index()
    # dff = dff_gw.append([dff_wwtp, dff_surface, dff_desal], sort=False, ignore_index=True)
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


def plot_map(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data):
    layout_map = layout.copy()

    data = scatterpointmap(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data)

    layout_map["mapbox"] = {"center": {"lon": -8.9, 'lat': 30.2}, 'zoom': 7,
                            'style': "light", 'accesstoken': token}
    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'event+select'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)
    map = dict(data=data, layout=layout_map)
    return map

def water_delivered_plot(data, water_delivered):
    name = 'Water delivered (Mm3)'
    dff_delivered = water_delivered.copy()
    dff_delivered.loc[dff_delivered['type'].str.contains('Agriculture'), 'category'] = 'Agriculture'
    dff_delivered.loc[dff_delivered['type'].str.contains('Domestic'), 'category'] = 'Domestic'
    dff_delivered.loc[dff_delivered['type'].str.contains('Aquifer'), 'category'] = 'Aquifer recharge'
    dff_delivered.loc[dff_delivered['type'].str.contains('SW|Aquifer'), 'type'] = 'Surface water'
    dff_delivered.loc[dff_delivered['type'].str.contains('GW'), 'type'] = 'Groundwater'
    dff_delivered.loc[dff_delivered['type'].str.contains('DS'), 'type'] = 'Desalinated water'
    dff_delivered = dff_delivered.loc[water_delivered['type'] != 'Transmission Pipeline']
    dff_delivered = dff_delivered.groupby(['Year', 'category', 'type'])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()

    dff_delivered_cat = dff_delivered.groupby(['Year', 'category']).sswd.sum()
    dff_delivered['share'] = dff_delivered.sswd / dff_delivered.set_index(['Year', 'category']).index.map(dff_delivered_cat)

    for df, name in zip([dff_delivered_cat.reset_index()], [name]):
        data[name] = [{'x': df.loc[df['category'] == category].Year,
                       'y': df.loc[df['category'] == category].sswd,
                       'name': category,
                       'stackgroup': 'one',
                       'mode': 'lines',
                       # 'text': dff_unmet.loc[dff_unmet['type'] == type].value,
                       'text': ["<br>".join([f'{row[1]["type"]}: {round(row[1]["share"]*100,2)}%' for row in group.iterrows()]) \
                                for year, group in dff_delivered.loc[dff_delivered['category'] == category].groupby('Year')],
                       'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}' +
                                        '<br><b>Supply</b>: %{text}'
                       } for category in sorted(df['category'].unique())]
    return data

def energy_demand_plot(data, water_delivered, wwtp_data, desal_data):
    emission_factor = 1.76
    dff_energy = water_delivered.copy()
    dff_energy.loc[dff_energy['type'].str.contains('GW'), 'type'] = 'Groundwater pumping'
    dff_energy.loc[dff_energy['type'].str.contains('SW|Pipeline'), 'type'] = 'Surface water conveyance'
    dff_energy.loc[dff_energy['type'].str.contains('DS'), 'type'] = 'Desalinated water conveyance'
    dff_energy = dff_energy.groupby(['Year', 'type'])['swpa_e'].sum() / 1000000
    dff_energy = dff_energy.reset_index()
    wwtp_data['type'] = 'Wastewater treatment'
    desal_data['type'] = 'Desalination energy'
    for df in [wwtp_data, desal_data]:
        dff = df.groupby(['Year', 'type'])['swpa_e'].sum() / 1000000
        dff = dff.reset_index()
        dff_energy = dff_energy.append(dff, sort=False)

    for df, name in zip([dff_energy], ['Energy demand (GWh)']):
        data[name] = [{'x': df.loc[df['type'] == type].Year,
                       'y': df.loc[df['type'] == type].swpa_e,
                       'name': type,
                       'stackgroup': 'one',
                       'mode': 'lines',
                       'text': df.loc[df['type'] == type].swpa_e * emission_factor / 1000,
                       'hovertemplate': '<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}' +
                                        '<br><b>Emissions</b>: %{text: 0.2f} MtCO2'
                       } for type in sorted(df['type'].unique())]
    return data

def lcoe_plot(data, ag_lcoe):
    dff_lcoe = ag_lcoe.dropna().copy()
    dff_lcoe['area_m2'] = dff_lcoe['Demand point'].map(cropland.area_m2)

    dff_lcoe.loc[dff_lcoe['least_cost_technology'].str.contains('Grid'), 'least_cost_technology'] = 'Grid electricity'
    count = len(dff_lcoe['Demand point'].unique())
    # total = dff_lcoe.loc[dff_lcoe['Year']==2020, 'area_m2'].sum()
    for df, name in zip([dff_lcoe], ['Least-cost power technology for Agriculture']):
        data[name] = [{'x': dff_lcoe.loc[dff_lcoe['least_cost_technology']==tech].groupby('Year')['least_cost_technology'].count().index,
                       'y': dff_lcoe.loc[dff_lcoe['least_cost_technology']==tech].groupby('Year')['least_cost_technology'].count() / count,
                       'name': tech,
                       'type': 'bar',
                       'hovertemplate': '<b>Value</b>: %{y:,.2%}' + '<br><b>Year</b>: %{x}'
                       } for tech in sorted(dff_lcoe['least_cost_technology'].unique())]
    return data

def get_graphs(data, water_delivered, water_required, wwtp_data, desal_data, ag_lcoe):
    data = water_delivered_plot(data, water_delivered)
    data = energy_demand_plot(data, water_delivered, wwtp_data, desal_data)
    data = lcoe_plot(data, ag_lcoe)
    # dff_required = water_required.groupby(['Year', 'type'])['value'].sum() / 1000000
    # dff_required = dff_required.reset_index()

    # dff_unmet = dff_required.copy()
    # dff_unmet['value'] = (dff_unmet.value - dff_unmet.set_index(['Year', 'type']).index.map(
    #     dff_delivered.set_index(['Year', 'type']).value)) / dff_unmet.value * 100
    # dff_unmet.loc[dff_unmet['value'] < 0.001, 'value'] = 0
    return data


@app.callback(
    [Output("current", "data"), Output('map', 'selectedData')],
    [
        Input("button-apply", "n_clicks"),
    ],
    [State('rate-wind', 'value'), State('rate-pv', 'value'), State('rate-grid', 'value'), State('rb-scenario', 'value'),
     State('climate-input', 'value'), State("map-options", 'value'), State('cho-map-drop', 'value')]
)
def update_current_data(n_1, rate_wind, rate_pv, rate_grid, scenario, climate, map_op, label):
    water_delivered, ag_lcoe, wwtp_data, desal_data = load_data(scenario, climate, f'W{rate_wind}_PV{rate_pv}_Grid{-rate_grid}')
    if map_op == 'cho-map':
        map = choroplethmap(provinces_chmap, provinces['Province'], label)
    elif map_op == 'sch-map':
        map = plot_map(water_delivered, None, None, None, wwtp_data, desal_data)
    else:
        map = {}
    data = {}
    graphs = get_graphs(data, water_delivered, None, wwtp_data, desal_data, ag_lcoe)

    data_dict = dict(map=map, graphs=graphs, scenario=scenario, level=climate, eff_end=None, eff_init=None)
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
    [Output("graphs", "children"), Output('resultsTitle', 'children')],
    [Input('map', 'selectedData'),
     Input("sidebar-toggle", "n_clicks"),
     Input('scenarioTitle', 'children'),
     Input('current', 'data')]
)
def update_results(selection, n_1, title, data_current):
    if data_current is None:
        raise PreventUpdate
    names = ['Water delivered (Mm3)']
    colors = {'water': "#59C3C3", 'energy': "#F9ADA0", 'food': "#849E68"}
    data = {}
    if selection is None:
        name = 'Souss-Massa basin'
        data = data_current['graphs']

    elif selection['points'][0]['customdata']['type'] in ['demand', 'supply']:
        name = selection['points'][0]['text']
        name_key = selection['points'][0]['customdata']['type']
        name_dict = {'demand': {'water': 'Water demand (Mm3)'}, 'supply': {'water': 'Water supplied (Mm3)'}}
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
            name_dict = {'supply': {'energy': 'Energy demand (GWh)'}}
            df = pd.DataFrame(selection['points'][0]['customdata']['water'][0])
            data[name_dict[name_key]['energy']] = [{'x': df.groupby('Year').swpa_e.sum().index,
                                                    'y': df.groupby('Year').swpa_e.sum(),
                                                    'fill': 'tozeroy',
                                                    'showlegend': False,
                                                    'mode': 'lines',
                                                    'line': {'color': colors['energy']}
                                                    }]

    # elif selection['points'][0]['customdata']['type'] in ['supply', 'WWTP', 'pipeline', 'Desalination']:
    #     name = selection['points'][0]['text']
    #     name_key = selection['points'][0]['customdata']['type']
    #     name_dict = {'supply': {'water': 'Water supplied (Mm3)', 'energy': 'Energy for pumping (GWh)'},
    #                  'WWTP': {'water': 'Treated wastewater (Mm3)', 'energy': 'Energy for treatment (GWh)'},
    #                  'pipeline': {'water': 'Water conveyed (Mm3)', 'energy': 'Energy for pumping (GWh)'},
    #                  'Desalination': {'water': 'Water desalinated (Mm3)', 'energy': 'Energy for desalination (GWh)'}}
    #
    #     df = pd.DataFrame(selection['points'][0]['customdata']['water'][0])
    #
    #     data[name_dict[name_key]['water']] = [{'x': group.Year,
    #                                            'y': group.sswd,
    #                                            'stackgroup': 'one',
    #                                            'mode': 'lines',
    #                                            'name': type,
    #                                            # 'line': {'color': colors['water']}
    #                                            } for type, group in df.groupby('type')]
    #     # data[name_dict[name_key]['energy']] = [{'x': selection['points'][0]['customdata']['energy'][1],
    #     #                                         'y': np.array(selection['points'][0]['customdata']['energy'][0]),
    #     #                                         'fill': 'tozeroy', 'mode': 'lines', 'showlegend': False,
    #     #                                         'line': {'color': colors['energy']}}]

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
        layout_plot['barmode'] = 'stack'
        layout_plot['font'] = dict(size=10, color="#7f7f7f")
        if 'Least-cost' in key:
            layout_plot['yaxis'] = {'tickformat': ',.0%','range': [0, 1]}
            layout_plot['annotations'] = [dict(xref='paper',
                                          yref='paper',
                                          x=0, y=1.12,
                                          showarrow=False,
                                          text='Percentage of Agricultural sites per technology')]
        plots.append(dcc.Graph(figure=dict(data=value, layout=layout_plot), config=dict(toImageButtonOptions=dict(
            format='png', filename=key, height=500,
            width=700, scale=2))))
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


@app.callback(
    Output('scenarioTitle', 'children'),
    [
        # Input("button-apply", "n_clicks"),
        Input("sidebar-toggle", "n_clicks"),
        Input("current", "data"),
    ],
)
def update_level_dropdown(n_2, data):
    if data is None:
        raise PreventUpdate

    scenario = data['scenario']
    level = data['level']
    if not level:
        level = ['historical trend']
    name = f' - {scenario} {level[0]} scenario'

    return name


@app.callback(
    Output("map", "figure"),
    [
        Input("current", "data"),
        # Input('scenarioTitle', 'children'),
    ],
)
def update_level_dropdown(data):
    if data is None:
        raise PreventUpdate
    map = data['map']
    return map


@app.callback(
    [Output('page-1', 'active'), Output('page-2', 'active'),
     Output('tools', 'hidden'), Output('visual-tools', 'hidden')],
    [Input('page-1', 'n_clicks_timestamp'), Input('page-2', 'n_clicks_timestamp')],
)
def selected_tab(n_1, n_2):
    if n_2 is None:
        n_1 = 1
        n_2 = 0
        tools = False
        visual = True
    elif n_1 is None:
        n_1 = 0
        n_2 = 1
        tools = True
        visual = False
    if n_1 > n_2:
        state_1 = True
        state_2 = False
        tools = False
        visual = True
    else:
        state_1 = False
        state_2 = True
        tools = True
        visual = False
    return state_1, state_2, tools, visual


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


@app.callback(
    [Output('cho-map-drop', 'disabled'),
     Output('cho-map-drop', 'value'),
     Output('cho-map-filter', 'disabled'),
     Output('cho-map-filter', 'value'),
     Output("cho-map-collapse", "is_open"),
     ],
    [Input('map-options', 'value')],
)
def enable_choromap_options(value):
    if value == 'cho-map':
        return False, None, False, None, True
    return True, None, True, None, False



if __name__ == "__main__":
    app.run_server(debug=True)
