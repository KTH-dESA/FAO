import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame
import plotly.graph_objects as go
import plotly.io as pio
from dash_extensions.snippets import send_file
import json
import pandas as pd

import geopandas as gpd
import os.path

from scripts import plotting

my_path = os.path.abspath(os.path.dirname(__file__))
# my_path = ''
spatial_data = os.path.join(my_path, 'spatial_data')
pio.templates.default = "plotly_white"

with open(os.path.join(spatial_data, 'Admin', 'governorates.geojson')) as response:
    governorates = json.load(response)
with open(os.path.join(spatial_data, 'Admin', 'borders.geojson')) as response:
    borders = json.load(response)
# demand_points = gpd.read_file(os.path.join(spatial_data, 'Demand_points.gpkg'))  # TEST: changed from geojson
# supply_points = gpd.read_file(os.path.join(spatial_data, 'Supply_points.gpkg'))
# pipelines = gpd.read_file(os.path.join(spatial_data, 'Pipelines.gpkg'))
WebMercator = 4326

# for gdf in [demand_points, supply_points, pipelines]:
#     gdf.to_crs(epsg=WebMercator, inplace=True)

points_coords = pd.read_csv(os.path.join(spatial_data, 'points_coords.csv'))
pipe_coords = pd.read_csv(os.path.join(spatial_data, 'pipe_coords.csv'))

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
                        size="lg"
                    )])


token = "pk.eyJ1IjoiY2FtaWxvcmciLCJhIjoiY2p1bTl0MGpkMjgyYjQ0b2E0anRibWJ1MSJ9.GhUUGD6gok1d36lvP17CQQ"

layout = dict(
    autosize=True,
    height=350,
    margin=dict(l=0, r=20, b=50, t=100),
    hovermode="closest",
    plot_bgcolor="#fff",
    paper_bgcolor="#fff",
    legend=dict(font=dict(size=10), orientation="h", title=''),
    xaxis={'title': ''},
    showlegend=True,
    hoverlabel={'align': 'left'},
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
        dbc.Col([html.H3("Jordan"), html.H6('NEXUS model', style={'color': 'gray'})]),
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
                   modal_content=[
                       html.P('Four scenarios have been analysed to date in order to explore nexus interactions and '
                              'the impact any given nexus solution or policy measure targeted to one of the systems '
                              '(i.e. water, energy or agriculture), may have upon the other systems. '
                              'The scenarios evaluated so far are a reference (Business as Usual), an improved '
                              'agricultural efficiency and a new water resources scenario '
                              '(i.e. New desalinated water), and are run for 30 years from 2020 to 2050.'),
                       html.H6('Reference scenario'),
                       html.P('This scenario takes a Business as Usual approach where the main trends '
                              '(in terms of demand, supply and growth) are kept unchanged. It assumes that '
                              'domestic demands will increase over time, with refugees staying (but no new '
                              'refugees coming), and agriculture and industry not growing over time.'),
                       html.H6('Improve agricultural efficiency scenario'),
                       html.P('This scenario considers 10 and 20 percent increase in irrigation efficiency '
                              'above current levels by year 2050.'),
                       html.H6('Reduce non-revenue water scenario'),
                       html.P('The WEAP model was set up by MWI to include a representation of non-revenue '
                              'water losses in each governorate. Non-revenue water losses (i.e. water lost '
                              'to leakage, under-billing, and theft) accounts for about half of all water in '
                              'Jordan, which makes it a critical factor in improving water services in Jordan.'
                              ' This scenario assumes a reduction of total non-revenue water by 20% and 40% '
                              'levels, to year 2050. The reduction is set as a goal for each municipality to '
                              'achieve by year 2050.'),
                       html.H6('New water resources scenario'),
                       html.P('This scenario assumes the realization of the Red Sea-Dead Sea project and new '
                              'desalination plant Red-Dead (110 MCM/yr).')]),
        html.Div(
            dbc.RadioItems(
                id="rb-scenario",
                options=[
                    {"label": "Reference", "value": 'Reference'},
                    {"label": "Improved agricultural efficiency", "value": 'Improve AG eff'},
                    {"label": "New water resources", "value": 'New Resources'},
                    {"label": "Reduce non-revenue water", "value": 'Reduce NRW'},
                    {"label": "Increased water productivity", "value": 'Increased Water Productivity'},
                ],
                value='Reference',
                className='checklist-selected-style',
            ),
        ),
    ],
    className='options'
)

eto_options = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(
                    dbc.Checklist(
                        options=[
                            {"label": "", "value": 'Climate Change'},
                        ],
                        value=['Climate Change'],
                        id="eto-input",
                        switch=True,
                    ),
                    width=2,
                ),
                dbc.Col(title_info(title='Increase Evapotranspiration (climate change)', info_id='Eto-info',
                                   info_title='Evapotranspiration information',
                                   modal_content='Historical reference evapotranspiration (ETo) data indicates an '
                                                 'increasing trend in the country, which will potentially impact '
                                                 'hydrology and crop water requirements. To assess the impact of '
                                                 'climate change, all scenarios were evaluated under two conditions: '
                                                 'without increasing trend in ETo and an increasing trend applied to '
                                                 'hydrology and irrigation requirements.')),
            ],
            no_gutters=True,
        )
    ],
    className='options'
)

level_options = html.Div(
    [
        title_info(title='Select level of variable', info_id='level-info', info_title='Variable level information',
                   modal_content=[
                       dcc.Markdown('This parameter allows to select the level of improvement wanted in the scenario. '
                                    'It only applies for the **Improved agricultural efficiency** and the '
                                    '**Reduce non-revenue water** scenarios.')]),
        html.Div(
            dcc.Dropdown(
                id="drop-level",
                options=[
                    {"label": "Select...", "value": ''},
                    # {"label": "Level 2", "value": 2},
                ],
                value='',
                clearable=False,
            ),
        ),
    ],
    className='options',
)

energy_options = html.Div(
    [
        title_info(title='Select energy pumping efficiency', info_id='pump-eff-info',
                   info_title='Energy pumping efficiency information', modal_content='This parameter allows to select '
                                                                                     'the current level of energy '
                                                                                     'efficiency and set an improvement '
                                                                                     'goal for year 2050. Those values '
                                                                                     'will affect pumping energy requirements '
                                                                                     'for groundwater extraction and '
                                                                                     'surface water conveyance '
                                                                                     'through the pipeline system.'),
        dbc.Row([
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
    ],
    className='options'
)


save_scenario = html.Div(
    [
        title_info(title='Save scenario', info_id='save-info',
                   info_title='Save scenario info', modal_content='This is the content of the modal'),
        dbc.InputGroup(
            [
                dbc.InputGroupAddon("Save scenario", addon_type="prepend"),
                dbc.Input(id='save-scenario',
                          placeholder="Name the scenario",
                          className='form-control'
                          ),
            ],
            size="sm",
            style={'marginBottom': '1em'}
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
                value=['current'],
                multi=True,
                className='checklist-selected-style',
            ),
        ),
    ],
    className='options'
)

scenario_tools = html.Div(
    [
        scenario_options,
        html.Hr(),
        eto_options,
        html.Hr(),
        level_options,
        html.Hr(),
        energy_options
    ],
    id='tools',
)

visual_tools = html.Div(
    [
        save_scenario,
        html.Hr(),
        compare_scenarios,
    ],
    id='visual-tools',
    style=dict(display='none')
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
            # html.Hr(),
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
        dbc.ButtonGroup([
            # dbc.Button([html.I(className='fa fa-save'), " Save"], color=button_color, className="mr-1",
            #            style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-save'),
            # dbc.Button([html.I(className='fa fa-chart-pie'), " Compare"], color=button_color, outline=True, className="mr-1",
            #            style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-compare'),
            dbc.Button([html.I(className='fa fa-download'), Download(id="download"), " Download"], color=button_color,
                       className="mr-1",
                       style={'fontSize': '0.85rem', 'fontWeight': '600'}, id='button-download'),
        ])
    ],
    align='center',
    justify="around",
    className='footer',
)

map = html.Div(
    [html.Div([
        dbc.Button(html.I(className='fa fa-layer-group'), color=button_color,
                   outline=True, className='map-ctrl-bottom mr-1',
                   id="popover-map-target"),
        html.Div(
            [
                # dbc.PopoverHeader("Header"),
                dbc.PopoverBody([
                    html.H6('Map visualization'),
                    dbc.RadioItems(
                        id="map-selection",
                        options=[
                            {"label": "Schematic", "value": 'schematic'},
                            {"label": "Governorates", "value": 'governorates'},
                        ],
                        value='schematic',
                        className='checklist-selected-style',
                    ),
                    html.Hr(),
                    html.H6('Map background'),
                    dbc.RadioItems(
                        id="map-background",
                        options=[
                            {"label": "White", "value": 'white-bg'},
                            {"label": "Light", "value": 'light'},
                            {"label": "Outdoors", "value": 'outdoors'},
                            {"label": "Basic", "value": 'basic'},
                            {"label": "Dark", "value": 'dark'},
                            {"label": "Satellite-streets", "value": 'satellite-streets'},
                        ],
                        value='light',
                        className='checklist-selected-style',
                    ),
                ]),
            ],
            className="popover-map",
            tabIndex='-1',
            id="popover-map",
            # style=dict(display='none'),
            # target="popover-map-target",
            # placement='top-start',
            # hide_arrow=True,
        )],
        className='map-controls',
    ),
        dcc.Loading(
            id="loading-map",
            type="default",
            children=dcc.Graph(id="map",
                         config=dict(showSendToCloud=True,
                                     toImageButtonOptions=dict(format='png',
                                                               filename='map',
                                                               height=700,
                                                               width=700, scale=2)))
        )
    ],
    id='map-div',
    className='col-xl-7 col-lg-12',
)

graphs = html.Div([results_header,
                   html.Div(dbc.Col(dcc.Loading(id='graphs', style=dict(height='10px'))), id='graphs-container'),
                   footer_results
                   ],
                  id='results-container',
                  className='col-xl-5 col-lg-12')

content = html.Div([map, graphs], id="page-content")

app.layout = html.Div([dcc.Store(id='current'), dcc.Store(id='compare', data={}), sidebar, content])


# Helper funtions
def plot_map(background, map_type):
    if background == 'white-bg':
        admin = governorates
    else:
        admin = borders

    layout_map = {}

    if map_type == 'schematic':
        fig = go.Figure()
        fig.add_traces(plotting.plot_borders(admin))
        fig.add_traces(plotting.plot_pipelines(pipe_coords))
        fig.add_traces(plotting.plot_points(points_coords))

    elif map_type == 'governorates':
        df = pd.DataFrame(governorates['features'])
        df['color'] = range(len(df['id']))
        fig = plotting.choroplet_map(governorates, df)

    layout_map["mapbox"] = {"center": {"lon": 36.8, 'lat': 31.2}, 'zoom': 6,
                            'style': background, 'accesstoken': token}

    layout_map["margin"] = {"r": 0, "t": 0, "l": 0, "b": 0}
    layout_map['clickmode'] = 'select+event'
    layout_map['legend'] = dict(font=dict(size=12), orientation="h", x=0, y=0)

    fig.update_layout(layout_map)
    return fig


def get_graphs(data, water_delivered, water_required, gw_pumped, pl_flow,
               wwtp_data, desal_data, crop_production):
    emission_factor = 0.643924449
    dff_delivered = water_delivered.groupby(['Year', 'type'])['value'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_required = water_required.groupby(['Year', 'type'])['value'].sum() / 1000000
    dff_required = dff_required.reset_index()

    dff_unmet = dff_required.copy()
    dff_unmet['value'] = (dff_unmet.value - dff_unmet.set_index(['Year', 'type']).index.map(
        dff_delivered.set_index(['Year', 'type']).value)) / dff_unmet.value

    df_production = crop_production.groupby(['Year', 'variable']).sum() / 1000000
    df_production.reset_index(inplace=True)

    dff_energy = pd.DataFrame()
    pl_flow['type'] = 'Water conveyance'

    for df in [wwtp_data, desal_data, gw_pumped, pl_flow]:
        dff = df.groupby(['Year', 'type'])['SWPA_E_'].sum() / (1000000)
        dff = dff.reset_index()
        dff_energy = dff_energy.append(dff, sort=False)
    dff_energy.rename(columns={'SWPA_E_': 'value'}, inplace=True)

    data['WaterDelivered'] = plotting.water_delivered(dff_delivered, layout, 'Water delivered (Mm3)')
    data['UnmetDemand'] = plotting.unmet_demand(dff_unmet.reset_index(), layout, 'Unmet water demand (%)')

    dff_wtd = gw_pumped.copy()
    dff_wtd['point'] = [x[1] for x in dff_wtd['point'].str.split('_')]
    dff_wtd = dff_wtd.groupby(['Year', 'point'])['wtd_m'].mean().reset_index()
    data['GWdepth'] = plotting.wtd_plot(dff_wtd, layout, 'Average depth to groundwater (mbgl)')

    dff = gw_pumped.copy()
    dff['point'] = [x[1] for x in dff['point'].str.split('_')]
    dff = dff.groupby(['Year', 'point'])['value'].sum() / 1000000
    dff = dff.reset_index()
    data['WaterExtraction'] = plotting.groundwater_extraction(dff, layout, 'Groundwater extraction (Mm3)')

    data['CropProduction'] = plotting.plot_production(df_production, layout,
                                                      'Annual cropland production by type (kton)',
                                                      'variable')

    data['EnergyDemand'] = plotting.energy_demand(dff_energy, layout, 'Energy demand (GWh)')

    return data


##### need to change the compare scenarios logic to the calback of the compare modal
#### also add logic to the save button and modal
@app.callback(
    [Output("current", "data"), Output('map', 'selectedData')],
    [
        Input("button-apply", "n_clicks"),
        # Input("popover-map-target", "n_clicks"),
    ],
    [State('pump-eff-init', 'value'), State('pump-eff-end', 'value'), State('rb-scenario', 'value'),
     State('eto-input', 'value'), State('drop-level', 'value')],
)
def update_current_data(n_1, eff_init, eff_end, scenario, eto, level):
    water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data, crop_production = plotting.load_data(
        my_path, scenario,
        eto, level, 'all')

    # map = plot_map(background, map_type)
    map = {}
    graphs = get_graphs({}, water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data, crop_production)
    data_dict = dict(gw_df=gw_pumped.to_dict(), map=map, graphs=graphs, scenario=scenario,
                     level=level, eto=eto, eff_end=eff_end, eff_init=eff_init)

    # compare_data[scenario_name] = dict(scenario=scenario, eto=eto, level=level, eff_end=eff_end, eff_init=eff_init)
    # else:
    #     water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data = load_data(compare_data[compare[0]]['scenario'],
    #                                                                                            compare_data[compare[0]]['eto'],
    #                                                                                            compare_data[compare[0]]['level'])
    #     map = plot_map(water_delivered, water_required, gw_pumped, pl_flow, wwtp_data, desal_data)
    #     graphs = {}
    #     data_dict = dict(map=map, graphs=graphs, scenario=scenario, level=level, eff_end=eff_end, eff_init=eff_init)
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
    [Input('map', 'selectedData'), Input('map-selection', 'value')],
    [State('current', 'data')],
    prevent_initial_call=True
)
def update_results(selection, map_type, data_current):
    if data_current is None:
        raise PreventUpdate
    colors = {'water': "#59C3C3", 'energy': "#F9ADA0", 'food': "#849E68"}
    data = {}

    if selection is None:
        if map_type == 'schematic':
            name = 'Jordan country'
            data = data_current['graphs']
        else:
            name = 'Jordan country'
            dff_delivered, crop_production = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                                                data_current['level'],
                                                                ['Water_delivered.gz', 'crop_production.gz'])

            df = dff_delivered.groupby(['Year', 'Governorate'])[['value']].sum() / 1000000
            df.reset_index(inplace=True)


            data['WaterDeliveredGov'] = plotting.plot_water_delivered_by_gov(df, layout,
                                                                             'Annual water delivered by Governorate (Mm3)',)

            df = crop_production.groupby(['Year', 'Governorate'])[['value']].sum() / 1000000
            df.reset_index(inplace=True)

            data['CropProduction'] = plotting.plot_production(df, layout,
                                                              'Annual cropland production by Governorate (kton)',
                                                              'Governorate')

            df = crop_production.groupby(['Governorate', 'Year', 'variable'])[['value']].sum() / 1000000
            df.reset_index(inplace=True)
            df = df.groupby(['Governorate', 'variable'])[['value']].mean().reset_index()
            data['CropProductionByType'] = plotting.plot_production_by_gov(df, layout,
                                                                                  'Annual average cropland production (kton)',
                                                                                  'variable', 'Governorate')

            # data['CropProductionByGovernorate'] = plotting.plot_production_by_gov(df, layout,
            #                                                                       'Total cropland production (Mton)',
            #                                                                       'Governorate', 'variable')

    elif selection['points'][0]['customdata'][0] in ['Municipality', 'Industry', 'Agriculture']:
        name = selection['points'][0]['customdata'][1]
        water_delivered, water_required = plotting.load_data(my_path, data_current['scenario'],
                                                             data_current['eto'], data_current['level'],
                                                             ['water_delivered.gz', 'water_requirements.gz'])

        dff, dff_unmet = plotting.get_demand_data(water_delivered, water_required, name)

        data[f'{name}WaterDemand'] = plotting.plot_water_delivered(dff, layout, 'Water demand (Mm3)')
        data[f'{name}UnmetDemand'] = plotting.plot_unmet_demand(dff_unmet, layout, 'Unmet demand (%)')

    elif selection['points'][0]['customdata'][0] in ['Groundwater supply']:
        name = selection['points'][0]['customdata'][1]
        gw_pumped = plotting.load_data(my_path, data_current['scenario'],
                                       data_current['eto'], data_current['level'],
                                       'Groundwater_pumping.gz')
        dff = gw_pumped.loc[gw_pumped['point'] == name]
        dff = dff.groupby(['Year', 'type', 'point']).agg({'value': lambda x: sum(x) / 1000000,
                                                          'SWPA_E_': lambda x: sum(x) / 1000000,
                                                          'wtd_m': 'mean'})
        dff = dff.reset_index()
        data['water'] = plotting.plot_water_supply(dff, [colors['water']], layout, 'Water supplied (Mm3)')

        dff_wtd = gw_pumped.loc[gw_pumped['point'] == name].copy()
        date = gw_pumped[['Year', 'Month']].copy()
        date['Day'] = 1
        dff_wtd['Date'] = pd.to_datetime(date)

        data[f'{name}GWdepth'] = plotting.plot_depth_groundwater(dff_wtd, layout, 'Average depth to groundwater (mbgl)')
        data[f'{name}EnergyDemand'] = plotting.plot_energy_for_pumping(dff, [colors['energy']], layout,
                                                                       'Energy for pumping (GWh)')

    elif selection['points'][0]['customdata'][0] in ['Wastewater plant']:
        name = selection['points'][0]['customdata'][1]
        wwtp_data = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                       data_current['level'], 'wwtp_data.gz')
        dff = wwtp_data.loc[wwtp_data['point'] == name]
        dff = dff.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
        dff = dff.reset_index()
        data[f'{name}TreatedWastewater'] = plotting.plot_water_supply(dff, [colors['water']], layout,
                                                                      'Treated wastewater (Mm3)')
        data[f'{name}EnergyForTreatment'] = plotting.plot_energy_for_pumping(dff, [colors['energy']],
                                                                             layout, 'Energy for treatment (GWh)')

    elif selection['points'][0]['customdata'][0] in ['Desalination']:
        name = selection['points'][0]['customdata'][1]
        desal_data = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                        data_current['level'], 'desal_data.gz')
        dff = desal_data.loc[desal_data['point'] == name]
        dff = dff.groupby(['Year', 'type', 'point'])[['value', 'SWPA_E_']].sum() / 1000000
        dff = dff.reset_index()
        data[f'{name}DesalWater'] = plotting.plot_water_supply(dff, [colors['water']],
                                                               layout, 'Water desalinated (Mm3)')
        data[f'{name}DesalEnergy'] = plotting.plot_energy_for_pumping(dff, [colors['energy']],
                                                                      layout, 'Energy for desalination (GWh)')

    elif selection['points'][0]['customdata'][0] in ['pipeline']:
        name = selection['points'][0]['customdata'][1]
        pl_flow = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                     data_current['level'], 'Pipelines_data.gz')
        dff = pl_flow.loc[pl_flow['pipeline'] == name]
        dff = dff.groupby(['Year', 'pipeline']).agg({'water_use': lambda value: sum(value) / 1000000,
                                                     'SWPA_E_': lambda value: sum(value) / 1000000,
                                                     'pipeline_length': 'first'}).reset_index()
        dff.rename(columns={'water_use': 'value'}, inplace=True)
        data[f'{name}WaterConveyed'] = plotting.plot_water_supply(dff, [colors['water']],
                                                                  layout, 'Water conveyed (Mm3)')
        data[f'{name}PumpingEnergy'] = plotting.plot_energy_for_pumping(dff, [colors['energy']],
                                                                        layout, 'Energy for pumping (GWh)')


    elif selection['points'][0]['customdata'][0] in ['River/pipeline supply']:
        name = selection['points'][0]['customdata'][1]
        pl_flow = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                     data_current['level'], 'Pipelines_data.gz')
        dff = pl_flow.loc[pl_flow['point'] == name]
        dff = dff.groupby(['Year', 'type', 'point'])[['water_use']].sum() / 1000000
        dff.rename(columns={'water_use': 'value'}, inplace=True)
        dff = dff.reset_index()
        data[f'{name}WaterSupplied'] = plotting.plot_water_supply(dff, [colors['water']], layout,
                                                                  'Water supplied (Mm3)')

    elif selection['points'][0]['customdata'][0] in ['Other supply']:
        name = selection['points'][0]['customdata'][1]
        water_delivered = plotting.load_data(my_path, data_current['scenario'],
                                             data_current['eto'], data_current['level'],
                                             ['Water_delivered.gz'])
        dff = water_delivered.loc[water_delivered['point'] == name]
        dff = dff.groupby(['Year', 'type', 'point'])[['value']].sum() / 1000000
        dff.reset_index(inplace=True)
        data[f'{name}WaterSupplied'] = plotting.plot_water_supply(dff, [colors['water']],
                                                                  layout, 'Water conveyed (Mm3)')

    else:
        name = selection['points'][0]['customdata'][0]

        dff_delivered, crop_production = plotting.load_data(my_path, data_current['scenario'], data_current['eto'],
                                             data_current['level'], ['Water_delivered.gz', 'crop_production.gz'])

        df = dff_delivered.loc[dff_delivered['Governorate']==name]
        df = df.groupby(['Year', 'type'])['value'].sum() / 1000000
        df = df.reset_index()

        data[f'{name}WaterDelivered'] = plotting.water_delivered(df, layout,
                                                                 'Annual water delivered (Mm3)')

        df = crop_production.loc[crop_production['Governorate']==name]
        df = df.groupby(['Year', 'variable'])[['value']].sum() / 1000000
        df.reset_index(inplace=True)
        data[f'{name}CropProduction'] = plotting.plot_production(df, layout,
                                                                 'Annual cropland production by type (kton)',
                                                                 'variable')

    plots = []
    for key, value in data.items():
        plots.append(dcc.Graph(figure=value, config=dict(toImageButtonOptions=dict(
            format='png', filename=key, height=400,
            width=400, scale=2))))

    return plots, name


@app.callback(
    [Output('drop-level', 'options'), Output('drop-level', 'disabled'),
     Output('drop-level', 'value')],
    [Input('rb-scenario', 'value')]
)
def update_level_dropdown(scenario):
    level_dict = {'Reference': {'Select...': ''},
                  'Improve AG eff': {'Select...': '',
                                     'by 10 percent': ' by 10percent',
                                     'by 20 percent': ' by 20percent'},
                  'New Resources': {'Select...': ''},
                  'Reduce NRW': {'Select...': '',
                                 'to 30 percent': ' to 30 percent',
                                 'to 20 percent': ' to 20 percent'},
                  'Increased Water Productivity': {'Select...': ''}}
    options = [{"label": key, 'value': value} for key, value in level_dict[scenario].items()]
    disable = False
    if len(level_dict[scenario].keys()) == 1:
        disable = True
    return options, disable, ''


@app.callback(
    Output('scenarioTitle', 'children'),
    [
        Input("current", "modified_timestamp")
    ],
    [State('current', 'data')]
)
def update_level_dropdown(ts, data):
    if data is None:
        raise PreventUpdate
    if ts is None:
        raise PreventUpdate
    level_dict = {'Reference': {'': ''},
                  'Improve AG eff': {'': '',
                                     ' by 10percent': 'by 10 percent',
                                     ' by 20percent': 'by 20 percent'},
                  'New Resources': {'': ''},
                  'Reduce NRW': {'': '',
                                 ' to 30 percent': 'to 30 percent',
                                 ' to 20 percent': 'to 20 percent'},
                  'Increased Water Productivity': {'': ''}}
    scenario = data['scenario']
    level = data['level']
    name = f' - {scenario} {level_dict[scenario][level]} scenario'

    return name


@app.callback(
    [Output("map", "figure"), Output("loading-map", "children")],
    [
        Input('map-background', 'value'),
        Input('map-selection', 'value'),
    ],
)
def update_level_dropdown(background, map_type):
    map = plot_map(background, map_type)
    return {}, dcc.Graph(figure=map, id="map",
                         config=dict(showSendToCloud=True,
                                     toImageButtonOptions=dict(format='png',
                                                               filename='map',
                                                               height=700,
                                                               width=700, scale=2)))


# @app.callback(
# [Output('page-1', 'active'), Output('page-2', 'active'),
# Output('tools', 'hidden'), Output('visual-tools', 'hidden')],
# [Input('page-1', 'n_clicks_timestamp'), Input('page-2', 'n_clicks_timestamp')],
# )
# def selected_tab(n_1, n_2):
# if n_2 is None:
# n_1 = 1
# n_2 = 0
# tools = False
# visual = True
# elif n_1 is None:
# n_1 = 0
# n_2 = 1
# tools = True
# visual = False
# if n_1 > n_2:
# state_1 = True
# state_2 = False
# tools = False
# visual = True
# else:
# state_1 = False
# state_2 = True
# tools = True
# visual = False
# return state_1, state_2, tools, visual


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


# @app.callback(
#     [Output('cho-map-drop', 'disabled'),
#      Output('cho-map-drop', 'value'),
#      Output('cho-map-filter', 'disabled'),
#      Output('cho-map-filter', 'value'),
#      Output("cho-map-collapse", "is_open"),
#      ],
#     [Input('map-options', 'value')],
# )
# def enable_choromap_options(value):
#     if value == 'cho-map':
#         return False, None, False, None, True
#     return True, None, True, None, False

@app.callback(
    [
        Output('rb-scenario', 'value'),
        Output('eto-input', 'value'),
        Output('pump-eff-init', 'value'),
        Output('pump-eff-end', 'value'),
        # Output('map-options', 'value')
        # Output('compare-options', 'value'),
    ],
    [Input('button-reset', 'n_clicks')],
)
def reset_output(n):
    return 'Reference', ['Climate Change'], 0.45, 0.45


@app.callback(Output("download", "data"), [Input("button-download", "n_clicks")], [State('current', 'data')])
def func(n_nlicks, data):
    if data is None:
        raise PreventUpdate
    df = pd.DataFrame(data['gw_df'])
    return send_data_frame(df.to_csv, "gw.csv")


if __name__ == "__main__":
    app.run_server(debug=True)
