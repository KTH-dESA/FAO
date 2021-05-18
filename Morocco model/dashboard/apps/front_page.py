# -*- coding: utf-8 -*-

import os
import pickle

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.io as pio

from app import app, my_path
from apps import about
from apps.explorer import charts_layout
from scripts import plotting
from scripts.read_data import get_language, load_summary_data

pio.templates.default = "plotly_white"

content_folder = os.path.join(my_path, 'assets', 'report_content')


def results_element(body, chart):
    element = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(dcc.Markdown(body), lg=5, md=12, style=dict(marginTop=30)),
                    dbc.Col(dcc.Graph(figure=chart), lg=7, md=12),
                ],
            ),
        ]
    )
    return element


content_menu = dbc.Nav(
    [
        dbc.Col(dbc.DropdownMenu(
            [html.A(dbc.DropdownMenuItem(className="drop-items", id='first-section'), href='#first'),
             html.A(dbc.DropdownMenuItem(className="drop-items", id='second-section'), href='#second'),
             html.A(dbc.DropdownMenuItem(className="drop-items", id='third-section'), href='#third')],
            label="Content", id='content-menu', nav=True, className="content"),
        ),
    ],
    # align="center",
    # no_gutters=True,
)

options_nav = dbc.Row(
    [
        about.layout,
        dbc.Col(
            dbc.Button("Explore the model", id='explore-model', color="primary",
                       className="ml-2", href='/explorer'),
            width="auto",
        ),
    ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)

navbar = dbc.Navbar(
    children=[
        content_menu,
        dbc.NavbarBrand(id='model', className="ml-2"),
        options_nav
    ],
    id='navbar',
    color="light",
    dark=False,
    sticky='top'
)

content_layout = dbc.Row(
    [
        dbc.Col(width=11,
                id='content-front-page')
    ],
    justify="center",
    no_gutters=True
)

layout = html.Div([
    navbar,
    content_layout
], dir='ltr')


@app.callback(
    [Output('explore-model', 'children'),
     Output('model', 'children'),
     Output('content-menu', 'label')],
    [Input('current-language', 'modified_timestamp')],
    [State('current-language', 'data')],
)
def update_language(ts, language):
    if not language:
        raise PreventUpdate
    language_dic = get_language(language)

    return [html.Span(className="fa fa-chart-bar", style=dict(marginRight='6px')),
            language_dic['front page']['explore']], \
           f"Souss-Massa {language_dic['sidebar']['title']}", \
           language_dic['front page']['content menu']


@app.callback(
    [Output('first-section', 'children'),
     Output('second-section', 'children'),
     Output('third-section', 'children')],
    [Input('current-language', 'modified_timestamp')],
    [State('current-language', 'data')],
)
def update_language(ts, language):
    if not language:
        raise PreventUpdate
    language_dic = get_language(language)

    return language_dic['front page']['sections']


@app.callback(
    Output("popover", "is_open"),
    [Input("popover-target", "n_clicks")],
    [State("popover", "is_open")],
)
def toggle_popover(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output('content-front-page', "children"),
    [Input('current-language', 'modified_timestamp')],
    [State('current-language', 'data')],
)
def toggle_popover(ts, language):
    if not language:
        raise PreventUpdate
    with open(os.path.join(content_folder, 'charts', 'my_file.pkl'), 'rb') as f:
        fig = pickle.load(f)
    language_dic = get_language(language)
    content = [html.H4(language_dic['front page']['sections'][0], id='first', className='header'),
               results_element(
                   open(os.path.join(content_folder, language, '1. first_section.txt'), encoding="utf-8").read(),
                   {}),
               dcc.Graph(figure=fig),
               html.Br(),
               html.Hr(),
               results_element(
                   open(os.path.join(content_folder, language, '2. first_section.txt'), encoding="utf-8").read(),
                   fig),
               html.Br(),
               html.Hr(className='section-hr', id='second'),
               html.H4(language_dic['front page']['sections'][1], className='header'),
               results_element(
                   open(os.path.join(content_folder, language, '2. first_section.txt'), encoding="utf-8").read(),
                   {}),
               html.Br(),
               html.Hr(className='section-hr', id='third'),
               html.H4(language_dic['front page']['sections'][2], className='header')
               ]
    return content
