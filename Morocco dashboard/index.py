# -*- coding: utf-8 -*-

import os

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.io as pio

from app import app
from apps import explorer, front_page

pio.templates.default = "plotly_white"

app.layout = html.Div([
    dcc.Store(id='current-language'),
    dcc.Location(id='url', refresh=False),
    html.Div(id='app-content'),
])

# app.validation_layout = html.Div([
#     explorer.layout,
#     front_page.layout,
#     app.layout
# ])


@app.callback(Output('app-content', 'children'),
              Input('url', 'pathname'))
def display_page(pathname):
    if pathname == '/explorer':
        return explorer.layout
    else:
        return front_page.layout


if __name__ == '__main__':
    app.run_server(debug=True)
