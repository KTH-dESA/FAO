import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

layout = html.Div([
    dcc.Link('Navigate to Explorer', href='/explorer'),
])
