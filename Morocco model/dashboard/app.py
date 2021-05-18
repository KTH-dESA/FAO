import os

import dash
import dash_bootstrap_components as dbc

my_path = os.path.join(os.path.abspath(os.path.dirname(__file__)))

app = dash.Dash(__name__, suppress_callback_exceptions=True,
                external_stylesheets=[dbc.themes.BOOTSTRAP,
                                      'https://use.fontawesome.com/releases/v5.12.0/css/all.css'],
                # these meta_tags ensure content is scaled correctly on different devices
                # see: https://www.w3schools.com/css/css_rwd_viewport.asp for more
                meta_tags=[
                    {"name": "viewport", "content": "width=device-width, initial-scale=1"}
                ],
                )
server = app.server
server.secret_key = os.environ.get('secret_key', 'secret')
