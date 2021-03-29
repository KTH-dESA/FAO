import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import numpy as np

from app import app
from scripts.read_data import get_language

layout = dbc.Nav([
    dbc.DropdownMenu(
        [dbc.DropdownMenuItem(["English ", html.I(className='fa fa-language')], className="drop-items", id="english"),
         dbc.DropdownMenuItem(["Spanish ", html.I(className='fa fa-language')], className="drop-items", id="spanish"),
         dbc.DropdownMenuItem(["French ", html.I(className='fa fa-language')], className="drop-items", id="french")],
        label="Language", id='language', nav=True, className="ml-2", disabled=True,
    ),
    dbc.NavItem(dbc.NavLink("About", id='about'), className="ml-2"),
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
])


@app.callback(
    Output("about-modal", "is_open"),
    [Input("about", "n_clicks"),
     Input("about-close", "n_clicks")],
    [State("about-modal", "is_open")],
)
def toggle_popover(n_1, n_2, is_open):
    if n_1 or n_2:
        return not is_open
    return is_open


@app.callback(
    Output('current-language', 'data'),
    [Input(language, 'n_clicks_timestamp') for language in ['english', 'spanish', 'french']],
    [State('current-language', 'data')],
    prevent_initial_call=True
)
def current_language(n1, n2, n3, language):
    language_list = ['english', 'spanish', 'french']
    n_list = []
    for n in [n1, n2, n3]:
        if n is None:
            n_list.append(0)
        else:
            n_list.append(n)
    if (n1 == n2) and (n1 == n3) and language:
        language = language
    else:
        language_index = np.array(n_list).argmax()
        language = language_list[language_index]
    return language


@app.callback(
    [Output('about', 'children'),
     Output('about-title', 'children'),
     Output('about-body', 'children'),
     Output('about-close', 'children'),
     Output('language', 'label')],
    [Input('current-language', 'modified_timestamp')],
    [State('current-language', 'data')],
)
def update_language(ts, language):
    if not language:
        raise PreventUpdate
    language_dic = get_language(language)

    about_content = [html.P(language_dic['about']['body'][0]), html.P(language_dic['about']['body'][1]),
                     html.P(language_dic['about']['body'][2]), html.P(language_dic['about']['body'][3]),
                     dbc.Row([dbc.Col(html.A(html.Img(src='../assets/kth.png', style={'height': '130px'}),
                                             href='https://www.energy.kth.se/energy-systems/about-the-division-of-energy-systems-1.937036'),
                                      width=3),
                              dbc.Col(html.A(html.Img(src='../assets/sei.png', style={'height': '130px'}),
                                             href='https://www.sei.org/'), width=4),
                              dbc.Col(html.A(html.Img(src='../assets/fao.png', style={'height': '130px'}),
                                             href='http://www.fao.org/home/en/'), width=2)], justify="center")
                     ]

    return language_dic['about']['header'], \
           language_dic['about']['title'], \
           about_content, language_dic['about']['close'], \
           language_dic['language']
