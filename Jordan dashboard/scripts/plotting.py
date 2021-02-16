import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import os
import numpy as np


def load_data(path, scenario, eto, level, files='all'):
    data_folder = os.path.join(path, 'data')
    if not eto:
        eto = ['Without Eto trend']
    data = os.path.join(data_folder, scenario, eto[0], level)

    if files == 'all':
        files = ['Water_delivered.csv', 'Water_requirements.csv',
                 'Groundwater_pumping.csv', 'Pipelines_data.csv',
                 'wwtp_data.csv', 'desal_data.csv', 'crop_production.csv']
    if isinstance(files, str):
        files = [files]

    if len(files) == 1:
        output = pd.read_csv(os.path.join(data, files[0]))
    else:
        output = []
        for file in files:
            output.append(pd.read_csv(os.path.join(data, file)))
    return output


def data_merging(demand_points, supply_points, pipelines):
    df1 = demand_points.groupby('point').agg({'type': 'first',
                                              'elevation_m': 'first',
                                              'geometry': 'first'}).reset_index()

    df2 = supply_points.groupby('point').agg({'type': 'first',
                                              'elevation_m': 'first',
                                              'geometry': 'first'}).reset_index()

    df_pipelines = pipelines.groupby('index').agg({'pipeline': 'first',
                                                   'segment_length_m': 'first',
                                                   'geometry': 'first'}).reset_index()

    df = df1.append(df2, ignore_index=True)
    df['lon'] = [point.xy[0][0] for point in df.geometry]
    df['lat'] = [point.xy[1][0] for point in df.geometry]

    pipe_coords = pd.DataFrame({'lon': [], 'lat': []})
    for name, point in zip(df_pipelines.pipeline, df_pipelines.geometry):
        lon = list(point.xy[0]) + [None]
        lat = list(point.xy[1]) + [None]
        df_temp = pd.DataFrame({'lon': lon, 'lat': lat})
        df_temp['name'] = name
        pipe_coords = pipe_coords.append(df_temp, ignore_index=True)

    pipe_coords['type'] = 'pipeline'
    return df, pipe_coords


def choroplet_map(geojson, df):
    fig = px.choropleth_mapbox(df, geojson=geojson, locations='id',
                               color='color',
                               color_continuous_scale=px.colors.sequential.Viridis,
                               custom_data=['id'])

    fig.update_layout(coloraxis_colorbar=dict(
        len=0.5,
        xanchor="right", x=1,
        yanchor='bottom', y=0.1,
        thickness=10,
    ))

    # fig = go.Figure(go.Choroplethmapbox(geojson=geojson, locations=df['id'],
    #                            z=df['color'],
    #                            colorscale='Viridis'))
    return fig


def plot_borders(geojson):
    opacity = 0.1
    line = dict(width=1, color='gray')
    colorscale = ((0, 'gray'), (1, 'gray'))
    df = pd.DataFrame(geojson['features'])
    df['color'] = 0

    trace = go.Choroplethmapbox(geojson=geojson, locations=df.id, z=df['color'],
                                colorscale=colorscale,
                                marker=dict(opacity=opacity,
                                            line=line),
                                hoverinfo='skip', showscale=False,
                                showlegend=False)
    return trace


def plot_points(df):
    trace = px.scatter_mapbox(df, lat="lat", lon="lon", color="type",
                              color_discrete_sequence=px.colors.qualitative.Dark2,
                              custom_data=['type', 'point'], hover_name='point').data
    return trace


def plot_pipelines(df):
    color = 'rgb(100,100,100)'
    trace = px.line_mapbox(df, lat="lat", lon="lon", color='type',
                           color_discrete_map={'pipeline': color},
                           custom_data=['type', 'name'], hover_name='name')
    trace.update_traces(mode='markers+lines', marker={'opacity': 0})
    return trace.data


def water_delivered(df, layout, title):
    fig = px.bar(df, x='Year', y='value', color='type',
                  color_discrete_sequence=px.colors.qualitative.Dark2)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, title=title)
    return fig


def unmet_demand(df, layout, title):
    fig = px.line(df, x='Year', y='value', color='type',
                  color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:,.2%}' + '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, yaxis=dict(tickformat='%'), title=title)
    return fig


def wtd_plot(df, layout, title):
    fig = px.line(df, x='Year', y='wtd', color='point',
                  color_discrete_sequence=px.colors.qualitative.Vivid)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' + '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, yaxis=dict(range=[520, 0]), title=title)
    return fig


def groundwater_extraction(df, layout, title):
    fig = px.bar(df, x='Year', y='value', color='point',
                  color_discrete_sequence=px.colors.qualitative.G10)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, title=title)
    return fig


def energy_demand(df, layout, title):
    # eff_val = 0.65
    # if eff_val > 0.45:
    #     eff = np.arange(0.45, eff_val+0.00001, (eff_val-0.45) / (2050-2020))
    #     df = df.copy()
    #     df = pd.merge(df, pd.DataFrame({'Year': [i for i in range(2020, 2051)],
    #                                     'Eff': eff}), on='Year')
    #     df['valueWithEff'] = df['value'] / df['Eff']
    # else:
    #     df['valueWithEff'] = df['value'] / 0.45
    fig = px.bar(df, x='Year', y='value', color='type',
                  color_discrete_sequence=px.colors.qualitative.T10)
    fig.update_traces(hovertemplate='<b>Value</b>: %{y:.2f}' +
                                    '<br><b>Year</b>: %{x}')
    fig.update_layout(layout, title=title)
    return fig


def get_demand_data(water_delivered, water_required, name):
    dff_delivered = water_delivered.loc[water_delivered['point'] == name]
    dff_delivered = dff_delivered.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_delivered['label'] = 'Water delivered'
    dff_required = water_required.loc[water_required['point'] == name]
    dff_required = dff_required.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_required = dff_required.reset_index()
    dff_required['label'] = 'Water required'
    dff = dff_delivered.append(dff_required, ignore_index=True)

    unmet = round((dff_required['value'] - dff_delivered['value']) / dff_required['value'], 4)
    dff_unmet = pd.DataFrame({'Year': dff_required['Year'],
                              'Unmet demand': unmet})
    dff_unmet['Unmet demand'].fillna(1, inplace=True)
    return dff, dff_unmet


def plot_water_delivered(df, layout, title):
    """
    :param df: tidy dataframe containing year and water delivered value
    :param layout: dictionary containing layout options for Plotly
    :param title: title of the plot
    :return: Plotly figure object
    """
    df['value'] = round(df['value'], 4)
    fig = px.line(df, x='Year', y='value', color='label',
                  color_discrete_sequence=px.colors.qualitative.Dark2)
    fig.update_traces(fill='tonexty')
    fig.update_layout(layout, title=title)
    return fig


def plot_unmet_demand(df, layout, title):
    """
    :param df: tidy dataframe containing year and unmet demand value (fraction)
    :param layout: dictionary containing layout options for Plotly
    :param title: title of the plot
    :return: Plotly figure object
    """
    fig = px.line(df, x='Year', y='Unmet demand',
                  hover_data={'Unmet demand': ':.2%'})
    fig.update_layout(layout, title=title, yaxis={'tickformat': '%'})
    return fig


def plot_water_supply(df, colors, layout, title):
    """
    :param df: tidy dataframe containing year and value of water supplied
    :param colors: list of colors to use
    :param layout: dictionary containing layout options for Plotly
    :param title: title of the plot
    :return: Plotly figure object
    """
    df['value'] = round(df['value'], 4)
    fig = px.area(df, x='Year', y='value', color_discrete_sequence=colors)
    fig.update_layout(layout, title=title)
    return fig


def plot_depth_groundwater(df, layout, title):
    """
    :param df: tidy dataframe containing year and water table dpth (wtd) value
    of the aquifer
    :param layout: dictionary containing layout options for Plotly
    :param title: title of the plot
    :return: Plotly figure object
    """
    fig = px.line(df, x='Date', y='wtd')
    fig.update_layout(layout, title=title, yaxis={'range': [df.wtd.max() * 1.1,
                                                            df.wtd.min() * 0.9]})
    return fig


def plot_energy_for_pumping(df, colors, layout, title):
    """
    :param df: tidy dataframe containing year and energy demand (SWPA_E_) for
    pumping
    :param colors: list of colors to use
    :param layout: dictionary containing layout options for Plotly
    :param title: title of the plot
    :return: Plotly figure object
    """
    df['value'] = round(df['SWPA_E_'], 4)
    fig = px.area(df, x='Year', y='SWPA_E_',
                  color_discrete_sequence=colors)
    fig.update_layout(layout, title=title)
    return fig


def plot_production(df, layout, title, color):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']
    fig = px.bar(df, x='Year', y='production',
                 color=color,
                 title=title,
                 color_discrete_sequence=colors)
    fig.update_layout(layout, height=450,
                      # bargap=0
                      )
    return fig

def plot_production_by_gov(df, layout, title, y, color):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']
    fig = px.bar(df, x='production', y=y, color=color,
                 title=title, orientation="h",
                 color_discrete_sequence=colors)
    fig.update_layout(layout, height=500, yaxis={'categoryorder': 'total ascending'},
                      legend=dict(
                          orientation="h",
                          # yanchor="top",
                          # y=0,
                          # xanchor="left",
                          x=-0.2
                      ))
    return fig


def plot_water_delivered_by_gov(df, layout, title):
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']
    fig = px.area(df, x='Year', y='value', color='Governorate',
                  title=title,
                  color_discrete_sequence=colors)
    fig.update_layout(layout)
    return fig
