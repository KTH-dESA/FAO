import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import itertools

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

def plot_borders(geojson, df=None):
    if df is None:
        df = pd.DataFrame(geojson['features'])
        df['color'] = 0
        colorscale = ((0, 'gray'), (1, 'gray'))

    trace = go.Choroplethmapbox(geojson=geojson, locations=df.id, z=df['color'],
                                colorscale=colorscale, marker=dict(opacity=0.1),
                                hoverinfo='skip', showscale=False,
                                showlegend=False, marker_line_width=0.5)
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
                           hover_name='name').data
    return trace
    
def get_municipal_data(water_delivered, water_required):
    dff_delivered = water_delivered.loc[water_delivered['point'] == name]
    dff_delivered = dff_delivered.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_delivered['label'] = 'Water delivered'
    dff_required = water_required.loc[water_required['point'] == name]
    dff_required = dff_required.groupby(['Year', 'type', 'point'])['value'].sum() / 1000000
    dff_required = dff_required.reset_index()
    dff_required['label'] = 'Water required'

    dff = dff_delivered.append(dff_required, ignore_index=True)
    dff['Unmet demand'] = round((dff_required['value'] - dff_delivered['value']) / dff_required['value'], 4)
    return dff