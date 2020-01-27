WebMercator = 4326
governorates.to_crs(epsg=WebMercator, inplace=True)

governorates.to_file("governorates.json", driver="GeoJSON")

token = "pk.eyJ1IjoiY2FtaWxvcmciLCJhIjoiY2p1bTl0MGpkMjgyYjQ0b2E0anRibWJ1MSJ9.GhUUGD6gok1d36lvP17CQQ"

import json
with open('governorates.json') as response:
    counties = json.load(response)
    for feature in counties['features']:
        feature['id'] = feature['properties']['id']

import pandas as pd
df = governorates

import plotly.graph_objects as go

fig = go.Figure(go.Choroplethmapbox(geojson=counties, locations=df.id, z=df.GLCV2sum,
                                    colorscale="Viridis", marker_line_width=1))
fig.update_layout(mapbox_style="light", mapbox_accesstoken=token,
                  mapbox_zoom=6, mapbox_center = {"lat": 31.4, "lon": 36.5})
fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
fig.show()

from urllib.request import urlopen
import json
with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
    counties = json.load(response)

import pandas as pd
df = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/fips-unemp-16.csv",
                   dtype={"fips": str})

import plotly.graph_objects as go

fig = go.Figure(go.Choroplethmapbox(geojson=counties, locations=df.fips, z=df.unemp,
                                    colorscale="Viridis", zmin=0, zmax=12, marker_line_width=0))
fig.update_layout(mapbox_style="light", mapbox_accesstoken=token,
                  mapbox_zoom=3, mapbox_center = {"lat": 37.0902, "lon": -95.7129})
fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
fig.show()