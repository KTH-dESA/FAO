import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import contextily as ctx
from matplotlib import cm
import os

cropland = gpd.read_file(str(snakemake.params.cropland), encoding='utf-8')
results = str(snakemake.params.results)
scenario = str(snakemake.params.scenario)
climate = str(snakemake.params.climate)
w = (snakemake.params.w_rate)
pv = str(snakemake.params.pv_rate)
grid = str(snakemake.params.grid_rate)

if climate == 'Trend':
    climate_title = f'Historical climate {climate}'
else:
    climate_title = 'Increased droughts'

folder = results.split(os.path.basename(results))[0]
df = pd.read_csv(results)
provinces = gpd.read_file('Data/GIS/Admin/Provinces.gpkg', encoding='utf-8')
provinces.to_crs('epsg:3857', inplace=True)

for year in range(2020, 2051):
    dff = df.loc[(df.Year==year)].copy()
    dff['geometry'] = dff['Demand point'].map(cropland.geometry)
    dff = gpd.GeoDataFrame(dff, crs='epsg:26192')
    dff.to_crs('epsg:3857', inplace=True)
    fig, ax = plt.subplots(1, 1, figsize=(10,7))

    provinces.plot(color=(0.4,0.5,0.4,0.2), edgecolor='black', ax=ax)
    dff.plot(column='least_cost_technology', categorical=True, ax=ax, legend=True)
    ctx.add_basemap(ax)
    ax.set_axis_off()
    fig.suptitle(f'Least-cost technology options {scenario} scenario year {year}',
                 fontsize=16)
    plt.figtext(.5,.9,f"{climate_title}, Wind {int(float(w)*100)}%, PV {int(float(pv)*100)}%, Grid {int(-float(grid)*100)}%", fontsize=14, ha='center')
    fig.savefig(f'{folder}/least-cost-{year}.png', dpi=150)
    plt.close('all')