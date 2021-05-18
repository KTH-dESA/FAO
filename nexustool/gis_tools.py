import itertools
import os
import urllib.request
import zipfile
import pandas as pd
import geopandas as gpd
import numpy as np
from rasterstats import zonal_stats


def homogenize_provinces(df, case):
    if case == 'morroco':
        df.loc[df['Province'] == 'Chtouka-Ait Baha', 'Province'] = 'Chtouka-Aït Baha'
        df.loc[df['Province'] == 'Inezgane-Ait Melloul', 'Province'] = 'Inezgane-Aït Melloul'
        df.loc[df['Province'] == 'Agadir-Ida-Ou-Tanane', 'Province'] = 'Agadir-Ida ou Tanane'


def disaggregate_data(data, spatial_layer, case='morocco'):
    homogenize_provinces(data, case)
    df = spatial_layer.loc[spatial_layer.Date.isin(data.Date.unique())].set_index(
        ['province', 'Date'])  # this should be the results dataframe

    summary_provinces_agri = data.loc[data['type'].str.contains('Agriculture')].groupby(['Province', 'Date'])[
        ['sswd', 'pwd', 'pa_e', 'pp_e']].sum()
    temp_cropland_provinces = df.reset_index()[['province', 'Date']].copy()
    temp_cropland_provinces.loc[
        temp_cropland_provinces['province'] == 'Inezgane-Aït Melloul', 'province'] = 'Taroudannt'
    for feature in list(summary_provinces_agri):
        df[feature] = temp_cropland_provinces.set_index(['province', 'Date']).index.map(
            summary_provinces_agri[feature]) * df['area_share']
    df.reset_index(inplace=True)

    return df


def get_zonal_stats(vector, path, stats, all_touched=False):
    # Run zonal statistics, store result in geopandas dataframe
    # with rasterio.open(path) as src:
    result = zonal_stats(vector, path, stats=stats,
                         geojson_out=True, all_touched=all_touched)
    geostats = gpd.GeoDataFrame.from_features(result)
    return geostats


def create_time_data(data, first_year, end_year):
    df = pd.DataFrame()
    df['Year'] = list(itertools.chain.from_iterable([[year]*12*data.shape[0] for year in range(first_year,end_year+1)]))
    y_number = df.Year.unique().shape[0]
    dff = pd.DataFrame(data.filter(regex='srad|province|Demand point|area|wtd').melt(id_vars=['Demand point','province','area_m2','wtd']))
    dff.rename(columns={'variable': 'Month', 'value': 'srad'}, inplace=True)
    # dff = dff.join(data.filter(like='_wind').melt())
    # dff.rename(columns={'value': 'wind'}, inplace=True)
    # dff.drop(columns=['variable'], inplace=True)
    dff['Month'] = dff['Month'].str.replace('srad', '').astype(int)
    dff.sort_values(['province','Month'], inplace=True)
    dff.reset_index(drop=True, inplace=True)
    df = df.join(pd.concat([dff]*y_number, ignore_index=True))
    df['Date'] = pd.to_datetime(df[['Year','Month']].join(pd.DataFrame({'day': [1]*df.shape[0]})))
    return df


def get_area_share(df, by, value):
    dff = df.groupby(by)[value].sum().copy()
    return df[value] / df.set_index(by).index.map(dff)


def download_data(url, file_path):
    if not os.path.isfile(file_path):
        print('Downloading...')
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        urllib.request.urlretrieve(url, file_path)

        if os.path.splitext(file_path)[1] == '.zip':
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(os.path.dirname(file_path))
        print('Done.')
    else:
        print('The data has been already downloaded.')

