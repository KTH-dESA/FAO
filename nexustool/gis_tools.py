import os
import pandas as pd
import numpy as np

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


