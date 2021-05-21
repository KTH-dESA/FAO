import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
import nexustool.weap_tools as wp
import re
from functools import reduce
import os
from shutil import copyfile


def extract_points_data(results, spacial_data, points, melt_on, regex):
    df = pd.DataFrame()
    for key, row in spacial_data.groupby(points):
        df_temp = results.filter(regex=regex.format(key)).copy()
        df_temp.columns = [i.split('\\')[-1] for i in df_temp.columns]
        df_temp[points] = key
        for i in melt_on:
            df_temp[i] = results[i]
        melt_vars = melt_on.copy()
        melt_vars.append(points)
        df = df.append(df_temp.melt(id_vars=melt_vars), sort=False)
    return df.reset_index(drop=True)
    
    
def extract_within(s, left, rigth):
    mo = re.search(f'{left}(.*){rigth}', s)
    if mo:
        return mo.group(1)
    return ''
    
    
def special_conditions(df, sheet_name):
    if sheet_name == 'Pipelines':
        df["PL_ZaraMain2SZ09SZ06 0 \ Reach"] = df["PL_ZaraMain2SZ09SZ06 0 \ Headflow"]
    
    
def get_data(sheet_names, data, spatial_data, variable, regex, 
             rename={'Unnamed: 0': 'Date'}, melt_on=['Date', 'Year', 'Month'], 
             look_in_sector=False):
    merged_data = pd.DataFrame()
    for sheet_name, value in sheet_names.items():
        data_temp = data.parse(sheet_name, skiprows=3)
        data_temp = data_temp.loc[data_temp['Unnamed: 0']!='Sum']
        data_temp.rename(columns=rename, inplace=True)
        if 'Date' in data_temp.columns:
            data_temp['Date'] = pd.to_datetime(data_temp['Date'])
            data_temp['Year'] = data_temp['Date'].dt.year
            data_temp['Month'] = data_temp['Date'].dt.month
        data_temp.columns = data_temp.columns.str.replace('"', '')
        data_temp.columns = data_temp.columns.str.replace('  ', ' ')
        data_temp.columns = data_temp.columns.str.strip()
        special_conditions(data_temp, sheet_name)
        if look_in_sector:
            merged_data_temp = extract_points_data(data_temp, 
                                                   spatial_data.loc[spatial_data['type']==value], 
                                                   variable, melt_on, regex)
        else:
            merged_data_temp = extract_points_data(data_temp, 
                                                   spatial_data, 
                                                   variable, melt_on, regex)
        merged_data_temp['variable'] = [i.split('[')[0].strip() for i in merged_data_temp.variable]
        merged_data_temp['type'] = value
        merged_data_temp.loc[merged_data_temp['type']=='variable', 'type'] = merged_data_temp.loc[merged_data_temp['type']=='variable', 'variable']

        merged_data = merged_data.append(merged_data_temp, sort=False)
    return merged_data
    
    
def enumerate_segments(df):
    x = len(df.Year.unique())*12
    df['n'] = 0
    for key, group in df.groupby('pipeline'):
        n = [[i]*x for i in range(int(group.shape[0]/x))]
        n = reduce(lambda a, b: a + b, n)
        df.loc[df.pipeline==key,'n'] = n
        

def get_elevation_delta(df):  
    df['elevation_delta'] = 0
    for pipeline in df.pipeline.unique():
        for n in range(df.loc[df.pipeline==pipeline, 'n'].max()):
            elevation_1 = df.loc[(df.pipeline==pipeline) & (df.n==n), 'elevation']
            elevation_2 = df.loc[(df.pipeline==pipeline) & (df.n==(n+1)), 'elevation']
            df.loc[(df.pipeline==pipeline) & (df.n==n), 'elevation_delta'] = np.array(elevation_2) - np.array(elevation_1)


def save_dataframe(df, init_year, end_year, output_folder, file_name):
    df.loc[(df.Year>=init_year) & (df.Year<=end_year)].to_csv(os.path.join(output_folder, file_name), index=False)