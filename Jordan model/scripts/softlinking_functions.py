import sys
sys.path.append("..") #this is to add the avobe folder to the package directory
import geopandas as gpd
import pandas as pd
import numpy as np
import nexus_tool.weap_tools as wp
import re
from functools import reduce
import os
from shutil import copyfile


def extract_weap(results, spacial_data, points, melt_on, regex):
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
    
    
def extract(s, left, rigth):
    mo = re.search(f'{left}(.*){rigth}', s)
    if mo:
        return mo.group(1)
    return ''
    
    
def get_demand_data(sheet_names, data, spatial_data, variable, regex):
    merged_data = pd.DataFrame()
    for sheet_name, value in sheet_names.items():
        data_temp = data.parse(sheet_name, skiprows=5)
        data_temp.rename(columns={'$Columns = Year': 'Year', ' Timestep': 'Month'}, inplace=True)
        data_temp.columns = data_temp.columns.str.replace('"', '')
        merged_data_temp = extract_weap(data_temp, spatial_data, variable, ['Year','Month'], regex)
        merged_data_temp['variable'] = [i.split('[')[0] for i in merged_data_temp.variable]
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
