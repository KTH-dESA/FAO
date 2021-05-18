import os.path

import pandas as pd
import yaml
import boto3, gzip
from decouple import config

from_server = True
server = 'souss-massa-project'


AWS_ACCESS_ID = config('AWS_ACCESS_ID')
AWS_SECRET_KEY = config('AWS_SECRET_KEY')
AWS_REGION = config('AWS_REGION')

resource = boto3.resource(
    's3',
    aws_access_key_id=AWS_ACCESS_ID,
    aws_secret_access_key=AWS_SECRET_KEY,
    region_name=AWS_REGION
)


def get_path(path, from_server):
    if from_server:
        return ('/').join(path)
    else:
        return os.path.join(*path)


def load_summary_data(path, name, from_server=from_server):
    if from_server:
        path = server
    return pd.read_csv(get_path([path, 'data', name], from_server))


# def load_data(path, scenario, climate, phaseout_year, pv_level,
#               files='all', from_server=from_server):
#     if from_server:
#         path = server
#     init_year = 2020
#     end_year = 2050
#     butane_scenario = f'{phaseout_year}' if phaseout_year != 2050 else 'None'
#     if not climate:
#         climate = ['Trend']
#     data = get_path([path, 'data', scenario, climate[0]], from_server)
#     # lcoe = os.path.join(data_folder, scenario, climate[0], level)
#
#     if files == 'all':
#         files = ['results.gz', 'wwtp_data.gz', 'desal_data.gz',
#                  'butane.gz', 'production_data.gz']
#
#     if isinstance(files, str):
#         files = [files]
#
#     if len(files) == 1:
#         if files[0] == 'butane.gz':
#             dff = pd.read_csv(get_path([data,
#                                         'Butane Calculations',
#                                         butane_scenario,
#                                         f'{pv_level}',
#                                         files[0]], from_server))
#         else:
#             dff = pd.read_csv(get_path([data, files[0]], from_server))
#         # dff = dff.loc[(dff.Year >= init_year) & (dff.Year <= end_year)]
#         output = dff
#     else:
#         output = []
#         for file in files:
#             if file == 'butane.gz':
#                 dff = pd.read_csv(get_path([data,
#                                             'Butane Calculations',
#                                             butane_scenario,
#                                             f'{pv_level}',
#                                             file], from_server))
#             else:
#                 dff = pd.read_csv(get_path([data, file], from_server))
#             # dff = dff.loc[(dff.Year >= init_year) & (dff.Year <= end_year)]
#             output.append(dff)
#     return output
    
def load_data(path, scenario, climate, phaseout_year, pv_level,
              files='all', from_server=from_server):
    if from_server:
        path = server

    butane_scenario = f'{phaseout_year}' if phaseout_year != 2050 else 'None'
    if not climate:
        climate = ['Trend']
    data = get_path(['data', scenario, climate[0]], from_server)

    if files == 'all':
        files = ['results.gz', 'wwtp_data.gz', 'desal_data.gz',
                 'butane.gz', 'production_data.gz']

    if isinstance(files, str):
        files = [files]

    if len(files) == 1:
        if files[0] == 'butane.gz':
            obj = resource.Object(server, get_path([data,
                                        'Butane Calculations',
                                        butane_scenario,
                                        f'{pv_level}',
                                        files[0]], from_server))

            with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipfile:
                dff = pd.read_csv(gzipfile)
        else:
            obj = resource.Object(server, get_path([data, files[0]], from_server))

            with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipfile:
                dff = pd.read_csv(gzipfile)
        #dff = dff.loc[(dff.Year >= init_year) & (dff.Year <= end_year)]
        output = dff
    else:
        output = []
        for file in files:
            if file == 'butane.gz':
                obj = resource.Object(server, get_path([data,
                                                        'Butane Calculations',
                                                        butane_scenario,
                                                        f'{pv_level}',
                                                        file], from_server))

                with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipfile:
                    dff = pd.read_csv(gzipfile)
            else:
                obj = resource.Object(server, get_path([data, file], from_server))

                with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipfile:
                    dff = pd.read_csv(gzipfile)
            #dff = dff.loc[(dff.Year >= init_year) & (dff.Year <= end_year)]
            output.append(dff)
    return output


def get_language(language):
    file = f"assets/{language}.yaml"
    with open(file, 'rt', encoding='utf8') as yml:
        language_dic = yaml.load(yml, Loader=yaml.FullLoader)
    return language_dic