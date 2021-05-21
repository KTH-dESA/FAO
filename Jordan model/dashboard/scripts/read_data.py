import os

import pandas as pd


def load_data(path, scenario, eto, level, files='all'):
    data_folder = os.path.join(path, 'data')
    if not eto:
        eto = ['Historical Trend']
    data = os.path.join(data_folder, f'{scenario}{level}', eto[0])

    if files == 'all':
        files = ['water_delivered.gz', 'water_requirements.gz',
                 'groundwater_pumping.gz', 'pipelines_data.gz',
                 'wwtp_data.gz', 'desal_data.gz', 'crop_production.gz']
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
                                              'geometry': 'first'}).reset_index()

    df2 = supply_points.groupby('point').agg({'type': 'first',
                                              'geometry': 'first'}).reset_index()

    df_pipelines = pipelines.groupby('index').agg({'geometry': 'first'}).reset_index()

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


def get_demand_data(water_delivered, water_required, name):
    dff_delivered = water_delivered.loc[water_delivered['point'] == name]
    dff_delivered = dff_delivered.groupby(['Year', 'type', 'point'])['sswd'].sum() / 1000000
    dff_delivered = dff_delivered.reset_index()
    dff_delivered['label'] = 'Water delivered'
    dff_required = water_required.loc[water_required['point'] == name]
    dff_required = dff_required.groupby(['Year', 'type', 'point'])['sswd'].sum() / 1000000
    dff_required = dff_required.reset_index()
    dff_required['label'] = 'Water required'
    dff = dff_delivered.append(dff_required, ignore_index=True)

    unmet = round((dff_required['sswd'] - dff_delivered['sswd']) / dff_required['sswd'], 4)
    dff_unmet = pd.DataFrame({'Year': dff_required['Year'],
                              'Unmet demand': unmet})
    dff_unmet['Unmet demand'].fillna(1, inplace=True)
    return dff, dff_unmet