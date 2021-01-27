import os.path

import pandas as pd

from_server = False
server = 's3://souss-massa-dev'

def get_path(path, from_server):
    if from_server:
        return ('/').join(path)
    else:
        return os.path.join(*path)


def load_summary_data(path, name, from_server=from_server):
    if from_server:
        path = server
    return pd.read_csv(get_path([path, 'data', name], from_server))


def load_data(path, scenario, climate, phaseout_year, pv_level,
              files='all', from_server=from_server):
    if from_server:
        path = server
    init_year = 2020
    butane_scenario = f'phaseout_{phaseout_year}' if phaseout_year != 2050 else 'phaseout_None'
    if not climate:
        climate = ['Trend']
    data = get_path([path, 'data', scenario, climate[0]], from_server)
    # lcoe = os.path.join(data_folder, scenario, climate[0], level)

    if files == 'all':
        files = ['results.gz', 'wwtp_data.gz', 'desal_data.gz',
                 'summary_results.gz', 'production.gz']

    if isinstance(files, str):
        files = [files]

    if len(files) == 1:
        dff = pd.read_csv(get_path([data, files[0]], from_server))
        dff = dff.loc[dff.Year >= init_year]
        output = dff
    else:
        output = []
        for file in files:
            if file == 'summary_results.gz':
                dff = pd.read_csv(get_path([path,
                                            'data',
                                            'Butane_calculations',
                                            butane_scenario,
                                            f'{pv_level}_PV',
                                            file], from_server))
            else:
                dff = pd.read_csv(get_path([data, file], from_server))
            dff = dff.loc[dff.Year >= init_year]
            output.append(dff)
    return output