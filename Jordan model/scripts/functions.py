import os

import pandas as pd
from shutil import copyfile


def save_results(jordan, jordan_gw, jordan_ww, jordan_desal, 
                 folder, template=None):
    os.makedirs(folder, exist_ok=True)
    
    if template:
        copyfile(os.path.join(template, 'crop_production.gz'), os.path.join(folder, 'crop_production.gz'))
        copyfile(os.path.join(template, 'water_delivered.gz'), os.path.join(folder, 'water_delivered.gz'))
        copyfile(os.path.join(template, 'water_requirements.gz'), os.path.join(folder, 'water_requirements.gz'))

    jordan.df.to_csv(os.path.join(folder, 'pipelines_data.gz'), index=False)
    jordan_gw.df.to_csv(os.path.join(folder, 'groundwater_pumping.gz'), index=False)
    jordan_ww.df.to_csv(os.path.join(folder, 'wwtp_data.gz'), index=False)
    jordan_desal.df.to_csv(os.path.join(folder, 'desal_data.gz'), index=False)


def merge_scenario_data(path, scenarios):
    all_results = [pd.DataFrame()] * 7
    for scenario in scenarios:
        results = load_scenarios(os.path.join(path, scenario))

        for dff in results: 
            dff['Scenario'] = scenario

        for i, dff in enumerate(results):
            all_results[i] = all_results[i].append(dff, ignore_index=True)

    df_delivered = all_results[0]
    df_required = all_results[1]
    df_gw = all_results[2]
    df_pipelines = all_results[3]
    df_wwtp = all_results[4]
    df_desal = all_results[5]
    df_crop = all_results[6]
    
    return df_delivered, df_required, df_gw, df_pipelines, df_wwtp, df_desal, df_crop
    
    
def load_scenarios(path):
    files = ['water_delivered.gz', 'water_requirements.gz',
             'groundwater_pumping.gz', 'pipelines_data.gz',
             'wwtp_data.gz', 'desal_data.gz', 'crop_production.gz']
    output = []
    for file in files:
        output.append(pd.read_csv(os.path.join(path, file)))
    
    return output