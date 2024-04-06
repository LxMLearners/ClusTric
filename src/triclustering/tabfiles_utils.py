import pandas as pd
import constants
import preprocessing.als_dataframe as alsd
from preprocessing.als_preprocess import label_encoder_als

def tab_from_baseline_temporal(baseline_file, outfile, n, features):
    data = pd.read_csv(baseline_file)

    try:
        data.drop(columns=['Evolution'], inplace=True)
    except:
        pass
    data.fillna(0, inplace=True)

    # features = ['ALSFRSb', 'ALSFRSsUL',
    #             'ALSFRSsLL', 'R', 'ALSFRS-R', '%FVC', 'MITOS-stage']
    alsd.write_tab_file_temp(data, outfile, features, n)

    return len(data)

def tab_from_baseline_static(baseline_file, outfile, features):
    data = pd.read_csv(baseline_file)
    data.fillna('0', inplace=True)
    data = label_encoder_als(data, features)
    data = data.round(1)
    # features = ['Gender', 'BMI', 'MND familiar history', 'Age at onset', 'Disease duration',
    #           'El Escorial reviewed criteria', 'UMN vs LMN', 'Onset form', 'C9orf72']
    alsd.write_tab_file(data, outfile, features, 1)

