import preprocessing.als_preprocess as als

import pandas as pd
import constants
import sys
from pathlib import Path


def load_data_baselines(features):
    infile = constants.DATA_FILE

    data = pd.read_csv(infile, low_memory=False)

    data['Evolution'] = 'No'
    data = data[features]

    print("y INI : ", len(data))

    return data


def compute_temporal_table(data, n, features):

    data = data[features]
    data_dict = als.df_to_dict(data)

    sps = als.compute_consecutive_snapshots_n(
        data_dict, n, 'Evolution')

    mats, y = als.create_matrix_temporal(data_dict, sps, n)

    mats.fillna(0, inplace=True)

    baseline_temporal = constants.BASELINE_DIR_T +"{}TPS/".format(n)

    Path(baseline_temporal).mkdir(parents=True, exist_ok=True)
    mats['Evolution'] = 'No'
    mats = mats.groupby('Patient_ID').first().reset_index()
    mats.to_csv(baseline_temporal +
                "{}TPS_baseline_temporal.csv".format(n), index=False)


n = int(sys.argv[1])
constants.get_config(sys.argv[2])
features = [constants.REF_FEATURE] + list(constants.TEMPORAL_FEATURES.keys()) + ['Evolution']

data = load_data_baselines(features)

t = compute_temporal_table(
    data, n, [constants.REF_FEATURE] + list(constants.TEMPORAL_FEATURES.keys()) + ['Evolution'])

