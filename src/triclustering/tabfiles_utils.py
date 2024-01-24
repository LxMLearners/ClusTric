import pandas as pd
import constants
import preprocessing.als_dataframe as alsd

def tab_from_baseline_temporal(baseline_file, outfile, n, features):
    data = pd.read_csv(baseline_file)

    try:
        data.drop(columns=[constants.TARGET], inplace=True)
    except:
        pass
    data.fillna(0, inplace=True)

    # features = ['ALSFRSb', 'ALSFRSsUL',
    #             'ALSFRSsLL', 'R', 'ALSFRS-R', '%FVC', 'MITOS-stage']
    alsd.write_tab_file(data, outfile, features, n)

    return len(data)
