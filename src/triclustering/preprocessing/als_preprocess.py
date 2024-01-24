import pandas as pd
import datetime as dt
import math
import numpy as np
import constants


def df_to_dict(data, discretize_prog_rate=False):
    """
    Transforms a DataFrame into a Dict
    """

    data_dict = data.to_dict('index')
    final_dict = dict()

    id_paciente_glob = 0
    time_counter = 0
    for k in data_dict.keys():
        ref = data_dict[k][constants.REF_FEATURE]  # ID Paciente da Iteração
        if ref != id_paciente_glob:
            id_paciente_glob = ref
            time_counter = 0

        del data_dict[k][constants.REF_FEATURE]
        if ref not in final_dict:
            final_dict[ref] = {time_counter: data_dict[k]}
        else:
            final_dict[ref][time_counter] = data_dict[k]

        time_counter += 1

    return final_dict

def compute_consecutive_snapshots_n(data, n, label, yes_label='Y'):
    """

    Parameters
    ----------
    data: is a dict with ALS data with the format returned by `df_to_dict`
    n: is the number of consecutive snapshots to consider, ie. the size of snapshots set
        the size of snapshots set could be defined 
    label: is the target problem
    strategy: (default) `flexible` - sets of snapshots have a maximum size `n`
                `strict` - sets of snapshots have a strict size of `n`

    """

    final = dict()
    for (p, t) in data.items():
        if len(t.keys()) >= n:
            fd = dict()
            for (key, val) in t.items():
                fd[key] = val
                final[p] = fd

    snaps = dict()
    for (p, ts) in data.items():
        for t in ts.keys():

            size_t = len(ts.keys())
            #size_n = min(n, size_t)
            size_n = n
            if t < size_t - (size_n-1) and all(map(lambda c: c != yes_label, [data[p][t+y][label] for y in range(0, size_n-1)])):
                if p not in snaps:
                    snaps[p] = list()
                snaps[p].append([(t+j, data[p][t+j][label])
                                for j in range(0, size_n)])
    return snaps

def create_matrix_temporal(data, sps, n):
    y = list()
    values = list()
    cols = list()
    cols.append("Patient_ID")
    for p in sps.keys():
        tp = data[p]
        for snaps in sps[p]:
            l = list()
            l.append(p)

            for e in snaps:
                i = e[0]
                l.extend([tp[i][feature]
                         for feature in tp[i].keys() if feature != "Evolution"])

            values.append(l)
            y.append(e[1])

    cols.extend([f"{ti}{feature}" for ti in range(n)
                 for feature in tp[i].keys() if feature != "Evolution"])

    mats = pd.DataFrame(data=values,
                        columns=cols)

    return mats, y

