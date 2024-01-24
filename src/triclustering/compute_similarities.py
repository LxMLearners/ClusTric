import pandas as pd
import constants
import os
from tricluster import Tricluster
from pathlib import Path
from similarity_utils import get_triclusters, compute_similarity_matrix, write_matrix

def sim_matrix_tric(n, baseline_file, path_trics, path_matr, last=False, tri=False):
    """
    Computes similarity matrices between triclusters and original data

    Parameters
    ----------
    n: set of snapshots size
    baseline_file: 
    path_trics:
    path_matr:
    last: True if (TODO)
    tri: `True` if similarities are computed using triclusters, `False` (default) if using biclusters
    """
    mats = pd.read_csv(baseline_file)
    mats = mats.fillna(0)
    y_res = mats["Evolution"]
    mats.drop(columns=["Evolution"], inplace=True)
    X_res = mats.loc[:, ].values

    n_feats = int((len(mats.columns)-1)/n)

    ps_tr = list()
    for e in X_res:
        p_tric = Tricluster(n, n_feats, 1)
        i = 0
        f = 0

        for v in e[1:]:
            p_tric.addValue("T-"+str(i), "S-"+str(f), "G-" + str(e[0]), v)
            if f == (n_feats-1):
                f = 0
                i += 1
            else:
                f += 1
        ps_tr.append(p_tric)

    # f_cont = ["S-5"]

    # f_cat = ["S-0", "S-1", "S-2", "S-3", "S-4", "S-6"]
    # # if not last:
    # #     f_cat.append("S-11")
    # #     f_cat.append("S-12")
    f_continuos = list(map(lambda f: f"S-{list(constants.TEMPORAL_FEATURES.keys()).index(f[0])}", filter(
        lambda x: x[1] == 'continuos', constants.TEMPORAL_FEATURES.items())))
    f_categorical = list(map(lambda f: f"S-{list(constants.TEMPORAL_FEATURES.keys()).index(f[0])}", filter(
        lambda x: x[1] == 'categorical', constants.TEMPORAL_FEATURES.items())))

    compute_and_write_matrices(
        path_trics, path_matr, ps_tr, f_categorical, f_continuos, y_res, tri=tri, bic=False)

    # Path(path_matr).mkdir(parents=True, exist_ok=True)
    # directory = os.fsencode(path_trics)
    # print('Distance Matrices')
    # for file in sorted(os.listdir(directory)):
    #     filename = os.fsdecode(file)
    #     if filename.endswith(".txt") and not filename.endswith("summary.txt"):
    #         print(filename)
    #         triclusters = get_triclusters(path_trics + '/' + filename)
    #         bin_matrix, cols = compute_similarity_matrix(ps_tr, triclusters, f_cat, f_cont, corr=False, tri=tri)
    #         write_matrix(bin_matrix, path_matr + filename[:-4] + "_DistanceMatrix.csv", y_res, cols)

    # directory = os.fsencode(path_trics)
    # print('Correlation Matrices')
    # for file in sorted(os.listdir(directory)):
    #     filename = os.fsdecode(file)
    #     if filename.endswith(".txt") and not filename.endswith("summary.txt"):
    #         print(filename)
    #         triclusters = get_triclusters(path_trics + '/' + filename)
    #         bin_matrix, cols = compute_similarity_matrix(ps_tr, triclusters, f_cat, f_cont, corr=True, tri=tri)
    #         write_matrix(bin_matrix, path_matr + filename[:-4] + "_CorrelationMatrix.csv", y_res, cols)


# Auxiliar

def compute_and_write_matrices(path, path_matr, ps_tr, f_cat, f_cont, y_res, tri, bic=False):

    Path(path_matr).mkdir(parents=True, exist_ok=True)
    directory = os.fsencode(path)
    print('Distance Matrices')
    for file in sorted(os.listdir(directory)):
        filename = os.fsdecode(file)
        if filename.endswith(".txt") and not filename.endswith("summary.txt"):
            print(filename)
            triclusters = get_triclusters(path + '/' + filename)
            bin_matrix, cols = compute_similarity_matrix(
                ps_tr, triclusters, f_cat, f_cont, corr=False, tri=tri)
            if bic:
                cols = list(map(lambda p: p.replace("Tric", "Bic"), cols))
            write_matrix(bin_matrix, path_matr +
                         filename[:-4] + "_DistanceMatrix.csv", y_res, cols)

    directory = os.fsencode(path)
    print('Correlation Matrices')
    for file in sorted(os.listdir(directory)):
        filename = os.fsdecode(file)
        if filename.endswith(".txt") and not filename.endswith("summary.txt"):
            print(filename)
            triclusters = get_triclusters(path + '/' + filename)
            bin_matrix, cols = compute_similarity_matrix(
                ps_tr, triclusters, f_cat, f_cont, corr=True, tri=tri)
            if bic:
                cols = list(map(lambda p: p.replace("Tric", "Bic"), cols))
            write_matrix(bin_matrix, path_matr +
                         filename[:-4] + "_CorrelationMatrix.csv", y_res, cols)
