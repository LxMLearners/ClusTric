import sys
import subprocess
from pathlib import Path
from tabfiles_utils import tab_from_baseline_temporal
from compute_similarities import sim_matrix_tric
import constants

if __name__ == "__main__":

    n = int(sys.argv[1])
    constants.get_config(sys.argv[2])
    group = ""
    top_folder = constants.TOP_FOLDER

    baseline_folder_temporal = constants.BASELINE_DIR_T + "{}TPS/".format(n)

    tabfiles_folder = constants.TAB_DIR + "{}TPS/".format(n)
    triclusters_folder = constants.TRICLUSTERS_DIR + \
        "{}TPS{}".format(n, group)

    matrices_folder_t = constants.MATRICES_DIR_T + \
        "{}TPS{}/".format(n, group)
    results_folder = constants.RESULTS_DIR
    results_folder_tric = constants.RESULTS_DIR_TRIC

    Path(tabfiles_folder).mkdir(parents=True, exist_ok=True)
    Path(triclusters_folder).mkdir(parents=True, exist_ok=True)
    Path(matrices_folder_t).mkdir(parents=True, exist_ok=True)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    Path(results_folder_tric).mkdir(parents=True, exist_ok=True)

    print("*** Starting Triclustering Process ***")
    print("Parameters: ")

    print("Number of Consecutive Snapshots:", n)

    print()
    print("Creating TAB FILES")

    out_tri = "{}TAB_{}TPS{}.tab".format(tabfiles_folder, n, group)
    baseline_file_temporal = baseline_folder_temporal + \
        "{}TPS{}_baseline_temporal.csv".format(n, group)

    nP = tab_from_baseline_temporal(
        baseline_file_temporal, out_tri, n, list(constants.TEMPORAL_FEATURES.keys()))
    nP = int(nP*0.1)

    print("TAB FILES created!")

    # 1. Run TCtriCluster (biclusters and triclusters)
    # 1.1 Triclustering first
    print("Running Triclustering")
    cmd = "python3 src/triclustering/TCtriCluster.py -f {}TAB_{}TPS{}.tab -sT 2 -sS 1 -sG {} -w 0.01 -o 1 -mv 0.5 > {}/out_1.txt".format(
        tabfiles_folder, n, group, nP, triclusters_folder)
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Triclustering Completed")

    path_trics = top_folder + "triclusters/{}TPS".format(n)
    path_matr_t = top_folder + "matrices/trics/{}TPS/".format(n)

    sim_matrix_tric(n, baseline_file_temporal, triclusters_folder, matrices_folder_t,
                    last=False)  # tri - Patterns 3D
