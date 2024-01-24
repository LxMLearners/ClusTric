from utils import hierarchical_clustering, tsne, pacmap_func,simple_trajectories,parse_data
import pandas as pd 
import subprocess
import triclustering.constants as constants
import sys

if __name__ == "__main__":
    constants.get_config(sys.argv[1])
    n = constants.MIN_APP
    #create longitudinal tables
    cmd = "python3 src/triclustering/longitudinal_tables_strat.py {} {}".format(n, sys.argv[1]) 
    print(cmd)
    subprocess.call(cmd, shell=True)

    #create triclustering representations
    cmd = "python3 src/triclustering/main.py {} {}".format(n,sys.argv[1])
    print(cmd)
    subprocess.call(cmd, shell=True)   

    #find the clusters in data
    print('*** CLUSTERING ***')
    data, patients, filtered_patients = parse_data()
    labels = hierarchical_clustering(data)


    #visualize the representations
    tsne(data, labels)
    pacmap_func(data,labels)

    df = pd.DataFrame() 
    df['Patient_ID'] = filtered_patients[constants.REF_FEATURE].values
    df['Labels'] = labels
    df.to_csv(constants.TOP_FOLDER + '/labels.csv') 

    #plot the trajectories of each cluster
    patients = pd.merge(patients, df[['Patient_ID', 'Labels']].rename(columns={'Patient_ID':constants.REF_FEATURE}), on = constants.REF_FEATURE)
    patients = patients[patients['Labels'].notna()]
    clusters = []
    for i in range(constants.N_CLUST):
        clusters.append(patients.loc[patients['Labels'] == i])
    simple_trajectories(clusters) 