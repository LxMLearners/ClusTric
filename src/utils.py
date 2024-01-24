import sys
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from scipy.cluster.hierarchy import ward, fcluster, cophenet
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import os
import errno
import pacmap
from itertools import zip_longest
import triclustering.constants as constants
from pathlib import Path

def parse_data():
    try:
        n = constants.MIN_APP
        matrix_name = constants.MATRICES_DIR_T + "{}TPS/out_1_DistanceMatrix.csv".format(n)
        snapshots_name = constants.DATA_FILE # snapshots file
        
        Path(constants.TRAJECTORY_DIR).mkdir(parents=True, exist_ok=True) 
        Path(constants.VISUALIZATION_DIR).mkdir(parents=True, exist_ok=True) 
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

    data = pd.read_csv(matrix_name)
    data.drop(data.columns[[-1,]], axis=1, inplace=True) # drop evolution column
    n_trics = len(data.columns) + 1
    patients = pd.read_csv(snapshots_name)

    # remove patietns with less than MIN_APP appointments
    counts = patients[constants.REF_FEATURE].value_counts()
    mask = counts >= constants.MIN_APP
    filtered_patients = patients[patients[constants.REF_FEATURE].isin(counts[mask].index)]
    filtered_patients = filtered_patients.groupby(constants.REF_FEATURE).first().reset_index()
    
    return data, patients,filtered_patients

def get_color_list():
    colormap = plt.cm.get_cmap('rainbow')
    colors = [colormap(i/constants.N_CLUST ) for i in range(constants.N_CLUST )]
    return colors

def tsne(tempData, labels):
    """
    Computes the tsne dimensionality reduction 

    Parameters
    ----------
    tempData: data to plot
    labels: labels of each patient in data
    output_name: name of the output folder
    """
    colors = get_color_list()
    tsne = TSNE(n_components=2, random_state=0)
    X_2d = tsne.fit_transform(tempData)

    new = tempData.copy()
    new['tsne-2d-one'] = X_2d[:,0]
    new['tsne-2d-two'] = X_2d[:,1]
    fig = plt.figure(figsize=(16,10))
    sns.scatterplot(
        x = "tsne-2d-one", y = "tsne-2d-two",
        hue = labels,
        palette = colors,
        data = new,
        legend = "full"
)
    fig.savefig(constants.VISUALIZATION_DIR + 'tsne.pdf')

def pacmap_func(tempData, labels):
    """
    Computes the pacmap dimensionality reduction 

    Parameters
    ----------
    tempData: data to plot
    labels: labels of each patient in data
    output_name: name of the output folder
    """
    colors = get_color_list()
    embedding = pacmap.PaCMAP(n_components=2, random_state=0)
    X_2d = embedding.fit_transform(tempData.values)

    new = tempData.copy()
    new['pacmap-2d-one'] = X_2d[:,0]
    new['pacmap-2d-two'] = X_2d[:,1]

    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x = "pacmap-2d-one", y = "pacmap-2d-two",
        hue = labels,
        palette = colors,
        data = new,
        legend = "full"
    )

    plt.savefig(constants.VISUALIZATION_DIR + 'pacmap.pdf')

def hierarchical_clustering(data):
    """
    Computes the agglomerative clustering 

    Parameters
    ----------
    data: data to cluster
    """
    plt.ylabel('distance')
    clusters = shc.linkage(data, method="ward", metric="euclidean")

    shc.dendrogram(clusters, p = 20, truncate_mode = 'lastp', # show only the last p merged clusters
                show_leaf_counts = False) 
    plt.gcf()
    plt.savefig(constants.TOP_FOLDER + '/dendrogram.pdf')

    Ward_model = AgglomerativeClustering(n_clusters= constants.N_CLUST, metric='euclidean', linkage='ward')
    Ward_model.fit(data)

    print('Silhouette Score: ', silhouette_score(data, Ward_model.labels_, metric='euclidean'))
    print('Calinski Harabasz Score: ', calinski_harabasz_score(data, Ward_model.labels_))
    print('Davies Bouldin Score: ', davies_bouldin_score(data, Ward_model.labels_))

    return Ward_model.labels_


def slope(x1, y1, x2, y2):
    m = (y2-y1)/(x2-x1)
    return m

def format_mogp_axs(ax, max_x=8, x_step=1.0, y_label=[0,24,48], y_minmax=(-3, 53)):
    ax.set_xlim([0, max_x])
    ax.set_xticks(np.arange(0, max_x + 1, x_step))
    ax.set_yticks(y_label)
    ax.set_ylim(y_minmax)
    return ax


def simple_trajectories(clusters):
    """
    Computes the trajectories of each clustering in the temporal features

    Parameters
    ----------
    clusters: list of n_clust lists with each list comprising the snapshots of the patients in the corresponding cluster
    """
    colors = get_color_list()

    for feature in list(constants.TEMPORAL_FEATURES.keys()):
        fig, ax = plt.subplots()

        max_val =0
        for j in range(constants.N_CLUST):
            prg = clusters[j].groupby(constants.REF_FEATURE)[feature]           

            lst = []
            for _, group in prg:
                lst.append(group.values)
                
            max_len = max(len(i) for i in lst)
            transposed = [[i[o] for i in lst if len(i) > o] for o in range(max_len)]
            #transposed = [list(filter(None,i)) for i in zip_longest(*lst)]

            n_samples = []
  
            means = []
            ci = []
            max_val = 0
            for values in transposed:
                if values != []:
                    if max_val < np.nanmax(values):
                        max_val = np.nanmax(values)
                    n_samples.append(len(values))
                    ci = np.append(ci, 1.96 * np.nanstd(values)/np.sqrt(len(values)))
                    means = np.append(means,    np.nanmean(values))
            app = np.arange(0,len(means))
            #num_pat = 'n = {}'.format(len(clusters[j].groupby('REF')))
            
            ax.plot(app, means, marker = '.', color = colors[j], label = 'Cluster ' + str(j+1))
            ax.fill_between(app, (means-ci), (means+ci), color=colors[j], alpha=.1)

        
            slope_value=[]
            for i in range(1,6):
                v=slope(i-1,means[i-1], i, means[i])
                slope_value.append(v)
                plt.text( i-0.5 , (means[i] + means[i-1])/2, str(round(v,2)), fontsize=5, color = colors[j])
                
                plt.text(i, means[i] + 0.001, str(n_samples[i-1]), fontsize=8, color = colors[j], fontweight= 'bold')

        format_mogp_axs(ax, 5, 1, y_label=[0,max_val/2,max_val], y_minmax = (0, max_val+1))

        plt.xlabel("Appointments")
        plt.ylabel(str(feature))
        plt.legend()
        fig.savefig(constants.TRAJECTORY_DIR + str(feature) + '.pdf')