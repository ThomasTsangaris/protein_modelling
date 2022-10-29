# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 20:47:31 2022

@author: Thomas Tsangaris
"""

"""Purpose: input the number of diffusion components and reduce the dimension
of your data using the diffusion maps framework. Then, cluster in the reduced
dimensions using agglomerative hierarchical clustering.

"""
import MDAnalysis as md
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import time
from sklearn.cluster import DBSCAN
from sklearn.metrics import davies_bouldin_score, calinski_harabasz_score, silhouette_score
from yellowbrick.cluster import KElbowVisualizer
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import pydiffmap
from pydiffmap import diffusion_map as dm
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator
import os
import shutil


ti=time.perf_counter()

seed = 42
np.random.seed(seed)

# determine this number using a different script
# N_DC = 516 # number of diffusion components
# N_DC = 4997 # number of diffusion components
N_DC = 14
# N_DC = 4

def inter_res_pairwise_euclidean_dis_metric(Ca_ens1: np.array, Ca_ens2: np.array):
    """Return the Euclidean distance between mean Ca-CA pairwise distance 
    matrices.
    
    """
    return np.linalg.norm(Ca_ens1 - Ca_ens2)
    

# TraDES 100% RC N = 5000 (April 2022)
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_32_91/'
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/'
# FFT NP N = 5000 (April 2022) 
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/'
# FFT NP with canonical binding helix region
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/ss_region_ensemble/'
# FFT NP without canonical binding helix region
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/'
# FFT 5P N = 5000 (April 2022)
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'

# BME FFT Weight Sampled Ensembles (April 2022) (N=5000):
# NP
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/'
# 5P 
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/'

files=glob.glob(mypath+"*.pdb")
u_ens = md.Universe(files[0], files, in_memory=True)

# number of alpha carbons = number of residues
N = len(u_ens.residues)
# number of frames in trajectory
M=len(files)

# initialize Ca pairwise distance matrix
X = np.zeros((M, int((N*(N-1))/2)))

i = 0
for frame in u_ens.trajectory: # loop through each frame
    CA_atoms = u_ens.select_atoms('name CA') # get Ca pairwise distances
    # pairwise Ca distances for a frame
    CA_distances = distances.self_distance_array(CA_atoms.positions)
    X[i, :] = CA_distances # append to array
    i += 1 # increment index
# now, entry i in X contains all Ca pairwise distances for the ith pdb file

epsilon = 0.05
n_evecs = N_DC
alpha = 0.5
nearest_neighbours = 5000
# nearest_neighbours = 500
# nearest_neighbours = 64

mydmap = dm.DiffusionMap.from_sklearn(n_evecs = n_evecs, epsilon = 'bgh', alpha = alpha, k = nearest_neighbours, metric = inter_res_pairwise_euclidean_dis_metric)
dmap = mydmap.fit_transform(X)

evecs = mydmap.evecs
evals = mydmap.evals

k_tup = (0.5, 20) # eps range to test

print(dmap)
print(np.shape(dmap))


X = dmap # use diffusion components for similarity measure in clustering

# simple code, but the following is from here: https://python-bloggers.com/2021/06/davies-bouldin-index-for-k-means-clustering-evaluation-in-python/
sil_results = {}
wss_results = {}
cal_results = {}
dbi_results = {}

for i in np.arange(min(k_tup), max(k_tup), 0.5):
    np.random.seed(seed)
    dbscan = DBSCAN(eps = i, min_samples = (2*N_DC))
    labels = dbscan.fit_predict(X)
    print(max(labels))
    print(list(labels).count(0))
    print(len(labels))
    db_index = davies_bouldin_score(X, labels)
    sil_index = silhouette_score(X, labels)
    cal_index = calinski_harabasz_score(X, labels)
    
    dbi_results.update({i: db_index})
    sil_results.update({i: sil_index})
    cal_results.update({i: cal_index})

plt.figure()
plt.plot(list(dbi_results.keys()), list(dbi_results.values()), linestyle='--', marker='o', color='b')
plt.xlabel("$\epsilon$")
plt.ylabel("davies-bouldin index")

plt.figure()
plt.plot(list(sil_results.keys()), list(sil_results.values()), linestyle='--', marker='o', color='b')
plt.xlabel("$\epsilon$")
plt.ylabel("mean silhouette score")

plt.figure()
plt.plot(list(cal_results.keys()), list(cal_results.values()), linestyle='--', marker='o', color='b')
plt.xlabel("$\epsilon$")
plt.ylabel("calinski harabasz score")


# reference for below: http://sefidian.com/2020/12/18/how-to-determine-epsilon-and-minpts-parameters-of-dbscan-clustering/#:~:text=For%202-dimensional%20data%2C%20use%20DBSCAN%E2%80%99s%20default%20value%20of,value%2C%20you%20can%20move%20on%20to%20determining%20%CE%B5.
# for determining epsilon (knee on the following plot):
neighbors = NearestNeighbors(n_neighbors=(N_DC*2))
neighbors_fit = neighbors.fit(X)
distances, indices = neighbors_fit.kneighbors(X)

distances = np.sort(distances, axis=0)
distances = distances[:,1]
print(distances[4980:])
kneedle = KneeLocator(range(len(distances)), distances, S=1.0, curve="convex", direction="increasing")
print("Knee occurs at x=", kneedle.knee, "and y=", kneedle.knee_y)

plt.figure()
plt.plot(distances)
plt.axvline(x=kneedle.knee, linestyle='--', color='k')
plt.xlabel("Unordered Datapoints")
plt.ylabel("Diffusion Distance from Nearest Neighbor")
plt.xlim([4950, 5000])
plt.ylim([-10, 100])



# epsilon chosen from above analysis:
dbscan = DBSCAN(eps = 4.279, min_samples = (2*N_DC))
labels = dbscan.fit_predict(X)

print("THE NUMBER OF CLUSTERS IS:", max(labels) + 1)

# now store cluster ensembles in directories
mapping = zip(files, labels)
# a list of tuples containing two entries: path of conformer and cluster 
# assignment (1 or 2 for two clusters)
pdb_to_cluster_num = list(mapping)

cluster_dir_name = "dbscan_4_dc_cluster_"

pdb_frame_num = 0 # keep track of frame number

# IMPORTANT NOTE: be careful here; make sure you empty ensemble cluster 
# directories before running so you don't double count
for i in range(len(pdb_to_cluster_num)):
    cluster_num = pdb_to_cluster_num[i][1]
    
    output_folder = mypath + cluster_dir_name + str(cluster_num)

    if (os.path.exists(output_folder)): # if the cluster dir exists
        # add to it
        shutil.copy(files[pdb_frame_num], output_folder + "/" + \
                    str(pdb_frame_num + 1) + ".pdb")

    if not os.path.exists(output_folder): # if cluster dir does not exist
        # create the appropriate cluster dir
        os.makedirs(output_folder)
        shutil.copy(files[pdb_frame_num], output_folder + "/" + \
                    str(pdb_frame_num + 1) + ".pdb")
    
    pdb_frame_num += 1


# min_clust_size_lst = [] # list of minimum cluster sizes per number of clusters
# for k in range(min(k_tup), max(k_tup)): # loop through each number of clusters
#     # run agglomerative hierarchical clustering
#     linkage = hierarchy.linkage(X, method="ward")
#     # extract k clusters from dendrogram
#     clusters = fcluster(linkage, k, criterion='maxclust')

#     clust_size_dict = {} # dict of cluster number to size of cluster
#     for cluster_assign in clusters: # loop through each cluster assignment
#          if cluster_assign not in clust_size_dict:
#              clust_size_dict[cluster_assign] = 0
#          clust_size_dict[cluster_assign] += 1
#     min_clust_size_lst.append(min(clust_size_dict.values())) # append minimum

# plt.figure()
# plt.plot(range(min(k_tup), max(k_tup)), min_clust_size_lst, linestyle='--', marker='o', color='b')
# plt.xlabel("k")
# plt.ylabel("minimum cluster size")

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")