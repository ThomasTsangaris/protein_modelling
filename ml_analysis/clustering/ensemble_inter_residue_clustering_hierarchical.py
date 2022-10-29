# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 02:55:44 2022

@author: Thomas Tsangaris
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster
import sklearn
import sys
import glob
import time
import MDAnalysis as md
from MDAnalysis.analysis import distances
import os
import shutil
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import davies_bouldin_score
# from yellowbrick.cluster import KElbowVisualizer

""" Purpose: perform hierarchical clustering of protein structures based on
 inter-residue distances by inputting the number of clusters desired.
 The distance metric employed is Euclidean.
"""

ti=time.perf_counter()

# Input your Ensemble Path
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
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'

# BME FFT Weight Sampled Ensembles (April 2022) (N=5000):
# NP
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/'
# 5P 
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/'


files=glob.glob(mypath+"*.pdb")[:]
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

# distance cutoff for clustering
# max_dis = 24000
max_dis = 22780

# number of clusters
k = 4

# run agglomerative hierarchical clustering
linkage = hierarchy.linkage(X, method="ward") # should try multiple methods and see differences
# print("LINKAGE:", linkage[:,2])

cutoffs = [] # list of cutoffs
# now, loop through the linkage matrix and plot cutoff distance vs number of clusters
for row in linkage:
    cutoff = row[2] # cutoff distance
    num_confs = row[3] # number of conformers that make up said cluster
    cutoffs.append(cutoff)

cutoffs = np.array(cutoffs)

# visualize the number of clusters determined by a biophysical cutoff distance
plt.figure()
cutoffs_sorted = -np.sort(-cutoffs) # 1st entry corresponds with top most cluster
plt.plot(cutoffs_sorted/7260, range(1, len(cutoffs_sorted)+1))
plt.xlabel('Normalized Cutoff Distance (Ang)', fontsize=12)
plt.ylabel('Number of Clusters', fontsize=12)
plt.ylim(0,12)
    

# extract k clusters from dendrogram
clusters = fcluster(linkage, k, criterion='maxclust')

""" for coloring each cluster on a plot, only makes sense for 
low-dimeensional data 
"""
# # plot k clusters
# # plt.figure(figsize=(10, 8))
# # plt.scatter(X[:,0], X[:,1], c=clusters, cmap='prism')  # plot points with cluster dependent colors
# # plt.show()

# # plot the dendrogram
plt.figure(figsize=(20,6))
dendrogram = hierarchy.dendrogram(linkage, color_threshold=max_dis, orientation="top",leaf_font_size=9, leaf_rotation=360)
# dendrogram = hierarchy.dendrogram(linkage, color_threshold=11500, orientation="top",leaf_font_size=9, leaf_rotation=360)
plt.axhline(max_dis)
plt.ylabel('Euclidean Distance', fontsize=14)
plt.xlabel('Conformers', fontsize=14)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
# save figure as a PNG file
# plt.savefig('C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/spencer_2022_dendrogram_NP.png', dpi=300, transparent=False, bbox_inches='tight')

print("clusters array:", clusters)
print("clusters array shape:", np.shape(clusters))

# # now store cluster ensembles in directories
# mapping = zip(files, clusters)
# # a list of tuples containing two entries: path of conformer and cluster 
# # assignment (1 or 2 for two clusters)
# pdb_to_cluster_num = list(mapping)

# cluster_dir_name = "ip_bio_cluster_check_"

# pdb_frame_num = 0 # keep track of frame number

# # IMPORTANT NOTE: be careful here; make sure you empty ensemble cluster 
# # directories before running so you don't double count
# for i in range(len(pdb_to_cluster_num)):
#     cluster_num = pdb_to_cluster_num[i][1]
    
#     output_folder = mypath + cluster_dir_name + str(cluster_num)

#     if (os.path.exists(output_folder)): # if the cluster dir exists
#         # add to it
#         shutil.copy(files[pdb_frame_num], output_folder + "/" + \
#                     str(pdb_frame_num + 1) + ".pdb")

#     if not os.path.exists(output_folder): # if cluster dir does not exist
#         # create the appropriate cluster dir
#         os.makedirs(output_folder)
#         shutil.copy(files[pdb_frame_num], output_folder + "/" + \
#                     str(pdb_frame_num + 1) + ".pdb")
    
#     pdb_frame_num += 1


# k_tup = (2, 12) # cluster range to test


 # Use Various Statistical Metrics to Determine Optimal Cluster Number
# plt.figure()
# model = AgglomerativeClustering()
# visualizer = KElbowVisualizer(model, k=k_tup, timings= False, metric='silhouette')
# visualizer.fit(X)
# plt.ylabel('<silhouette score>')
# visualizer.show()

# plt.figure()
# model = AgglomerativeClustering()
# visualizer = KElbowVisualizer(model, k=k_tup, timings= False, metric='calinski_harabasz')
# visualizer.fit(X)
# plt.ylabel('calinski-harabasz index')
# visualizer.show()

# plt.figure()
# model = AgglomerativeClustering()
# visualizer = KElbowVisualizer(model, k=k_tup, timings= False, metric='distortion')
# visualizer.fit(X)
# plt.ylabel('<sum of squared distances to centers>')
# visualizer.show()

# # simple code, but the following is from here: https://python-bloggers.com/2021/06/davies-bouldin-index-for-k-means-clustering-evaluation-in-python/
# results = {}
# for i in range(min(k_tup), max(k_tup)):
#     kmeans = AgglomerativeClustering(n_clusters=i)
#     labels = kmeans.fit_predict(X)
#     db_index = davies_bouldin_score(X, labels)
#     results.update({i: db_index})
# plt.figure()
# plt.plot(list(results.keys()), list(results.values()), linestyle='--', marker='o', color='b')
# plt.xlabel("k")
# plt.ylabel("davies-bouldin index")


# min_clust_size_lst = [] # list of minimum cluster sizes per number of clusters
# for k in range(min(k_tup), max(k_tup)): # loop through each number of clusters
#     # run agglomerative hierarchical clustering
#     linkage = hierarchy.linkage(X, method="ward")
#     # extract k clusters from dendrogram
#     clusters = fcluster(linkage, k, criterion='maxclust')
    
#     clust_size_dict = {} # dict of cluster number to size of cluster
#     for cluster_assign in clusters: # loop through each cluster assignment
#           if cluster_assign not in clust_size_dict:
#               clust_size_dict[cluster_assign] = 0
#           clust_size_dict[cluster_assign] += 1
#     min_clust_size_lst.append(min(clust_size_dict.values())) # append minimum

# plt.figure()
# plt.plot(range(min(k_tup), max(k_tup)), min_clust_size_lst, linestyle='--', marker='o', color='b')
# plt.xlabel("k")
# plt.ylabel("minimum cluster size")
# plt.ylim([-10, 2])

# for k in range(k_tup):
    
#     # run agglomerative hierarchical clustering
#     linkage = hierarchy.linkage(X, method="ward") # should try multiple methods and see differences
#     print("LINKAGE:", linkage[:,2])
    
#     # extract k clusters from dendrogram
#     clusters = fcluster(linkage, k, criterion='maxclust')
    
#     # # plot k clusters
#     # # plt.figure(figsize=(10, 8))
#     # # plt.scatter(X[:,0], X[:,1], c=clusters, cmap='prism')  # plot points with cluster dependent colors
#     # # plt.show()
    
#     # # plot the dendrogram
#     plt.figure(figsize=(20,6))
#     dendrogram = hierarchy.dendrogram(linkage, color_threshold=19500, orientation="top",leaf_font_size=9, leaf_rotation=360)
#     # dendrogram = hierarchy.dendrogram(linkage, color_threshold=11500, orientation="top",leaf_font_size=9, leaf_rotation=360)
#     plt.axhline(19500)
#     plt.ylabel('Euclidean Distance', fontsize=14)
#     plt.xlabel('Conformers', fontsize=14)
#     plt.tick_params(
#         axis='x',          # changes apply to the x-axis
#         which='both',      # both major and minor ticks are affected
#         bottom=False,      # ticks along the bottom edge are off
#         top=False,         # ticks along the top edge are off
#         labelbottom=False) # labels along the bottom edge are off
#     # save figure as a PNG file
#     # plt.savefig('C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/spencer_2022_dendrogram_NP.png', dpi=300, transparent=False, bbox_inches='tight')
    
#     print("clusters array:", clusters)
#     print("clusters array shape:", np.shape(clusters))
     
#     # now store cluster ensembles in directories
#     mapping = zip(files, clusters)
#     # a list of tuples containing two entries: path of conformer and cluster 
#     # assignment (1 or 2 for two clusters)
#     pdb_to_cluster_num = list(mapping)
    
#     cluster_dir_name = str(k) + "_bio_cluster_"
    
#     pdb_frame_num = 0 # keep track of frame number
    
#     # IMPORTANT NOTE: be careful here; make sure you empty ensemble cluster 
#     # directories before running so you don't double count
#     for i in range(len(pdb_to_cluster_num)):
#         cluster_num = pdb_to_cluster_num[i][1]
        
#         output_folder = mypath + cluster_dir_name + str(cluster_num)
    
#         if (os.path.exists(output_folder)): # if the cluster dir exists
#             # add to it
#             shutil.copy(files[pdb_frame_num], output_folder + "/" + \
#                         str(pdb_frame_num + 1) + ".pdb")
    
#         if not os.path.exists(output_folder): # if cluster dir does not exist
#             # create the appropriate cluster dir
#             os.makedirs(output_folder)
#             shutil.copy(files[pdb_frame_num], output_folder + "/" + \
#                         str(pdb_frame_num + 1) + ".pdb")
        
#         pdb_frame_num += 1

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")
