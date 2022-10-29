# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 00:41:35 2022

@author: Thomas Tsangaris
"""

"""Purpose: Plot all the diffusion map 2d plots to determine which eigenvalues
can be ignored to reduce dimensionality.

"""
import MDAnalysis as md
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import time
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import davies_bouldin_score
from yellowbrick.cluster import KElbowVisualizer
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import pydiffmap
from pydiffmap import diffusion_map as dm
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster

# Reference for diffusion maps: Ronald R. Coifman and Stéphane Lafon. Diffusion maps. Applied and Computational Harmonic Analysis, 21(1):5–30, July 2006. 02271. doi:10.1016/j.acha.2006.04.006.

# Reference for module: https://pydiffmap.readthedocs.io/en/master/jupyter%20notebook%20tutorials/Metrics/Metrics.html

ti=time.perf_counter()


def inter_res_pairwise_euclidean_dis_metric(Ca_ens1: np.array, Ca_ens2: np.array):
    """Euclidean distance between inter-residue C alpha pairwise distance
    matrices.
    
    """
    return np.linalg.norm(Ca_ens1 - Ca_ens2)
    


# TraDES 100% RC N = 5000 (April 2022)
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_32_91/'
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
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/'
# 5P 
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/'

# BME FFT Weight Sampled Ensembles (April 2022) (N=10000):
# NP
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50_10000/'
# 5P 
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1_10000/'


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

N_DC = 4 # number of diffusion components to compute

epsilon = 0.05
n_evecs = N_DC
alpha = 0.5
k = 5000

mydmap = dm.DiffusionMap.from_sklearn(n_evecs = n_evecs, epsilon = 'bgh', alpha = alpha, k = k, metric = inter_res_pairwise_euclidean_dis_metric)
dmap = mydmap.fit_transform(X)

# print(X)

evecs = mydmap.evecs
evals = mydmap.evals

num_clusters = 5

# # run agglomerative hierarchical clustering
# linkage = hierarchy.linkage(X, method="ward")
# # extract k clusters from dendrogram
# clusters = fcluster(linkage, num_clusters, criterion='maxclust')
# print("cluster array:", clusters[:200])

# colors = sns.color_palette("hls", num_clusters)
# cluster_colors =[]

# i = 0
# for cluster_num in clusters:
#     cluster_colors.append(colors[cluster_num-1]) # same cluster gets same color
#     i += 1

# print(cluster_colors[:100])

# plot the eigenvalue spectrum
N_evecs_visual = N_DC # number of eigenvectors (DCs) to visualize
i_j_lst = [] # keep track of i, j's already plotted
for i in range(1, N_evecs_visual+1):
    for j in range(1, N_evecs_visual+1):
        if (i, j) not in i_j_lst:
            i_j_lst.append((i, j))
            i_j_lst.append((j, i)) # order doesn't matter
            
            plt.figure()
            # plt.scatter(x=dmap[:,i-1], y=dmap[:,j-1], \
            #                    color=cluster_colors)
            plt.scatter(x=dmap[:,i-1], y=dmap[:,j-1])
            if i - 10 >= 0 and j - 10 >= 0:
                plt.xlabel(r'$\psi_{i:02d}$'.format(i=i))
                plt.ylabel(r'$\psi_{j:02d}$'.format(j=j))
            elif i - 10 >= 0:
                plt.xlabel(r'$\psi_{i:02d}$'.format(i=i))
                plt.ylabel(r'$\psi_{j:d}$'.format(j=j))
            elif j - 10 >= 0:
                plt.xlabel(r'$\psi_{i:d}$'.format(i=i))
                plt.ylabel(r'$\psi_{j:02d}$'.format(j=j))
            else:
                plt.xlabel(r'$\psi_{i:d}$'.format(i=i))
                plt.ylabel(r'$\psi_{j:d}$'.format(j=j))
            
            # save figure as a PNG file
            plt.savefig('C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_trades_rc_dc_plots/dcs_' \
                        + str(i) + '_' + str(j) + '.png', transparent=False, \
                            bbox_inches='tight')

# plt.xlim(-10, 25)
# plt.ylim(-75, 100)

plt.figure()
plt.scatter(range(1,n_evecs+1), evals)
plt.xlabel(r'i')
plt.ylabel(r'$\lambda_i$')
# plt.xticks(range(0,n_evecs+2,2))


tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")
