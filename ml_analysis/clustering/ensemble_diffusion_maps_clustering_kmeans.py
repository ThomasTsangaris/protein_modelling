# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 20:44:11 2022

@author: Thomas Tsangaris
"""

"""Purpose: input the number of diffusion components and reduce the dimension
of your data using the diffusion maps framework. Then, cluster in the reduced
dimensions using K-Means. Evaluate resulting clusters using various metrics.

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


ti=time.perf_counter()

seed = 42
np.random.seed(seed)

# determine this number using a different script
N_DC = 14 # number of diffusion components

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
nearest_neighbours = 4999

mydmap = dm.DiffusionMap.from_sklearn(n_evecs = n_evecs, epsilon = 'bgh', alpha = alpha, k = nearest_neighbours, metric = inter_res_pairwise_euclidean_dis_metric)
dmap = mydmap.fit_transform(X)

evecs = mydmap.evecs
evals = mydmap.evals

k_tup = (2, 30) # cluster range to test

print(dmap)
print(np.shape(dmap))


X = dmap # use diffusion components for similarity measure in clustering

plt.figure()
model = KMeans()
visualizer = KElbowVisualizer(model, k=k_tup, timings= True, metric='silhouette', random_state=seed)
visualizer.fit(X)
plt.ylabel('<silhouette score>')
visualizer.show()

plt.figure()
model = KMeans()
visualizer = KElbowVisualizer(model, k=k_tup, timings= True, metric='calinski_harabasz', random_state=seed)
visualizer.fit(X)
plt.ylabel('calinski-harabasz index')
visualizer.show()

plt.figure()
model = KMeans()
visualizer = KElbowVisualizer(model, k=k_tup, timings= True, metric='distortion', random_state=seed)
visualizer.fit(X)
plt.ylabel('<sum of squared distances to centers>')
visualizer.show()

# simple code, but the following is from here: https://python-bloggers.com/2021/06/davies-bouldin-index-for-k-means-clustering-evaluation-in-python/
results = {}
for i in range(min(k_tup), max(k_tup)):
    np.random.seed(seed)
    kmeans = KMeans(n_clusters=i, random_state=seed)
    labels = kmeans.fit_predict(X)
    db_index = davies_bouldin_score(X, labels)
    results.update({i: db_index})
plt.figure()
plt.plot(list(results.keys()), list(results.values()), linestyle='--', marker='o', color='b')
plt.xlabel("k")
plt.ylabel("davies-bouldin index")

np.random.seed(seed)
# Gap Statistic for K means
def optimalK(data, nrefs=3, maxClusters=15):
    """ from https://towardsdatascience.com/cheat-sheet-to-implementing-7-methods-for-selecting-optimal-number-of-clusters-in-python-898241e1d6ad
    Calculates KMeans optimal K using Gap Statistic 
    Params:
        data: ndarry of shape (n_samples, n_features)
        nrefs: number of sample reference datasets to create
        maxClusters: Maximum number of clusters to test for
    Returns: (gaps, optimalK)
    """
    gaps = np.zeros((len(range(1, maxClusters)),))
    resultsdf = pd.DataFrame({'clusterCount':[], 'gap':[]})
    for gap_index, k in enumerate(range(1, maxClusters)):
    # Holder for reference dispersion results
        refDisps = np.zeros(nrefs)
    # For n references, generate random sample and perform kmeans getting resulting dispersion of each loop
        for i in range(nrefs):
            
            np.random.seed(seed)
            
            # Create new random reference set
            randomReference = np.random.random_sample(size=data.shape)
            
            # Fit to it
            km = KMeans(k)
            km.fit(randomReference)
            
            refDisp = km.inertia_
            refDisps[i] = refDisp
    # Fit cluster to original data and create dispersion
        km = KMeans(k)
        km.fit(data)
            
        origDisp = km.inertia_
    # Calculate gap statistic
        gap = np.log(np.mean(refDisps)) - np.log(origDisp)
    # Assign this loop's gap statistic to gaps
        gaps[gap_index] = gap
            
        resultsdf = resultsdf.append({'clusterCount':k, 'gap':gap}, ignore_index=True)
    return (gaps.argmax() + 1, resultsdf)


plt.figure()
score_g, df = optimalK(X, nrefs=5, maxClusters=max(k_tup))
plt.plot(df['clusterCount'], df['gap'], linestyle='--', marker='o', color='b')
plt.xlabel(r'k')
plt.ylabel(r'Gap Statistic')
plt.show()

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")