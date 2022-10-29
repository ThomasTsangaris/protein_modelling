# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:34:51 2022

@author: Thomas Tsangaris
"""

""" Purpose: perform UMAP on an input ensemble of protein structures to visualize
the conformers reduced to 2-dimensions. This program uses an inter-residue pairwise
distance Euclidean metric for determining similarity of conformers.
"""
import umap
import MDAnalysis as md
from MDAnalysis.analysis import distances
import time
import glob
import numpy as np
import matplotlib.pyplot as plt


ti=time.perf_counter()

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

umap_2d = umap.UMAP(n_components=2, n_neighbors=100000).fit_transform(X)

plt.scatter(umap_2d[:,0], umap_2d[:,1], marker='.', color='blue')
plt.show()

print(np.shape(X))

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")