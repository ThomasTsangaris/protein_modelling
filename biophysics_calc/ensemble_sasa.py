# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:09:06 2022

@author: Thomas Tsangaris
"""

""" Purpose: calculate the SASA and RSASA of two inputted conformational 
ensembles. Plot the SASA and RSASA of each (with errors approximated as SEMs)
and the SASA and RSASA differences (with errors propagated appropriately).

""" 
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import time
import glob
from scipy.stats import sem
import csv


def get_weights(weights_file: str) -> list:
    """Return an array of weights, corresponding to each conformer in a pdb
    ensemble from a valid inputted path to a weights file (weights_file)
    where each row is a single floating point number, representing the weight
    of the conformer in the same order that the pdb ensemble is organized.
    
    Precondition: the weights file contains exactly one header line.
    
    """
    weights = []
    hl = 1
    
    with open(weights_file, 'r') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for j in range(hl):
            # skip header
            next(reader)
        for row in reader:
            #Read remaining lines (row at a time)
            weights.append(float(row[0]))
        f.close()
    
    return weights


def prop_error_diff(xerr:np.array, yerr:np.array) -> np.array:
    """Return the propagated errors of subtracting the two arrays x and y with
    errors contained in xerr and yerr arrays respectively.
    
    Precondition: np.shape(xerr) = np.shape(yerr) 
    """
    return np.sqrt((xerr**2)+(yerr**2))


tic = time.perf_counter()

# paths
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/' # NP
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/' # 5P

# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/'
# mypath1 = mypath2
# NP cluster 1
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_1/'
# NP cluster 2
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_2/'
# NP cluster 3
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_3/'
# 5P cluster 1
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_1/'
# 5P cluster 2
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_2/'
# 5P cluster 3
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_3/'


# Bioclusters, Clustering on Initial Pool:
# NP FFT Cluster 1
mypath2 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_1/"
# NP FFT Cluster 2
mypath1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_2/"
# NP FFT Cluster 3
# mypath1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_3/"
# NP FFT Cluster 4
# mypath1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_4/"

# 5P cluster 1
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_1/'
# 5P cluster 2
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_2/'
# 5P cluster 3
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_3/'




# Bio clusters (clusters determined with the help of a biophysical cutoff):
# NP cluster 1
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/bio_cluster_1/'
# NP cluster 2
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/bio_cluster_2/'
# NP cluster 3
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/bio_cluster_3/'
# NP cluster 4
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/bio_cluster_4/'

# 5P cluster 1
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_1/'
# 5P cluster 2
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_2/'
# 5P cluster 3
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_3/'
# 5P cluster 4
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_4/'
# 5P cluster 5
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_5/'
# 5P cluster 6
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_6/'
# 5P cluster 7
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_7/'
# 5P cluster 8
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/bio_cluster_8/'

# 121 no mutation TraDES RC (N=5000) June 2022:
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/bp2_rc/'


# weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/ss_region_ensemble/ss_region_weights.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/no_ss_region_weights.dat'
# Bioclusters, Clustering on Initial Pool:
# NP FFT Cluster 1
weights_path2 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_1/NP_FFT_cluster_1_weights.dat"
# NP FFT Cluster 2
weights_path1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_2/NP_FFT_cluster_2_weights.dat"
# NP FFT Cluster 3
# weights_path1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_3/NP_FFT_cluster_3_weights.dat"
# NP FFT Cluster 4
# weights_path1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_4/NP_FFT_cluster_4_weights.dat"
    
# 5P cluster 1
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_1/5P_FFT_cluster_1_weights.dat'
# 5P cluster 2
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_2/5P_FFT_cluster_2_weights.dat'
# 5P cluster 3
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_3/5P_FFT_cluster_3_weights.dat'



pdb_files1 = glob.glob(mypath1 + '*.pdb')
print(len(pdb_files1))
traj1 = md.load(pdb_files1)

pdb_files2 = glob.glob(mypath2 + '*.pdb')
print(len(pdb_files2))
traj2 = md.load(pdb_files2)


# for equal weights:
numfiles = len(pdb_files1) # for first ensemble, they may have different Nconf
equal_weights1 = []
for i in range(numfiles):
    equal_weights1.append(1/numfiles)
equal_weights1 = np.array(equal_weights1)

numfiles = len(pdb_files2) # for second ensemble, they may have different Nconf
equal_weights2 = []
for i in range(numfiles):
    equal_weights2.append(1/numfiles)
equal_weights2 = np.array(equal_weights2)


# extract weights from files
weights1 = get_weights(weights_path1)
# weights1 = equal_weights1
# weights2 = get_weights(weights_path2)
weights2 = equal_weights2

# phosphosites
phosphosites = [37, 46, 65, 70, 83]
# amino acid sequence
seq = "MSSSAGSGHQPSQSRAIPTRTVAISDAAQLPHDYCTTPGGTLFSTTPGGTRIIYDRKFLLDRRNSPMAQTPPCHLPNIPGVTSPGTLIEDSKVEVNNLNNLNNHDRKHAVGDDAQFEMDIC"
# MSASA dictionary
res_to_msasa = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}
# compute SASAs
sasa_compute_1 = md.shrake_rupley(traj1, mode='residue')  
    
sasa_compute_2 = md.shrake_rupley(traj2, mode='residue')  

sasa_arr1 = np.average(sasa_compute_1, weights=weights1, \
                                    axis=0)*100 # per residue
sasa_arr2 = np.average(sasa_compute_2, weights=weights2, \
                                    axis=0)*100 # per residue
    
sasa_arr1_hist = np.average(sasa_compute_1, \
                                    axis=1)*100 # per conformer
sasa_arr2_hist = np.average(sasa_compute_2, \
                                    axis=1)*100 # per conformer


print(np.shape(sasa_arr1))
print(np.shape(sasa_arr2))
# compute RSASAs
rsasa_arr1 = []
rsasa_arr2 = []
for i in range(len(seq)):
    msasa = res_to_msasa[seq[i]]
    rsasa_arr1.append(sasa_arr1[i] / msasa)
    rsasa_arr2.append(sasa_arr2[i] / msasa)
# compute SASA uncertainties (standard deviation of the mean)
sasa_sem_arr1 = sem(sasa_compute_1)*100
sasa_sem_arr2 = sem(sasa_compute_2, axis=0)*100
print(np.shape(sasa_sem_arr1))
print(np.shape(sasa_sem_arr2))

# compute RSASAs uncertainties
rsasa_sem_arr1 = []
rsasa_sem_arr2 = []
for i in range(len(seq)):
    msasa = res_to_msasa[seq[i]]
    rsasa_sem_arr1.append(sasa_sem_arr1[i] / msasa)
    rsasa_sem_arr2.append(sasa_sem_arr2[i] / msasa)

# make sure every array is a numpy array
sasa_arr1 = np.array(sasa_arr1)
sasa_arr2 = np.array(sasa_arr2)
rsasa_arr1 = np.array(rsasa_arr1)
rsasa_arr2 = np.array(rsasa_arr2)
sasa_sem_arr1 = np.array(sasa_sem_arr1)
sasa_sem_arr2 = np.array(sasa_sem_arr2)
rsasa_sem_arr1 = np.array(rsasa_sem_arr1)
rsasa_sem_arr2 = np.array(rsasa_sem_arr2)

# compute SASA differences
sasa_arr_diff = sasa_arr1 - sasa_arr2
# compute RSASA differences
rsasa_arr_diff = rsasa_arr1 - rsasa_arr2
# compute SASA differences' uncertainties
sasa_diff_unc = prop_error_diff(sasa_sem_arr1, sasa_sem_arr2)
# compute RSASA differences' uncertainties
rsasa_diff_unc = prop_error_diff(rsasa_sem_arr1, rsasa_sem_arr2)

# plot SASAs
# first ensemble
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), sasa_arr1, edgecolor='k', yerr=sasa_sem_arr1)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean SASA (Angstroms$^2$)", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')

# second ensemble
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), sasa_arr2, edgecolor='k', yerr=sasa_sem_arr2)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean SASA (Angstroms$^2$)", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')
# plot RSASAs
# first ensemble
print(np.shape(rsasa_arr1))
print(np.shape(rsasa_sem_arr1))
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), rsasa_arr1, edgecolor='k', yerr=rsasa_sem_arr1)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean Relative SASA", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')

# second ensemble
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), rsasa_arr2, edgecolor='k', yerr=rsasa_sem_arr2)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean Relative SASA", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')
# plot SASA differences
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), sasa_arr_diff, edgecolor='k', yerr=sasa_diff_unc)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean SASA Difference (Angstroms$^2$)", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')
# plot RSASA differences
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(seq)+1), rsasa_arr_diff, edgecolor='k', yerr=rsasa_diff_unc)
plt.xlabel("Residue", fontsize=16)
plt.ylabel("Mean Relative SASA Difference", fontsize=16)
for p in phosphosites:
    plt.axvline(x=p, linestyle='--', color='r')
    

# below needs to be double checked
# plot ratio of SASA 1 and SASA 2
plt.figure(figsize=(9,4)) # adjust figure size
plt.hist((sasa_arr1_hist), density=True, bins=30, weights=weights1)
plt.xlabel("SASA Ensemble Ratio", fontsize=16)
plt.ylabel("Probability Density", fontsize=16)
# plt.xlim(0, 15)
plt.xlim(80, 140)

print(np.average(sasa_arr1_hist, weights=weights1))
print(np.average(sasa_arr2_hist, weights=weights2))

plt.figure(figsize=(9,4)) # adjust figure size
plt.hist((sasa_arr2_hist), density=True, bins=30, weights=weights2)
plt.xlabel("SASA Ensemble Ratio", fontsize=16)
plt.ylabel("Probability Density", fontsize=16)
plt.xlim(80, 140)

    
toc = time.perf_counter()
print(f"Program ran in {toc-tic:0.4f} seconds")