# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 20:59:15 2022

@author: Thomas Tsangaris
"""

""" Purpose: Input an ensemble of pdb files, output normalized weights based on\
    occurence of the pdb in the ensemble. For example, a protein ensemble of 4
    conformations in which the first 3 are the same structure and the fourth is
    different, would result in weights of 0.75, 0.25; a weight for every 
    distinct conformer.
    
    In the case of the inputted ensemble being sampled from specific weights,
    the original ensemble can be inputted along with the original weights and 
    the weights of both ensembles will be compared via a residuals plot.

"""
import matplotlib.pyplot as plt
import numpy as np
import time
import glob
import re
import csv
import time
import os



tic = time.perf_counter()


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


# input path to ensembles
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50/' # NP
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/' # 5P
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1_10000/' # 5P
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1_15000/' # 5P

# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/bme_pdb_ensemble_NP_2FRET_SAXS_CS_PRE_intensity_validation_theta_65_omega_50_20000/'

# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'
# mypath1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/'
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


# Bioclusters, Clustering on Initial Pool:
# NP FFT Cluster 1
# mypath1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_1/"
# NP FFT Cluster 2
# mypath1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_2/"
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



# testing convergence
mypath2 = "D:/testing_weight_conv/" # testing sampling convergence
mypath1 = mypath2 + '800_1_conf/' # sampled with weights
weights_path = mypath2 + "test_50.dat" # testing sampling convergence



# Bioclusters, Clustering on Initial Pool:
# NP FFT Cluster 1
# weights_path = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_1/NP_FFT_cluster_1_weights.dat"
# NP FFT Cluster 2
# weights_path = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_2/NP_FFT_cluster_2_weights.dat"
# NP FFT Cluster 3
# weights_path = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_3/NP_FFT_cluster_3_weights.dat"
# NP FFT Cluster 4
# weights_path = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_4/NP_FFT_cluster_4_weights.dat"
    
# 5P cluster 1
# weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_1/5P_FFT_cluster_1_weights.dat'
# 5P cluster 2
# weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_2/5P_FFT_cluster_2_weights.dat'
# 5P cluster 3
# weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_3/5P_FFT_cluster_3_weights.dat'


pdb_files1 = [os.path.basename(pdb) for pdb in glob.glob(mypath1 + '*.pdb')]
N1 = len(pdb_files1)
print(N1)

# weights paths
# weights_path = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_path = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'

weights = get_weights(weights_path)

regex = '(\d+)_' # depends on how your pdbs are named
pdb_num_to_count = {}

for pdb in pdb_files1:
    pdb_num = int(re.search(regex, pdb).group()[:-1])
    if pdb_num not in pdb_num_to_count:
        pdb_num_to_count[pdb_num] = 0
    pdb_num_to_count[pdb_num] += 1

weight_residuals = [] # conformer weight differences
weight_quotients = [] # conformer weight quotients
weights1 = [] # weights for first inputted ensemble
conf_nums = range(len(weights))

for pdb_num in range(1, len(weights)+1):
    if pdb_num not in pdb_num_to_count: # conformer is not in the first ensemble
        pdb_num_to_count[pdb_num] = 0
        
    weight1 = pdb_num_to_count[pdb_num] / N1
    weight2 = weights[pdb_num-1]
    print(pdb_num)
    weight_residual = (weight1) - weight2
    weight_residuals.append(weight_residual)
    
    weights1.append(weight1)
    
    if (weight1==0) or (weight2==0):
        weight_quotients.append(0)
    elif weight1 > weight2:
        weight_quotients.append(weight1/weight2)
    else:
        weight_quotients.append(weight2/weight1)


# for equal/original weights:
equal_weights = []
for i in range(len(weights)):
    equal_weights.append(1/len(weights))
equal_weights = np.array(equal_weights)




# compute the relative shannon entropy
s_rel = 0
i = 0
for pdb_num in pdb_num_to_count:
    print(pdb_num_to_count)
    weight1 = pdb_num_to_count[pdb_num] / N1
    equal_weight = equal_weights[i]
    # print(weight1)
    # print(weight1==0)
    if weight1 != 0: # reference justification in this link: https://stats.stackexchange.com/questions/57069/alternative-to-shannons-entropy-when-probability-equal-to-zero
        s_rel += weight1 * np.log(weight1/equal_weight)
    # print(weight1)
    i += 1

s_rel *= -1
Neff = np.exp(s_rel)
# print the Neff value
print("Neff for first inputted ensemble is:", Neff) # first is sampled

# compute the relative shannon entropy
s_rel2 = 0
i2 = 0
for pdb_num in conf_nums:
    weight2 = weights[pdb_num]
    equal_weight = equal_weights[i2]
    # print(weight2)
    # print(weight2==0)
    if weight2 != 0: # reference justification in this link: https://stats.stackexchange.com/questions/57069/alternative-to-shannons-entropy-when-probability-equal-to-zero
        s_rel2 += weight2 * np.log(weight2/equal_weight)
    # print(weight2)
    i2 += 1
s_rel2 *= -1
Neff2 = np.exp(s_rel2)
# print the Neff value
print("Neff for second inputted ensemble is:", Neff2) # first is sampled


plt.figure(figsize=(9,4)) # adjust figure size
plt.plot(conf_nums, weight_residuals)
plt.xlabel("Conformer Number")
plt.ylabel("Weight Difference")

plt.figure(figsize=(9,4)) # adjust figure size
plt.plot(conf_nums, weights1)
# plt.plot(conf_nums, weights, linestyle='--')
plt.xlabel("Conformer Number")
plt.ylabel("Weights")

plt.figure(figsize=(9,4)) # adjust figure size
plt.plot(conf_nums, weights, linestyle='--')
plt.xlabel("Conformer Number")
plt.ylabel("Weights")

print(np.mean(weight_quotients))

# plt.axhline(max(weights)*0.1)
# plt.axhline(-max(weights)*0.1)

# get values for plot:
pdb_count_lst = []
for pdb_num in pdb_num_to_count:
    for count in range(pdb_num_to_count[pdb_num]):
        pdb_count_lst.append(pdb_num)
        
plt.figure(figsize=(9,4)) # adjust figure size
plt.bar(range(1, len(weights)+1), weights1, color='r', alpha=0.5) # sampled ensemble
plt.bar(range(1, len(weights)+1), np.array(weights)/sum(weights), color='b', alpha=0.5) # exact
plt.xlabel("Conformer Number", fontsize=16)
plt.ylabel("Weight", fontsize=16)

toc = time.perf_counter()
print(f"Program ran in {toc-tic:0.4f} seconds")