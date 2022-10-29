# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 18:57:36 2022

@author: Thomas Tsangaris
"""
# Purpose: Parse and perform a weighted average of HYCUD results. Weights are
# assumed to be in a BME output format.
import csv
import numpy as np
from scipy import stats


def get_weights(weights_file: str) -> np.array:
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

weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path2 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/ss_region_ensemble/ss_region_weights.dat'
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/no_ss_region_weights.dat'

# NP cluster 1
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_1/NP_FFT_cluster_1_weights.dat'
# NP cluster 2
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_2/NP_FFT_cluster_2_weights.dat'
# NP cluster 3
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/cluster_3/NP_FFT_cluster_3_weights.dat'

# 5P cluster 1
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_1/5P_FFT_cluster_1_weights.dat'
# 5P cluster 2
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_2/5P_FFT_cluster_2_weights.dat'
# 5P cluster 3
# weights_path1 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/cluster_3/5P_FFT_cluster_3_weights.dat'

# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/'
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_32_91/'
# mypath1 = 'C:/Users/14169/Dropbox/Local BP2/np_rc_73_121/np_rc_73_121/'
# mypath2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/'
# ensemble2 = glob.glob(mypath2 + "*.pdb") # list of pdb file paths
# weights_path1 = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path2 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/no_ss_region_weights.dat'



# Bioclusters, Clustering on Initial Pool:
# NP FFT Cluster 1
# weights_path1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_1/NP_FFT_cluster_1_weights.dat"
# NP FFT Cluster 2
# weights_path1 = "C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/ip_bio_cluster_2/NP_FFT_cluster_2_weights.dat"
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


# for equal weights:
# numfiles = len(ensemble1) # for first ensemble, they may have different Nconf
# equal_weights1 = []
# for i in range(numfiles):
#     equal_weights1.append(1/numfiles)
# equal_weights1 = np.array(equal_weights1)


weights = get_weights(weights_path1)


hycud_results = 'D:/hycud_transfer/NP_FFT_IP_5000.csv'
hl = 2 # 2 header lines to skip
res_num_counter = 0
eta_dict = {}
r_corr_dict = {}
hm_corr_dict = {}
r_dict = {}
hm_dict = {}
with open(hycud_results, 'r') as f:                                                                                          
    reader = csv.reader(f, delimiter=',')
    for j in range(hl):
        # skip header
        next(reader)
    for row in reader:
        #Read remaining lines (row at a time)
        frag_num = float(row[1]) # all of these should be the same for the same residue span
        res_num_counter += 1
        res_num = res_num_counter # residue number
        eta = float(row[5])
        r_corr = float(row[6])
        hm_corr = float(row[7])
        r = float(row[8])
        hm = float(row[9])
        
        if frag_num not in eta_dict:
            eta_dict[frag_num] = []
        if frag_num not in r_corr_dict:
            r_corr_dict[frag_num] = []
        if frag_num not in hm_corr_dict:
            hm_corr_dict[frag_num] = []
        if frag_num not in r_dict:
            r_dict[frag_num] = []
        if frag_num not in hm_dict:
            hm_dict[frag_num] = []
        
        eta_dict[frag_num].append(eta)
        r_corr_dict[frag_num].append(r_corr)
        hm_corr_dict[frag_num].append(hm_corr)
        r_dict[frag_num].append(r)
        hm_dict[frag_num].append(hm)
        
    f.close()


for frag_num in hm_corr_dict:
    
    # Computing weighted Means:
    # just run without weights argument for unweighted
    # print("Average hm_corr of Fragment", str(frag_num)+':', str(np.average(np.array(hm_corr_dict[frag_num]), weights=weights, axis=0)))
    # print("Average hm of Fragment", str(frag_num)+':', str(np.average(np.array(hm_dict[frag_num]), weights=weights, axis=0)))
    print("Average r_corr of Fragment", str(frag_num)+':', str(np.average(np.array(r_corr_dict[frag_num]), weights=weights, axis=0)))
    # print("Average r of Fragment", str(frag_num)+':', str(np.average(np.array(r_dict[frag_num]), weights=weights, axis=0)))
    # print("Average eta of Fragment", str(frag_num)+':', str(np.average(np.array(eta_dict[frag_num]), weights=weights, axis=0)))


    # Computing weighted SDOM:
    # weighted_mean_hm_corr = np.average(np.array(hm_corr_dict[frag_num]), weights = weights, axis=0)
    # weighted_variance_hm_corr = np.average((np.array(hm_corr_dict[frag_num]) - weighted_mean_hm_corr)**2, weights=weights)
    # weighted_std_hm_corr = np.sqrt(weighted_variance_hm_corr)
    # weighted_sdom_hm_corr = weighted_std_hm_corr / np.sqrt(len(hm_corr_dict[frag_num]))
    # print("STDOM hm_corr of Fragment", str(frag_num)+':', str(weighted_sdom_hm_corr))
    
    # weighted_mean_hm = np.average(np.array(hm_dict[frag_num]), weights = weights, axis=0)
    # weighted_variance_hm = np.average((np.array(hm_dict[frag_num]) - weighted_mean_hm)**2, weights=weights)
    # weighted_std_hm = np.sqrt(weighted_variance_hm)
    # weighted_sdom_hm = weighted_std_hm / np.sqrt(len(hm_dict[frag_num]))
    # print("STDOM hm of Fragment", str(frag_num)+':', str(weighted_sdom_hm))
    
    weighted_mean_r_corr = np.average(np.array(r_corr_dict[frag_num]), weights = weights, axis=0)
    weighted_variance_r_corr = np.average((np.array(r_corr_dict[frag_num]) - weighted_mean_r_corr)**2, weights=weights)
    weighted_std_r_corr = np.sqrt(weighted_variance_r_corr)
    weighted_sdom_r_corr = weighted_std_r_corr / np.sqrt(len(r_corr_dict[frag_num]))
    print("STDOM r_corr of Fragment", str(frag_num)+':', str(weighted_sdom_r_corr))

    # weighted_mean_r = np.average(np.array(r_dict[frag_num]), weights = weights, axis=0)
    # weighted_variance_r = np.average((np.array(r_dict[frag_num]) - weighted_mean_r)**2, weights=weights)
    # weighted_std_r = np.sqrt(weighted_variance_r)
    # weighted_sdom_r = weighted_std_r / np.sqrt(len(r_dict[frag_num]))
    # print("STDOM r of Fragment", str(frag_num)+':', str(weighted_sdom_r))
    
    # weighted_mean_eta = np.average(np.array(eta_dict[frag_num]), weights = weights, axis=0)
    # weighted_variance_eta = np.average((np.array(eta_dict[frag_num]) - weighted_mean_eta)**2, weights=weights)
    # weighted_std_eta = np.sqrt(weighted_variance_eta)
    # weighted_sdom_eta = weighted_std_eta / np.sqrt(len(r_dict[frag_num]))
    # print("STDOM eta of Fragment", str(frag_num)+':', str(weighted_sdom_eta))