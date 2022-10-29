# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:03:35 2022

@author: Thomas Tsangaris
"""

"""Purpose: compute the reduced chi^2 in the same way that BME does from 
the inputs that BME takes (experimental and back-calculated values in BME
format). ALso, you can input BME weights where each row contains a weight for 
a conformer and the file has one header line that is not read. This script can
also compute the root-mean-square deviation instead of reduced chi^2, but only
for PRE intensity data where the back-calculated file and the experimental file has each row representing a datapoint, the 
first two columns specifying the two labelling positions, the third the PRE 
intensity ratio and the associated uncertainty as the fourth column. The 
experimental PRE file must have no header line.

Note: BME = Bayesian Maximum Entropy, a method developed by Lindorff-Larsen's
group.

"""
from typing import Tuple
import csv
import numpy as np


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


def parse_exp_bme(exp_bme: str) -> Tuple[np.array, np.array]:
    """Return an array containing the experimental values and uncertainties
    given a path to an experimental file in BME format.
    
    """
    hl = 1
    exp_vals = []
    exp_uncs = []
    with open(exp_bme, 'r') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for j in range(hl):
            # skip header
            next(reader)
        for row in reader: # Read remaining lines (row at a time)
            exp_vals.append(float(row[1]))
            exp_uncs.append(float(row[2]))
        f.close()
    return np.array(exp_vals), np.array(exp_uncs)
    

def parse_bc_bme(bc_bme: str) -> np.array:
    """Return a tuple containing the back-calculated values given a path to an 
    back-calculated file in BME format.
    
    """
    hl = 1
    bc_vals = []
    with open(bc_bme, 'r') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for j in range(hl):
            # skip header
            next(reader)
        for row in reader: # Read remaining lines (row at a time)
            bc_vals.append(row[1:])
        f.close()
    return np.array(bc_vals)


def parse_pre_exp_bme(exp_bme: str) -> Tuple[np.array, np.array]:
    """Return an array containing the experimental values and uncertainties
    given a path to a PRE experimental file that has each row 
    representing a datapoint, the first two columns specifying the two 
    labelling positions, the third the PRE intensity ratio and the associated
    uncertainty as the fourth column. The PRE experimental file must have no
    header line.
    
    """
    exp_vals = []
    exp_uncs = []
    with open(exp_bme, 'r') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for row in reader: # Read remaining lines (row at a time)
            exp_vals.append(float(row[2]))
            exp_uncs.append(float(row[3]))
        f.close()
    return np.array(exp_vals), np.array(exp_uncs)


def parse_pre_bc_bme(bc_bme: str) -> Tuple[np.array, np.array]:
    """Return an array containing the back-calculated values
    given a path to a PRE back-calculated file that has each row 
    representing a datapoint, the first two columns specifying the two 
    labelling positions and the third being the mean PRE intensity ratio. 
    The PRE back-calculated file must have no header line.
    
    """
    bc_vals = []
    with open(bc_bme, 'r') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for row in reader: # Read remaining lines (row at a time)
            bc_vals.append(float(row[2]))
        f.close()
    return np.array(bc_vals)


def compute_red_chi2(exp_bme: str, bc_bme: str, weights: list = []) -> float:
    """Return the reduced chi squared value given valid paths to experimental 
    and back-calculated files in BME format.
    
    """
    exp_vals, exp_uncs = parse_exp_bme(exp_bme)
    bc_vals = parse_bc_bme(bc_bme)
    sum_ = []
    
    if weights == []: # no weights supplied 
        for i in range(np.shape(bc_vals)[0]):
            weights.append(1/np.shape(bc_vals)[0]) # use equal weights (1/Nconf)
    
    for i in range(np.shape(bc_vals)[1]): # loop through each column, number of datapoints
        bc_mean_val = np.average(bc_vals[:,i].astype(np.float), weights=weights) # average over all confs for bc value
        exp_val = float(exp_vals[i])
        exp_unc = float(exp_uncs[i])
        sum_.append(((bc_mean_val - exp_val)**2) / (exp_unc**2))
    
    red_chi2 = sum(sum_) / len(exp_vals)
    return red_chi2


def compute_pre_rmsd(exp_bme: str, bc_bme: str) -> float:
    """Return the root-mean-square deviation value given valid paths to 
    experimental and back-calculated files where the back-calculated 
    file and the experimental file has each row 
    representing a datapoint, the first two columns specifying the two 
    labelling positions, the third the PRE intensity ratio and the associated
    uncertainty as the fourth column. The experimental PRE file must have no
    header line.
    
    """
    exp_vals, exp_uncs = parse_pre_exp_bme(exp_bme)
    bc_vals = parse_pre_bc_bme(bc_bme)
    squared_dif = []
    
    for i in range(np.shape(bc_vals)[0]): # loop through each column, number of datapoints
        bc_mean_val = float(bc_vals[i]) # average over all confs for bc value
        exp_val = float(exp_vals[i])
        exp_unc = float(exp_uncs[i])  # unused since rmsd does not use uncertainty
        squared_dif.append(((bc_mean_val - exp_val)**2))
    
    rmsd = np.sqrt(np.nanmean(squared_dif))
    return rmsd


# Paths:
# NP *
# exp_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP4E-BP2_SAXS_Ang_bme.dat'
# calc_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_BME_SAXS_PEPSI_dro018.dat' # has mutant cys
# # calc_SAXS = 'D:/np_fft_theta_65_omega_50_desktop_transfer/NP_FFT_BME_sampled_SAXS_PEPSI_dro018.dat' # BME sampled
# # calc_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/np_crysol_saxs_bc_bme.dat'
# exp_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/BME_exp_CS_CACB.dat'
# # calc_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/BME_calc_CS_CACB.dat'
# calc_CS = 'D:/np_fft_theta_65_omega_50_desktop_transfer/BME_calc_CS_CACB.dat' # BME sampled
# # 32-91 and 73-121:
# # exp_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/FRET3291_FRET73121_BME_exp_FRET_Omega_50.dat'
# # calc_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_calc_32_91_73_121_conjoined.dat'
# # 32-91 only:
# exp_32_91 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_32_91_FRET_omegas/FRET3291_BME_exp_FRET_Omega_50.dat'
# # calc_32_91 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_32_91_BME_calc_FRET_avtraj.dat'
# calc_32_91 = 'D:/np_fft_theta_65_omega_50_desktop_transfer/NP_FFT_32_91_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
# # 73-121 only:
# exp_73_121 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_73_121_FRET_omegas/FRET73121_BME_exp_FRET_Omega_50.dat'
# calc_73_121 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_73_121_BME_calc_FRET_avtraj.dat'
# # calc_73_121 = 'D:/np_fft_theta_65_omega_50_desktop_transfer/NP_FFT_73_121_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
# exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_PRE_Intensity_Height_Ratios.dat' # in intensity height ratios 
# # calc_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/NP_FFT_5000_bc.dat'
# calc_PRE = 'D:/np_fft_theta_65_omega_50_desktop_transfer/NP_FFT_5000_bc.dat' # BME sampled

# 32-91 only:
# exp_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_32_91_FRET_omegas/FRET3291_BME_exp_FRET_Omega_' \
#     + str(Omega) + '.dat'
# calc_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_32_91_BME_calc_FRET_avtraj.dat'
# # 73-121 only:
# exp_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_73_121_FRET_omegas/FRET73121_BME_exp_FRET_Omega_' \
#     + str(Omega) + '.dat'
# calc_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_73_121_BME_calc_FRET_avtraj.dat'

#  # # #PRE data
#  # # # exp_PRE = './Experimental/Sic1_PRE_rates.dat'
#  # # exp_PRE = '/mnt/c/NP_FFT_IP/NP_PRE_exp_rates.dat'

#     calc_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/NP_calc_PRE_dis_bme.dat' # uses CB as label
# #PRE data
# # exp_PRE = './Experimental/Sic1_PRE_rates.dat'
# exp_PRE = '/mnt/c/NP_FFT_IP/NP_PRE_exp_rates.dat'

# 5P *
# exp_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/p4E_SAXS_Ang_bme.dat'
# # # calc_SAXS = mypath + '5P_FFT_BME_SAXS_PEPSI_dro018.dat' # has mutant cys
# calc_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_FFT_BME_SAXS_PEPSI_dro018.dat'
# # calc_SAXS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_BME_sampled_SAXS_PEPSI_dro018.dat' # BME sampled
# exp_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_CS_no_comments_bme_no_18_62.dat' # used in BME optimization
# # exp_CS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/BME_exp_CS_CACB.dat' # full thing
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/' # path is only used to grab experimental files
# calc_CS = mypath + 'BME_calc_CS_CACB_no_18_62.dat'
# # calc_CS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/BME_calc_CS_CACB.dat' # BME sampled
# # # exp_FRET = mypath + '5P_32_91_73_121_FRET_omegas/FRET3291_FRET73121_BME_exp_FRET_Omega_1.dat'
# # # calc_FRET = mypath + '5P_FFT_calc_32_91_73_121_conjoined.dat'
# exp_32_91 = mypath + '5P_32_91_FRET_omegas/FRET3291_BME_exp_FRET_Omega_1.dat'
# calc_32_91 = mypath + '5P_FFT_32_91_BME_calc_FRET_avtraj.dat' 
# # calc_32_91 = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_32_91_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
# exp_73_121 = mypath + '5P_73_121_FRET_omegas/FRET73121_BME_exp_FRET_Omega_1.dat'
# calc_73_121 = mypath + '5P_FFT_73_121_BME_calc_FRET_avtraj.dat'
# # calc_73_121 = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_73_121_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
# # # # exp_FRET = exp_32_91 # only restrain with 32-91
# # # # calc_FRET = calc_32_91 # only restrain with 32-91
# # # # exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_PRE_Intensity_Height_Ratios.dat'
# exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_PRE_Intensity_Height_Ratios_pre_removed_folded.dat'
# calc_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/removed_phosphates/5P_FFT_5000_bc.dat'
# # # exp_PRE_dis = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/PREdist_wt_p4EBP2_bme.dat'
# # calc_PRE = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_5000_bc.dat' # BME sampled

# weights_path = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_path = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
# weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/ss_region_ensemble/ss_region_weights.dat'
# weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/no_ss_region_weights.dat'



# Bioclusters, Clustering on Initial Pool:
# NP:
# shortcut_path = 'D:/NP_clusters/ip_bio_cluster_1/'
# shortcut_path = 'D:/NP_clusters/ip_bio_cluster_2/'
# shortcut_path = 'D:/NP_clusters/ip_bio_cluster_3/'
# shortcut_path = 'D:/NP_clusters/ip_bio_cluster_4/'

# 5P:
# shortcut_path = 'D:/5P_clusters/ip_bio_cluster_1/'
# shortcut_path = 'D:/5P_clusters/ip_bio_cluster_2/'
shortcut_path = 'D:/5P_clusters/ip_bio_cluster_3/'

# NP
# exp_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP4E-BP2_SAXS_Ang_bme.dat'
# # calc_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_BME_SAXS_PEPSI_dro018.dat' # has mutant cys
# calc_SAXS = shortcut_path + 'NP_FFT_BME_SAXS_PEPSI_dro018.dat' # BME sampled
# # calc_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/np_crysol_saxs_bc_bme.dat'
# exp_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/BME_exp_CS_CACB.dat'
# # calc_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/BME_calc_CS_CACB.dat'
# calc_CS = shortcut_path + 'BME_calc_CS_CACB.dat' # BME sampled
# # 32-91 and 73-121:
# # exp_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/FRET3291_FRET73121_BME_exp_FRET_Omega_50.dat'
# # calc_FRET = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_calc_32_91_73_121_conjoined.dat'
# # 32-91 only:
# exp_32_91 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_32_91_FRET_omegas/FRET3291_BME_exp_FRET_Omega_50.dat'
# # calc_32_91 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_32_91_BME_calc_FRET_avtraj.dat'
# calc_32_91 = shortcut_path + 'NP_FFT_32_91_BME_calc_FRET_avtraj.dat' # BME sampled
# # 73-121 only:
# exp_73_121 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_73_121_FRET_omegas/FRET73121_BME_exp_FRET_Omega_50.dat'
# # calc_73_121 = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_73_121_BME_calc_FRET_avtraj.dat'
# calc_73_121 = shortcut_path + 'NP_FFT_73_121_BME_calc_FRET_avtraj.dat' # BME sampled
# exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_PRE_Intensity_Height_Ratios.dat' # in intensity height ratios
# # #calc_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/NP_FFT_5000_bc.dat'
# # calc_PRE = shortcut_path + 'NP_FFT_5000_bc.dat' # BME sampled # not ready yet


# 5P
exp_SAXS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/p4E_SAXS_Ang_bme.dat'
# # calc_SAXS = mypath + '5P_FFT_BME_SAXS_PEPSI_dro018.dat' # has mutant cys
calc_SAXS = shortcut_path + '5P_FFT_BME_SAXS_PEPSI_dro018.dat'
# calc_SAXS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_BME_sampled_SAXS_PEPSI_dro018.dat' # BME sampled
exp_CS = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_CS_no_comments_bme_no_18_62.dat' # used in BME optimization
# exp_CS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/BME_exp_CS_CACB.dat' # full thing
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/' # path is only used to grab experimental files
calc_CS = shortcut_path + 'BME_calc_CS_CACB.dat'
# calc_CS = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/BME_calc_CS_CACB.dat' # BME sampled
# # exp_FRET = mypath + '5P_32_91_73_121_FRET_omegas/FRET3291_FRET73121_BME_exp_FRET_Omega_1.dat'
# # calc_FRET = mypath + '5P_FFT_calc_32_91_73_121_conjoined.dat'
exp_32_91 = mypath + '5P_32_91_FRET_omegas/FRET3291_BME_exp_FRET_Omega_1.dat'
calc_32_91 = shortcut_path + '5P_FFT_32_91_BME_calc_FRET_avtraj.dat' 
# calc_32_91 = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_32_91_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
exp_73_121 = mypath + '5P_73_121_FRET_omegas/FRET73121_BME_exp_FRET_Omega_1.dat'
calc_73_121 = shortcut_path + '5P_FFT_73_121_BME_calc_FRET_avtraj.dat'
# calc_73_121 = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_73_121_BME_sampled_calc_FRET_avtraj.dat' # BME sampled
# # # exp_FRET = exp_32_91 # only restrain with 32-91
# # # calc_FRET = calc_32_91 # only restrain with 32-91
# # # exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_PRE_Intensity_Height_Ratios.dat'
exp_PRE = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_PRE_Intensity_Height_Ratios_pre_removed_folded.dat'
calc_PRE = shortcut_path + 'removed_phosphates/5P_FFT_5000_bc.dat' # not created yet
# # exp_PRE_dis = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/PREdist_wt_p4EBP2_bme.dat'
# calc_PRE = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/5P_FFT_5000_bc.dat' # BME sampled



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
weights_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/ip_bio_cluster_3/5P_FFT_cluster_3_weights.dat'
    

weights = get_weights(weights_path)
print(compute_red_chi2(exp_SAXS, calc_SAXS, weights))
print(compute_red_chi2(exp_CS, calc_CS, weights))
print(compute_red_chi2(exp_32_91, calc_32_91, weights))
print(compute_red_chi2(exp_73_121, calc_73_121, weights))
print(compute_pre_rmsd(exp_PRE, calc_PRE, weights))

# print(compute_red_chi2(exp_SAXS, calc_SAXS))
# print(compute_red_chi2(exp_CS, calc_CS))
# print(compute_red_chi2(exp_32_91, calc_32_91))
# print(compute_red_chi2(exp_73_121, calc_73_121))
# print(compute_pre_rmsd(exp_PRE, calc_PRE))