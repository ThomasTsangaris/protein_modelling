# author: Thomas Tsangaris

# Description: This script creates a new folder filled with the same pdbs you supply
# but any atom registered as "TPO" or "SEP" is rewritten as "THR" or "SER" 
# respectively. To be clear, this script simply renames TPOs to THRs and 
# SEPs to SERs

import glob
import os
import time

ti=time.perf_counter()

# Path to your ensemble (folder of pdb files):
# path = '/mnt/e/5P4E-BP2_FFT/FFT_1000/'

# path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'
# path = 'D:/5p_fft_theta_25_omega_1_desktop_transfer/'
# path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/bme_pdb_ensemble_5P_2FRET_SAXS_CS_PRE_intensity_validation_25_omega_1/' # 5P  

# 5P 2 (N=5000)for confidence:
# path = 'D:/5P_FFT_IP_2/'
# 5P 3
# path = 'D:/5P_FFT_IP_3/'
# 5P 4
path = 'D:/5P_FFT_IP_4/'


pdb_files = glob.glob(path + '*.pdb') # list of pdb files

# new_dir_name = 'HYCUD_PDBs/' # directory where the results will be placed
new_dir_name = 'removed_phosphates/' # directory where the results will be placed

if not os.path.exists(path + new_dir_name): # if the results folder does not exist
    os.makedirs(path + new_dir_name) # create the results folder
os.chdir(path + new_dir_name) # change current working directory 

j = 0
    
for pdb in pdb_files:

    new_pdbs_lst = []

    lines = open(pdb, 'r').readlines()
    for i in range(len(lines)):
        if lines[i][17:20] == 'SEP':
             lines[i] = lines[i][:17] + 'SER' + lines[i][20:]
        elif lines[i][17:20] == 'TPO':
             lines[i] = lines[i][:17] + 'THR' + lines[i][20:]
    new_pdbs_lst.append(lines)
    
    for new_pdb_lines in new_pdbs_lst:
        with open(path + new_dir_name + str(j) + '_edited.pdb', 'w') as new_pdb:
            new_pdb.writelines(new_pdb_lines)
        j += 1
        

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")