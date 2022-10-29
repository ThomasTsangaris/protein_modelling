# -*- coding: utf-8 -*-
"""
Created on Sun Aug  1 22:42:29 2021

@author: Thomas Tsangaris
"""

# Purpose: To calculate the Hydrodynamic Radius of a protein ensemble given 
# the path to its HYDROPRO results produced using the script batch_hydropro.py
# The user can also supply a weights file (in BME format) and calculate the 
# weighted value of their ensemble

import glob
import sys
import subprocess
import os
import shutil
import math
import statistics
import time
import csv
import numpy as np
import matplotlib.pyplot as plt

R_h_lst = []
results_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/' + 'HYDROPRO_Results/'
# results_path = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/' + 'HYDROPRO_Results/'
# results_path = '/mnt/c/NP_FFT_IP/' + 'HYDROPRO_Results/'
# results_path = 'C:/NP_FFT_IP/' + 'HYDROPRO_Results/' # directory where the results will be placed
os.chdir(results_path) # change current working directory 

txt_files = glob.glob(results_path + '*.txt') # list of all text files
res_files = [] # list of all res files



for txt in txt_files: # go through all text files
    if txt[-8:-4] == '-res': # is a res file
        res_files.append(txt) # add to list of res files
i = 1
for res in res_files: # loop through each res file (for each conformer)
    # if i < 501:
    # print(i)
    res_file = open(res[len(results_path):], 'r')
    res_lines = res_file.readlines() # list of lines in res file
    T = float(res_lines[18][41:46]) + 273.15 # Temperature in K
    # print("T: " + str(T)) #debug
    eta = float(res_lines[19][41:49]) / 10 # Solvent Viscosity in N*s/m^2
    # print("eta: " + str(eta)) #debug
    D_T = float(res_lines[24][41:50]) / 10000 # Translational Diffusion Coefficient in m^2/s
    # print("D_T: " + str(D_T)) #debug
    D_R = float(res_lines[27][41:50]) # Rotational Diffusion Coefficient in 1/s
    # print("D_R: " + str(D_R)) #debug
    k_B = 1.38065e-23 # Boltazmann Constant in m^2*kg/(s^2*K)
    
    R_h = (k_B * T) / (6 * math.pi * eta * D_T) # Hydrodynamic Radius in m
    # R_h = (k_B * T) / (6 * math.pi * eta * D_R) # Hydrodynamic Radius in m debug purposes
    R_h *= 1e10 # Hydrodynamic Radius in Angstroms
    
    R_h_lst.append(R_h) # add pdb's R_h value to the

    i += 1
    
if (len(R_h_lst) >= 1):
    mean = statistics.mean(R_h_lst) # Mean R_h in Angstroms 
if len(R_h_lst) > 1:
    std = statistics.stdev(R_h_lst) # Standard Deviation for R_h in Angstroms
    std_of_mean = std / math.sqrt(len(R_h_lst)) # Standard Deviation of the Mean for R_h in Angstroms
if (len(R_h_lst) >= 1):
    print('Mean Rh: ' + str(mean) + '\n')
if len(R_h_lst) > 1:
    print('STD Rh: ' + str(std) + '\n')
    print('SDOM Rh: ' + str(std_of_mean) + '\n')
    
print("After Reweighting:")
# Input Weights File Here:
weights_file = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_file = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'
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
print("Sum of Weights Should be One:", str(sum(weights)))

if len(R_h_lst) > 1:
    # Weighted Average:
    weighted_mean = np.average(R_h_lst, weights = weights)
    weighted_variance = np.average((R_h_lst - weighted_mean)**2, weights=weights)
    weighted_std = np.sqrt(weighted_variance)
    weighted_sdom = weighted_std / np.sqrt(len(R_h_lst))
    print('Mean Rh: ' + str(weighted_mean) + '\n')
    print('STD Rh: ' + str(weighted_std) + '\n')
    print('SDOM Rh: ' + str(weighted_sdom) + '\n')
    
    fig, ax = plt.subplots()
    ax.hist(R_h_lst, 25, weights=weights, density = True, edgecolor = 'black')
    # y, x = np.histogram(all_E, bins=25, weights=weights, density=True)
    # ax.plot(x[:-1], y)
    # fig.show()
    
    # print("x = ", x)
    # print("y = ", y)
    # plt.figure()
    # sns.histplot(x=x[:-1], y=y, bins = len(y), stat = 'density', color = 'b', kde = False)
    plt.xlabel(r"$R_h$")
    plt.ylabel("Probability Density")

print("DONE")