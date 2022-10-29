# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 02:58:42 2022

@author: Thomas Tsangaris
"""

"""Purpose: Compute weighted FRET efficiency from pre-calculated files. The 
calculations are typically performed from an accessible volume simulation
(see the avtraj Python package).

"""
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


# ensemble 1:
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_FFT_73_121_BME_calc_FRET_avtraj.dat'
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/5P_FFT_73_121_BME_calc_FRET_avtraj.dat'

# mypath = 'D:/NP_FFT_IP_2/NP_FFT_73_121_BME_calc_FRET_avtraj.dat'
# mypath = 'D:/5P_FFT_IP_2/5P_FFT_73_121_BME_calc_FRET_avtraj.dat'

weights_file = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/NP_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_65_omega_50.dat'
# weights_file = 'C:/Users/14169/Dropbox/4E-BP2 projects/BME/Weights/5P_2FRET_SAXS_CS_PRE_intensity_validation/weights_theta_25_omega_1.dat'

# weights_file = 'D:/NP_FFT_IP_2/NP_omega_results_2FRET_SAXS_CS_19_7_22/Weights/weights_theta_65_omega_50.dat'
# weights_file = 'D:/5P_FFT_IP_2/5P_omega_results_2FRET_SAXS_CS_19_7_22/Weights/weights_theta_50_omega_1.dat'
weights = get_weights(weights_file)

E_vals = []
hl = 1

with open(mypath, 'r') as f:                                                                                          
    reader = csv.reader(f, delimiter='\t')
    for j in range(hl):
        # skip header
        next(reader)
    for row in reader:
        #Read remaining lines (row at a time)
        E_vals.append(float(row[1]))
    f.close()
    
print(np.average(E_vals, weights=weights))