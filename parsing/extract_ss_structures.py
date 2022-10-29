# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 19:36:37 2022

@author: Thomas Tsangaris
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import mdtraj
import time
import shutil
import os
import stat

ti=time.perf_counter()

""" Purpose: extract structures that contain a certain secondary structural (ss) 
element (user inputs pdb ensemble, type of ss element and residue range of ss element).
The DSSP criteria is used to determine the ss of a residue range
"""

# input ensemble
# NP (N=5000):
mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/FFT_NP_Initial_Pool_cys/'

# 5P (N=5000):
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/5P_FFT_IP_5000/'

# TraDES (73 and 121 are cysteines) (N=5000):
# mypath = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/'
    
# IDP Conf Gen (ICG):
# mypath = ''

# define mdtraj trajectory
pdb_files = glob.glob(mypath+"*.pdb")
traj = mdtraj.load(pdb_files, top=pdb_files[0])

# type of ss element: H=helix, E=extended, C=coil (everything other than H and E)
ss_type = 'H'

# lower bound of residue range for ss element
lower = 54
# lower = 1

# higher bound of residue range for ss element
upper = 63
# upper = 3

# minimum number of residues in range to be classified as the ss region
threshold = 5

def extract_ss_segment(traj, ss_type, lower, upper, threshold = upper-lower) -> tuple: 
    """Return an ensemble which contains the secondary structural element
    defined by the secondary structure type ss_type and range of residues 
    participating in that structure from lower to upper and also return an 
    ensemble which does not contain the specified ss region A threshold can also 
    be defined, in that you can specify a number of residues in the range 
    that (inclusively) matches the criteria of having that ss region.
    
    Precondition: threshold <= upper - lower
    
    """
    # Back-calculate Secondary Structure:
    # array containing pdb structures that have the designated ss region 
    ss_arr = []
    # array containing pdb structures that do not have the designated ss region
    no_ss_arr = []
    
    i = 0
    for frame in traj: # loop through each pdb file in the ensemble
        # determine secondary structure code of specified range
        dssp_code = mdtraj.compute_dssp(frame)[0][lower-1:upper] # parse for just protein/peptide
        # number of helical residues in the specified range
        num_helix = list(dssp_code).count('H')
        
        if num_helix >= threshold: # if the ss region is present in the conformer
            ss_arr.append(pdb_files[i])
        else: # otherwise it does not have the region
            no_ss_arr.append(pdb_files[i])
        i += 1
    
    return ss_arr, no_ss_arr

ss_arr, no_ss_arr \
    = extract_ss_segment(traj, ss_type, lower, upper, threshold)

print("There are", len(ss_arr), "conformers with the designated ss region")
print("There are", len(no_ss_arr), "conformers without the designated ss region")

# make folders and store pdb ensembles in them
# output folder for conformations that have ss region
ss_output_folder =  'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/ss_region_ensemble/'
# output folder for conformations that do not have ss region
no_ss_output_folder = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/np_rc_73_121/no_ss_region_ensemble/'

# now we delete the output folder if they already exist
if (os.path.exists(ss_output_folder)): # ss output directory already exists
    os.removedirs(ss_output_folder) # remove the directory
if (os.path.exists(no_ss_output_folder)): # no ss output directory already exists
    os.removedirs(no_ss_output_folder) # remove the directory

# create the output folder
if not os.path.exists(ss_output_folder): # if the ss folder does not exist
    os.makedirs(ss_output_folder) # create the ss folder
os.chdir(ss_output_folder) # change current working directory
if not os.path.exists(no_ss_output_folder): # if the no ss folder does not exist
    os.makedirs(no_ss_output_folder) # create the no ss folder
os.chdir(no_ss_output_folder) # change current working directory 

# number of last characters (right to left) of the pdb file name to save
# as the name in the new ensemble, this will depend on your naming convention/preference
num_last_chars = 7

for i in range(len(ss_arr)): # create pdb ensemble for ss region conformers
    shutil.copy(str(ss_arr[i]), ss_output_folder + str(ss_arr[i])[len(ss_arr)-num_last_chars:]) # copy pdb file to new ensemble
for i in range(len(no_ss_arr)): # create pdb ensemble for no ss region conformers
    shutil.copy(str(no_ss_arr[i]), no_ss_output_folder + str(no_ss_arr[i][len(ss_arr)-num_last_chars:])) # copy pdb file to new ensemble

tf=time.perf_counter()
total_time=tf-ti
print("Program took " + str(total_time) + " seconds")