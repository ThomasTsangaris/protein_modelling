# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:48:20 2020

@author: Thomas Tsangaris
"""
import MDAnalysis as md
import numpy as np
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import math
import os.path
import shutil

from random import seed
from random import randint
from random import random
from random import uniform


# This script takes the compact population from 4E-BP2 ensemble from running 
# ENSEMBLE with only SAXS restraints and puts it into a folder called 
# "compact_conformers". Do the same thing for extended conformers in a folder
# called extended_conformers:


# From Greg's ENSEMBLE_polymer.py:
# 400 Conformer SAXS Only Relaxed:
# mypath = 'D:/0.0326_4_ensembles_Relaxed/'
# FFT Initial Pool of 8000 conf NP 4E-BP2 SAXS + PRE (N = 100):
# mypath = 'C:/Ensemble2/ENSEMBLE/NP_FFT_SAXS_PRE/100_conf_1/JOB01/Results/BEST.ENS_ANALYZE/Conformers/Weighted/'

# NP BP2 SAXS FFT 400 Conf:
mypath = '/mnt/c/NP_FFT_SAXS_400_conf/'

files=glob.glob(mypath+"*.pdb")
numfiles=len(files)
#Make sure to uncomment:
u_ens = md.Universe(files[0],files)
#print(u_ens.trajectory)
#print(files)
#print(u_ens)
#print(len(u_ens) == len(numfiles))
# print(numfiles)
# print(len(u_ens.trajectory)) #should be the number of frames (conformers) in
# the trajectory
# Print the topology:
#print(files[0])
# Print the first pdb analyzed that isn't the topology
#print(files[1]) 
bb = u_ens.select_atoms('protein and backbone')  # a selection (AtomGroup)
nterm = u_ens.select_atoms('protein and name CA')[1]
cterm = u_ens.select_atoms('protein and name CA')[-1]
Ree_ens=[]
Rg_ens=[]
A_ens=[]
Rhinv_ens=[]
ca=u_ens.select_atoms('protein and name CA')  # a selection (AtomGroup)
# for file in files:
#     print(file)
i = 0
compact_files = []
extended_files = []
for ts in u_ens.trajectory:     # iterate through all frames
    #print(ts)
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    # print(np.shape(bb[:].position))
    
    self_distances = distances.self_distance_array(ca.positions)
    Rhinv=(np.mean(self_distances**-1))
    # print(Rh)    
    
    if np.linalg.norm(r) < 100:
        Ree_ens.append(np.linalg.norm(r))  # end-to-end distance
    if bb.radius_of_gyration() < 40:
        Rg_ens.append( bb.radius_of_gyration())  # method of AtomGroup
    if bb.asphericity() < 0.5:
        A_ens.append(bb.asphericity())
    Rhinv_ens.append(Rhinv)
    
    # Append the pdbs that are in compact formation, assuming files and 
    # u_ens.select_atoms are parallel:
    if np.linalg.norm(r) < 100 and bb.radius_of_gyration() < 40 and \
        bb.asphericity() < 0.5:
            compact_files.append(files[i])
    elif np.linalg.norm(r) > 100 and bb.radius_of_gyration() > 40 and \
        bb.asphericity() > 0.5:
            extended_files.append(files[i])
    i += 1


# Check compact_files was created correctly:
numfiles=len(compact_files)
#Make sure to uncomment:
u_ens = md.Universe(compact_files[0], compact_files)
# print(numfiles)
# Print the topology:
print(files[0])
# Print the first pdb analyzed that isn't the topology
print(files[1]) 
bb = u_ens.select_atoms('protein and backbone')  # a selection (AtomGroup)
nterm = u_ens.select_atoms('protein and name CA')[1]
cterm = u_ens.select_atoms('protein and name CA')[-1]
Ree_ens=[]
Rg_ens=[]
A_ens=[]
Rhinv_ens=[]
ca=u_ens.select_atoms('protein and name CA')  # a selection (AtomGroup)
for ts in u_ens.trajectory:     # iterate through all frames
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    # print(np.shape(bb[:].position))
    
    self_distances = distances.self_distance_array(ca.positions)
    Rhinv=(np.mean(self_distances**-1))
    # print(Rh)    
    
    Ree_ens.append(np.linalg.norm(r))  # end-to-end distance
    Rg_ens.append( bb.radius_of_gyration())  # method of AtomGroup
    A_ens.append(bb.asphericity())
    Rhinv_ens.append(Rhinv)
    
# Now let's take away the "extended" region: 
compact_Rg_ens = []
compact_Ree_ens = []
compact_A_ens = []

for Rg in Rg_ens:
    if Rg < 40:
        compact_Rg_ens.append(Rg)
    # Should be no output:
    else:
        print(Rg)
for Ree in Ree_ens:
    if Ree < 100:
        compact_Ree_ens.append(Ree)
    # Should be no output:
    else:
        print(Ree)
for A in A_ens:
    if A < 0.5:
        compact_A_ens.append(A)
    # Should be no output:
    else:
        print(A)
        
# Probability Density Plot:
Ree_plot = sns.distplot(compact_Ree_ens, bins = 20, axlabel = 'Ree (A)', color = 'g', hist = True, kde=True)
print(compact_files)


# Important Note: You have to make a folder with the same name as the one 
# defined in save_path for this to work, otherwise you will make a file instead
# of a directory. 

# Folder to save pdb files:
save_path = '/mnt/e/NP_FFT_SAXS_compact_conformers'

# Put the compact conformers into the directory defined by save_path:
for file in compact_files:
    shutil.copy(file, save_path)
    
# Folder to save pdb files:
save_path = '/mnt/e/NP_FFT_SAXS_extended_conformers'

# Put the expanded conformers into the directory defined by save_path:
for file in extended_files:
    shutil.copy(file, save_path)