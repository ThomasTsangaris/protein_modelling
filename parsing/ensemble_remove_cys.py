# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:51:31 2020

@author: Thomas Tsangaris
"""
from typing import List

import MDAnalysis as md
import numpy as np
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import math

from random import seed
from random import randint
from random import random
from random import uniform



# Program was made to remove the mutant cysteine attached at the C-terminus 
# required for PRE measurements from the conformers generated from a PRE 
# restrained ENSEMBLE calculation



#mypath = 'D:/ENSEMBLE/p4E-BP2_PRE_and_SAXS_400_Conf/'
#mypath = 'D:/ENSEMBLE/p4E-BP2_PRE_400_Conf/'
mypath = 'D:/FloppyTail_Structures_500_3mer_list/'

folder_lst = []

for i in range(0, 499):
    folder_lst.append('D:/ENSEMBLE/FloppyTail_Structures_500_3mer_list/FloppyTail_out_' + str(i) + '/')


files = glob.glob(mypath + "*.pdb")

i = 0
print(files)
for file in files:
    pdb = open(file, 'r')
    with open('D:/FloppyTail_Structures_500_3mer_list/new_FloppyTail_out_' + str(i) + '.pdb', "w") as new_pdb:
        lines_lst = pdb.readlines()
        #print(len(lines_lst))
        new_pdb.writelines(lines_lst[:1796])
    i += 1

