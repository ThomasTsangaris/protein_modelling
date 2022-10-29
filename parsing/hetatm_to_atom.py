# author: Thomas Tsangaris

# Description: This script creates a new folder filled with the same pdbs you supply
# but any atom registered as "HETATM" is renamed as "ATOM".

import glob
import os

# Path to your ensemble (folder of pdb files):
path = '/mnt/e/5P4E-BP2_FFT/FFT_1000/'

for i in range(1, 21):

    path = 'D:/5P4E-BP2_FFT/FFT_5P4E-BP2_SPO_Type' + str(i) + '/'   
    
    pdb_files = glob.glob(path + '*.pdb') # list of pdb files
    
    new_dir_name = 'Edited_PDBs/' # directory where the results will be placed
    
    if not os.path.exists(path + new_dir_name): # if the results folder does not exist
        os.makedirs(path + new_dir_name) # create the results folder
    os.chdir(path + new_dir_name) # change current working directory 
    
    j = 0
    
    for pdb in pdb_files:
    
        new_pdbs_lst = []
    
        lines = open(pdb, 'r').readlines()
        for i in range(len(lines)):
            if lines[i][:6] == 'HETATM':
                 lines[i] = 'ATOM    ' + lines[i][8:]
        new_pdbs_lst.append(lines)
        
        for new_pdb_lines in new_pdbs_lst:
            with open(path + new_dir_name + str(j) + '_edited.pdb', 'w') as new_pdb:
                new_pdb.writelines(new_pdb_lines)
            j += 1
        
print('DONE')