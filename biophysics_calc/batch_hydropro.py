# author: Thomas Tsangaris

# Description: This script runs HYDROPRO on a folder of pdb files 
# supplied by the user as the first (and only) argument

# Sample use in linux terminal: python3 batch_hydropro.py /full_path_to_ensemble_folder/


import glob
import sys
import subprocess
import os
import shutil
import math
import statistics
import time

# Helper Function:
def write_pdb_to_hydropro_dat(pdb_name: str, dat_file: str) -> None:
    """Append the seventeen required lines that denote the pdb file 
    representing <pdb_name> to the file <dat_file>.

    Returns: None
    """
    with open(dat_file, mode = 'a+') as dat_file:
        lines = []

        # Change these to your liking: 

        # From Greg's Paper:
        lines.append(pdb_name + ' !Name of molecule\n')
        lines.append(pdb_name + ' !Name for output file\n')
        lines.append(pdb_name + '.pdb' + ' !Structural (PDB) file\n')
        lines.append('1' + ' !Type of calculation\n')
        lines.append('2.9,' + ' !AER, radius of primary elements\n')
        lines.append('6,' + ' !NSIG\n')
        lines.append('1.0,' + ' !Minimum radius of beads in the shell (SIGMIN)\n')
        lines.append('2.0,' + ' !Maximum radius of beads in the shell (SIGMAX)\n')
        lines.append('20.,' + ' !T (temperature, centigrade)\n')
        lines.append('0.01,' + ' !ETA (Viscosity of the solvent in poises)\n')
        lines.append('10000.,' + ' !RM (Molecular weight)\n')
        lines.append('0.702,' + ' !Partial specific volume, cm3/g\n')
        lines.append('1.0,' + ' !Solvent density, g/cm3\n')
        lines.append('-1' + ' !Number of values of Q\n')
        lines.append('-1' + ' !Number of intervals\n')
        lines.append('0,' + ' !Number of trials for MC calculation of covolume\n')
        lines.append('1' + ' !IDIF=1 (yes) for full diffusion tensors\n')
        

        # ENSEMBLE: 
        # lines.append(pdb_name + ' !Name of molecule\n')
        # lines.append(pdb_name + ' !Name for output file\n')
        # lines.append(pdb_name + '.pdb' + ' !Structural (PDB) file\n')
        # lines.append('1' + ' !Type of calculation\n')
        # lines.append('3.1,' + ' !AER, radius of primary elements\n')
        # lines.append('2,' + ' !NSIG\n')
        # lines.append('1.7,' + ' !Minimum radius of beads in the shell (SIGMIN)\n')
        # lines.append('2.0,' + ' !Maximum radius of beads in the shell (SIGMAX)\n')
        # lines.append('293.,' + ' !T (temperature, centigrade)\n')
        # lines.append('0.01,' + ' !ETA (Viscosity of the solvent in poises)\n')
        # lines.append('14320.,' + ' !RM (Molecular weight)\n')
        # lines.append('0.702,' + ' !Partial specific volume, cm3/g\n')
        # lines.append('1.0,' + ' !Solvent density, g/cm3\n')
        # lines.append('-1' + ' !Number of values of Q\n')
        # lines.append('-1' + ' !Number of intervals\n')
        # lines.append('10,' + ' !Number of trials for MC calculation of covolume\n')
        # lines.append('1' + ' !IDIF=1 (yes) for full diffusion tensors\n')
        
        # Testing Extreme Effects:
        # lines.append(pdb_name + ' !Name of molecule\n')
        # lines.append(pdb_name + ' !Name for output file\n')
        # lines.append(pdb_name + '.pdb' + ' !Structural (PDB) file\n')
        # lines.append('1' + ' !Type of calculation\n')
        # lines.append('3.1,' + ' !AER, radius of primary elements\n')
        # lines.append('6,' + ' !NSIG\n')
        # lines.append('1.0,' + ' !Minimum radius of beads in the shell (SIGMIN)\n')
        # lines.append('2.0,' + ' !Maximum radius of beads in the shell (SIGMAX)\n')
        # lines.append('20.,' + ' !T (temperature, centigrade)\n')
        # lines.append('0.1,' + ' !ETA (Viscosity of the solvent in poises)\n')
        # lines.append('15000.,' + ' !RM (Molecular weight)\n')
        # lines.append('1.2,' + ' !Partial specific volume, cm3/g\n')
        # lines.append('2.0,' + ' !Solvent density, g/cm3\n')
        # lines.append('-1' + ' !Number of values of Q\n')
        # lines.append('-1' + ' !Number of intervals\n')
        # lines.append('0,' + ' !Number of trials for MC calculation of covolume\n')
        # lines.append('1' + ' !IDIF=1 (yes) for full diffusion tensors\n')

        dat_file.writelines(lines) # write the seventeen prequiste lines
        dat_file.close() # close the dat file


# Main:
if __name__ == "__main__":
    start_time = time.time() # start time

    dir_path = sys.argv[1] # path to folder of pdb files

    # First, we must create the hydropro.dat file, iterate through 
    # each pdb file and add it to the file
    pdb_files = glob.glob(dir_path + '*.pdb') # list of pdb files
    results_path = dir_path + 'HYDROPRO_Results/' # directory where the results will be placed

    if not os.path.exists(results_path): # if the results folder does not exist
        os.makedirs(results_path) # create the results folder
    os.chdir(results_path) # change current working directory 

    hydropro_dat = 'hydropro.dat' # name of the dat file to create 
    line_num = 0

    if (os.path.exists(results_path + hydropro_dat)): # if the hydropro.dat file exists
        os.remove(results_path + hydropro_dat) # remove the hydropro.dat file

    for pdb in pdb_files: # iterate through all pdb files
        shutil.copy(pdb, results_path) # copy the pdb files into the results directory
        write_pdb_to_hydropro_dat(pdb[len(dir_path):-4], hydropro_dat) # helper function

    bash_cmd = ['hydropro10-lnx.exe'] # bash command, change according to the name of the executable you downloaded
    print('Running HYDROPRO on ' + str(len(pdb_files)) + ' pdb files')
    process = subprocess.run(bash_cmd, stdout = subprocess.PIPE) # run HYDROPRO
    print('HYDROPRO calculations completed')

    # Now, we must loop through every res file, parse it for the required values, calculate the hydrodynamic 
    # radius, then find the average and standard deviation of all those values 
    R_h_lst = []

    txt_files = glob.glob(results_path + '*.txt') # list of all text files
    res_files = [] # list of all res files
    
    for txt in txt_files: # go through all text files
        if txt[-8:-4] == '-res': # is a res file
            res_files.append(txt) # add to list of res files

    all_Rh_file = results_path + 'Rh_Values_for_' + str(len(pdb_files)) + '_conformers.txt'
    if (os.path.exists(all_Rh_file)): # if the file exists
        os.remove(all_Rh_file) # remove the file
    
    with open(all_Rh_file, mode = 'a+') as all_Rh_file:
        for res in res_files: # loop through each res file
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
    
            R_h_lst.append(R_h) # add pdb's R_h value to the list
            
            all_Rh_file.write(str(R_h) + '\n') # write Rh to the text file
            
        all_Rh_file.close() # close the text file

    Rh_ensemble_stats = results_path + 'ensemble_statistics.txt'
    
    if (os.path.exists(Rh_ensemble_stats)): # if the file exists
        os.remove(Rh_ensemble_stats) # remove the file
    
    with open(Rh_ensemble_stats, mode = 'a+') as Rh_ensemble_stats:
        print(R_h_lst) #debug
        print("Before Reweighting:")
        Rh_ensemble_stats.write('Before Reweighting:' + '\n')
        if (len(R_h_lst) >= 1):
            mean = statistics.mean(R_h_lst) # Mean R_h in Angstroms 
            Rh_ensemble_stats.write('Mean: ' + str(mean) + '\n')
        if len(R_h_lst) > 1:
            std = statistics.stdev(R_h_lst) # Standard Deviation for R_h in Angstroms
            Rh_ensemble_stats.write("Standard Deviation: " + str(std) + '\n')
            std_of_mean = std / math.sqrt(len(R_h_lst)) # Standard Deviation of the Mean for R_h in Angstroms
            Rh_ensemble_stats.write("Standard Deviation of the Mean: " + str(std_of_mean) + '\n\n')
        if (len(R_h_lst) >= 1):
            print('Mean Rh: ' + str(mean) + '\n')
        if len(R_h_lst) > 1:
            print('STD Rh: ' + str(std) + '\n')
            print('SDOM Rh: ' + str(std_of_mean) + '\n')
        Rh_ensemble_stats.close()

print("Program took %s seconds" % (time.time() - start_time))
