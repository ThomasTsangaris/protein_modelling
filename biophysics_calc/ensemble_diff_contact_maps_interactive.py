# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 15:57:54 2020

@author: Thomas Tsangaris
"""
from typing import List, Tuple, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.dihedrals import Ramachandran
from contact_map import ContactMap, ContactFrequency, ContactDifference, AtomMismatchedContactDifference, OverrideTopologyContactDifference
import time

start_time = time.time() # start time

# NOTE: THIS PROGRAM IS INTERACTIVE

# Important Note: This script only works if you have a folder with pdb files 
# that define your ensemble.




def ensemble_contacts(path: str, n: int) -> Optional[list]:
    """Return a list of tuples containing a list of the residue pair that 
    forms a contact as the first element in the tuple and the fraction of 
    frames where that contact occurs in the trajectory created from the 
    ensemble defined by <path> for the second element. Residues numbered i and 
    j will have their contact considered if their contact satisifes the 
    inequality: |i-j| > <n>.
    
    Precondition: <path> contains pdb files which make up an ensemble of a 
    protein.
    """
    # Get the pdb files that define the ensemble:
    try:
        files = glob.glob(path+"*.pdb")
        # Load the trajectory:
        traj = md.load(files)
    except:
        return None
        print("\nNot a valid path.\n")
    # Finds the fraction of frames where each contact exists:
    trajectory_contacts = ContactFrequency(traj, n_neighbors_ignored = n)
    # A list for the higher frequency contacts: 
    frequent_contacts = []
    # Fill up frequent_contacts:
    for contact in trajectory_contacts.residue_contacts.most_common():
        if contact[1] > 0:
            frequent_contacts.append(contact)
    return frequent_contacts




#################################
# Guide user through the process:
#################################

print("Would you like to analyze one ensemble (input 1) " + 
      "or compare two ensembles (input 2)?")
setting = input()

if setting == "1":
    print("This program will give you the most common contacts in your" +
          " ensemble. Please input the path to the folder containing the pdb" + 
              " files which define the ensemble you want to analyze " + 
              "(ex. D:/ENSEMBLE/ensemble_to_analyze/):")
    input_path = input()
    
    print("Input a number n such that residues with residue numbers i and j " + 
          "will have their contact considered if |i-j| > n is true " + 
          "(default is 10):")
    try: 
        n = int(input())
    except:
        n = 10
    print("Loading...")
    frequent_contacts = ensemble_contacts(input_path, n)
    # Check if the path is invalid:
    while frequent_contacts == None:
        print("The path to your folder is invalid. Try again:")
        input_path = input()
        print("Loading...")
        frequent_contacts = ensemble_contacts(input_path, n)
    
    
    
    # print("Input the name of the protein:")
    # name = input()
    
    if frequent_contacts == []:
        print("No contacts with a frequency greater than 0.1 found.")
        print("Loading...")
    else:
        print("Here are all contacts with a frequency greater than 0.1 " + 
          "(ie. the fraction that the contact occurs in the trajectory is " +
          "greater than 0.1):")
        print('')
        print(frequent_contacts)
        print("Loading...")
        
    
    # Make a "frequency plot":
    fraction_lst = []
    # Full trajectory:
    files = glob.glob(input_path + "*.pdb")
    traj = md.load(files, stride = 1)
    topology = traj.topology
    trajectory_contacts = ContactFrequency(traj, n_neighbors_ignored = n)
    
    # Get residues which make contact:
    residues_in_contact = []
    for contact in frequent_contacts:
        # Append both residues involved in the contact:
        residues_in_contact.append(int(str(contact[0][0])[3:]))
        residues_in_contact.append(int(str(contact[0][1])[3:]))
    
        # Append two frequencies for both residues:
        fraction_lst.append(float(contact[1]))
        fraction_lst.append(float(contact[1]))
    
    # Get most common contacts at a chosen residue:
    print('Would you like to analyze a residue for its most common contacts? ' + 
          '(Answer "yes" or "no", default is no)')
    analyze_contacts = input()
    if analyze_contacts.lower() == 'yes' or analyze_contacts.lower() == 'y':
        print("Input the residue number of the amino acid you would like to " +
              "analyze:")
        residue_num = input()
        while not residue_num.isnumeric():
            print("That is not a valid residue number. Try again:")
            residue_num = input()
        # Print all of the important contacts:
        residue = topology.residue(int(residue_num) - 1) # This is how the 
                                                         # numbering works for pdbs
        aa_contact_lst = [] 
        for contact in \
            trajectory_contacts.residue_contacts.most_common(residue):
            if contact[1] > 0.1:
                aa_contact_lst.append(contact)
        if aa_contact_lst == []:
            print("No contacts containing " + str(residue) + 
                  " with a frequency greater than 0.1 found.")
        else:
            print("Here are all contacts with a frequency greater than 0.1 that "+ 
              "involve " + str(residue) + " (ie. the fraction that the contact " + 
              "occurs in the trajectory is greater than 0.1):")
            print('')
            for contact in aa_contact_lst:
                print(contact)
                
    
    
    # sns.jointplot cannot plot parallel lists of 2 or less elements with kde, 
    # switch to scatter plot if this is the case:
    
    if len(fraction_lst) <= 2:
        print("First plot is to give you an idea which residues make contacts:")
        plot = sns.jointplot(residues_in_contact, fraction_lst, kind="scatter", \
                             color="g", xlim=[0, 120])
        plot.ax_joint.set_xlabel(r'Residue', fontweight='bold')
        plot.ax_joint.set_ylabel(r'Frequency', fontweight='bold')
        
    # Otherwise, make the plot with kde:
    else:
        print("First plot is to give you an idea which residues make contacts:")
        plot = sns.jointplot(residues_in_contact, fraction_lst, kind="kde", \
                             color="g", xlim=[0, 120])
        plot.ax_joint.set_xlabel(r'Residue', fontweight='bold')
        plot.ax_joint.set_ylabel(r'Frequency', fontweight='bold')
    
    print("Second plot is a contact map for your entire ensemble:")
    # Make a contact map of the full trajectory:
    fig, ax = trajectory_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, 
                                                      vmax=1)
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    
    print("Third plot is a contact map for the zeroth frame (first conformer):")
    #frame_contacts = ContactMap(traj[0])
    #frame_contacts = ContactMap(traj[1])
    #frame_contacts = ContactMap(traj[2])
    #frame_contacts = ContactMap(traj[3])
    frame_contacts = ContactMap(traj[0])
    fig, ax = frame_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, 
                                                      vmax=1)



else:
    print("This program will give you a comparison of the contacts found in" +
          " two ensembles:")
    print("Input the path to the first ensemble " + 
          "(ex. D:/ENSEMBLE/ensemble_to_analyze/):")
    path_1 = input()
    print("Input the path to the second ensemble:")
    path_2 = input()
    print("Input a number n such that residues with residue numbers i and j " + 
          "will have their contact considered if |i-j| > n is true " + 
          " for the first ensemble you entered (default is 10):")
    try: 
        n_1 = int(input())
    except:
        n_1 = 10
    print("Input a number n for the second ensemble: (default is 10)")
    try: 
        n_2 = int(input())
    except:
        n_2 = 10
    print("Loading...")
    co = 0.45 #cut off (default 0.45 nm)
    # co = 0.8 # for Calpha, "usually used": doi - 10.1007/978-1-4939-3572-7_24
    # Get contacts for the first ensemble:
    files_1 = glob.glob(path_1 + "*.pdb")
    traj_1 = md.load(files_1, top = files_1[0], stride = 1)
    topology_1 = traj_1.topology
    trajectory_contacts_1 = ContactFrequency(traj_1, n_neighbors_ignored = n_1, cutoff = co)
    print("trajectory_1 is:", trajectory_contacts_1)
    #Get contacts for the second ensemble:
    files_2 = glob.glob(path_2 + "*.pdb")
    traj_2 = md.load(files_2, top = files_1[0], stride = 1)
    topology_2 = traj_2.topology
    trajectory_contacts_2 = ContactFrequency(traj_2, n_neighbors_ignored = n_2, cutoff = co)
    print("trajectory_2 is:", trajectory_contacts_2)
    
    # Now get the difference between the contacts:
    # diff = trajectory_contacts_1 - trajectory_contacts_2
    # ignore residues that do not match exactly
    diff = AtomMismatchedContactDifference(trajectory_contacts_1, trajectory_contacts_2)
    # diff = OverrideTopologyContactDifference(trajectory_contacts_1, trajectory_contacts_2, files_1[0])
    
    print("Here is a list of the most common fractions found in the " + 
          "difference contact map (only positive differences greater than" + 
          " 0.05 are shown):")
    
    print('')
    
    diff_most_common = diff.residue_contacts.most_common()
    frequent_contacts = []
    # Fill up frequent_contacts:
    for contact in diff_most_common:
        if abs(contact[1]) > 0.05:
            frequent_contacts.append(contact)

    print(frequent_contacts)
    
    print("The first plot is a contact map of the difference of " + 
          "frequencies for the two ensembles. Red means that the " + 
          "contact occurs in the first ensemble you inputted but not the " +
          "second, blue is the opposite:")
    
    # Plot the contact map:
    # Change vmin and vmax to actually see what's going on
    # (fig, ax) = diff.residue_contacts.plot(cmap='seismic', vmin=-0.1, vmax=0.1)
    # plt.xlabel("Residue")
    # plt.ylabel("Residue")

    # (fig, ax) = diff.residue_contacts.plot(cmap='seismic', vmin=-0.05, vmax=0.05)
    # plt.xlabel("Residue")
    # plt.ylabel("Residue")

    
    S=diff.residue_contacts.sparse_matrix
    S=S.toarray()
    plt.figure()
    #plt.imshow(S, interpolation=None, origin='low',cmap='seismic')
    # plt.xlim([18, 62])
    # plt.ylim([63, 121])
    plt.imshow(S, interpolation=None, origin='lower',cmap='seismic',vmin=-0.05,vmax=0.05)
    plt.colorbar()
    plt.xlabel("Residue")
    plt.ylabel("Residue")


print("Program took %s seconds to finish" % (time.time() - start_time))
