# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 21:54:52 2020

@author: Thomas Tsangaris
"""
import numpy as np


# In the folder where you put this also specify the chemical shifts file used


cs_data = open("cs_raw_data.txt", "r")


cs_values = []


for line in cs_data:
    cs_end = line.find(" ", 63)
    cs_start = line.rfind(" ", 0, cs_end) + 1
    cs = line[cs_start : cs_end]
    cs_values.append(float(cs))
  
    


cs_data = open("cs_raw_data.txt", "r")

atom_lst = []
with open('CS_4E-BP2.con', "w", newline='\n') as file:
    for line in cs_data:
        cs_end = line.find(" ", 63)
        cs_start = line.rfind(" ", 0, cs_end) + 1
        cs = float(line[cs_start : cs_end])
        atom_end = line.find(" ", 46)
        atom_start = 46
        atom = line[atom_start : atom_end]
        residue_end = line.find(" ", 35)
        residue_start = line.rfind(" ", 0, residue_end) + 1
        residue = line[residue_start : residue_end]
        if atom not in atom_lst:
            print(atom)
        atom_lst.append(atom)
        if atom != "CD" and atom != "CG" and atom != "CG2" and atom != "CD2" \
            and atom != "CD1" and atom != "CE" and atom != "HA2" and \
                atom != "HA3":
            # Note that I'm using the same uncertainties as the ones for Sic1
            # (0.4 and 0.03) because 4E-BP2 has no uncertainties with the 
            # chemical shifts data
            if atom[0] == "N" or atom[0] == "C":
                file.write("%s    %s    %f    %f"%(residue, atom, cs, 0.4)) 
            elif atom[0] == "H":
                file.write("%s    %s    %f    %f"%(residue, atom, cs, 0.03)) 
            file.write("\n")
        