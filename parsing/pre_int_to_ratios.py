# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 13:05:27 2022

@author: Thomas Tsangaris
"""

""" Purpose: Input files with oxidized and reduced PRE intensities in a specific 
 format, output a file in a specific format containing the intensity ratios

 Input: Oxidized and Reduced are in the format (inside the angle brackets):
     column 1 = <#> where # is the residue number of the label
     column 2 = <A#_...> where A and # are the residue letter and number respectively.
                The ellipsis indicates any characters. Assumes the only digits
                correspond to the residue number.
     column 3 = Fit Height
 Output: (format required by DEERPREdict scripts)
     column 1 = <#> where # is the residue of the label
     column 2 = <#> where # is the residue number
     column 3 = <#> where # is the PRE intensity peak ratio (oxidized/reduced)
     column 4 = <#> where # is the corresponding error to column 3 (does not work currently)
"""

import csv

ox_file = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5p_pre_oxidized_intensities.dat' # oxidized experimental file
red_file = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5p_pre_reduced_intensities.dat' # reduced experimental file

ratio_file = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/5P_PRE_Intensity_Height_Ratios_32_raw.dat' # output file with intensity ratios

# Oxidized
ox_int_lst = [] # list of experimental oxidized peak intensities
ox_res_lst = [] # parallel list containing residue numbers
ox_label_lst = [] # parallel list containing label residue number
hl=0 # how many header lines
with open(ox_file, newline = '') as f:                                                                       
    reader = csv.reader(f, delimiter='\t')
    for j in range(hl):
        #just skip headerlines, already stored
        next(reader)
    for row in reader:
        #Read remaining lines (row at a time)
        ox_label_lst.append(int(row[0]))
        res_num = int(''.join(filter(str.isdigit,str(row[1])))) # source: https://stackoverflow.com/questions/28526367/get-only-numbers-from-string-in-python
        ox_res_lst.append(res_num)
        ox_int_lst.append(float(row[2]))
    f.close()

# Reduced
red_int_lst = [] # list of experimental oxidized peak intensities
red_res_lst = [] # parallel list containing residue numbers
red_label_lst = [] # parallel list containing label residue number
hl = 0 # how many header lines
with open(red_file, newline = '') as f:                                                                              
    reader = csv.reader(f, delimiter='\t')
    for j in range(hl):
        #just skip headerlines, already stored
        next(reader)
    for row in reader:
        #Read remaining lines (row at a time)
        red_label_lst.append(int(row[0]))
        res_num = int(''.join(filter(str.isdigit,str(row[1])))) # source: https://stackoverflow.com/questions/28526367/get-only-numbers-from-string-in-python
        red_res_lst.append(res_num)
        red_int_lst.append(float(row[2]))
    f.close()

# Sanity Checks:
if ox_label_lst != red_label_lst:
    print("Oxidized and reduced files do not represent the same PRE experiment."\
          + " Label numbers do not match.")
if ox_res_lst != red_res_lst:
    print("Oxidized and reduced files do not represent the same PRE experiment."\
          + " Residue numbers do not match.")
if len(ox_int_lst) != len(red_int_lst):
    print("Oxidized and reduced files do not represent the same PRE experiment."\
          + " Length of files do not match.")

# Write to output file:
with open(ratio_file,'w',newline='\n') as f:
    for i in range(len(red_res_lst)):
        res_num = int(ox_res_lst[i])
        label_num = int(ox_label_lst[i])
        int_ratio = ox_int_lst[i] / red_int_lst[i] # compute ratio
        if int_ratio > 1:
            print("Intensity ratio at residue number", str(res_num), \
                  "from label residue number", str(label_num), \
                      "has an intensity ratio above one!")
        f.write('%d \t %d \t %f \t %f \n' % (label_num,res_num,int_ratio,0.05))
    f.close()