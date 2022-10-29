# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 00:25:41 2021

@author: Thomas Tsangaris
"""
import csv

# Purpose: To create the CS BME back-calculated and experimental files from 
# X-EISD's.

xeisd_cs_bc = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_CS_pavithra.txt' # X-EISD back-calculated file
xeisd_cs_exp = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_4EBP2_CS_XEISD.txt' # X-EISD experimental file
bme_cs_bc = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_CS_BME_calc_pavithra.dat' # file to create
bme_cs_exp = 'C:/Users/14169/Dropbox/PC/Documents/4E-BP2 Project/con_files/NP_CS_BME_exp_pavithra.dat'
with open(bme_cs_bc, 'a+') as bme_bc, open(xeisd_cs_bc, 'r') as xeisd_bc:
    bme_bc.write('#Label\n')
    reader = csv.reader(xeisd_bc, delimiter=',')
    j = 0
    for row in reader:
        bme_bc.write('frame_' + str(j) + '\t')
        for i in range(len(row)): # loop through values in a row by index:
            if row[i] == ' ': # skip spaces at the end of a row
                bme_bc.write('\n')
                continue
            bc_val = str(row[i])
            if i != (len(row) - 1): # not the last value in a row
                bme_bc.write(bc_val + '\t')
            else: # last value in a row
                 bme_bc.write(str(bc_val) + '\n')
        j += 1
    bme_bc.close()
    xeisd_bc.close()
    
with open(bme_cs_exp, 'a+') as bme_exp, open(xeisd_cs_exp, 'r') as xeisd_exp:
    bme_exp.write('# DATA=CS PRIOR=GAUSS\n')
    reader = csv.reader(xeisd_exp, delimiter=',')
    for row in reader:
        res_num = str(row[1]) # residue number 
        atom = str(row[2]) # atom letter(s)
        cs = str(row[3]) # chemical shift value
        unc = str(row[4]) # uncertainty in chemical shfit value
        line = atom + '_' + res_num + '\t' + cs + '\t' + unc + '\n' # line to write in BME experimental data format
        bme_exp.write(line)
    bme_exp.close()
    xeisd_exp.close()


