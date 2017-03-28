#!/usr/bin/env python

import sys, csv, re

# ---------------- ---------------- ---------------- ---------------- #
# Function 
#
# Parameters:
#   
#
# Return:
#   
#
# The description.
# ---------------- ---------------- ---------------- ---------------- #
def as_coor(species):
    as_coor_dic = {}
    if species == "hs":
        species = "human"
    f_coor = "original_data/{0:s}_nm_rna_region.coord".format(species)

 
    # Open and read the file
    with open(f_coor, 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            data_list = re.split('\s+', line[0])
            cds_list = data_list[3].split('..')
            
            as_coor_dic[data_list[0]] = [data_list[1], int(cds_list[0]) - 1]
    fh.close()

    return as_coor_dic

# ---------------- ---------------- ---------------- ---------------- #
# Function 
#
# Parameters:
#   
#
# Return:
#   
#
# The description.
# ---------------- ---------------- ---------------- ---------------- #
def as_fna(as_coor_dic, species):
    as_fna_dic = {}
    if species == "hs":
        species = "human"
    f_fna = "original_data/{0:s}_nm_rna.fna".format(species)
    
    as_indic = ""

    # Open and read the file
    with open(f_fna, 'r') as fh:
        data = csv.reader(fh)
        
        for line in data:
            if line[0][0] == ">":
                temp = line[0][1:]
                as_fna_dic[line[0][1:]] = ""
                if as_indic != '':
                    as_fna_dic[as_indic] = as_fna_dic[as_indic][0:as_coor_dic[as_indic][1]]
            else:
                as_fna_dic[temp] += line[0]
                as_indic = temp
    as_fna_dic[as_indic] = as_fna_dic[as_indic][0:as_coor_dic[as_indic][1]]
    fh.close()

    return as_fna_dic 

# ---------------- ---------------- ---------------- ---------------- #
# Function 
#
# Parameters:
#   
#
# Return:
#   
#
# The description.
# ---------------- ---------------- ---------------- ---------------- #
def gene_list(coor_dic):
    temp_list = []
    for acc in coor_dic.keys():
        gene = coor_dic[acc][0]
        if not gene in temp_list:
            temp_list.append(gene)
    
    temp_list.sort()
    
    return temp_list