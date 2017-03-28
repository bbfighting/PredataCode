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
def replace_seq(str_seq):
    str_seq = str_seq.replace('T', 'U')
    str_seq = str_seq.replace('A', 'T')
    str_seq = str_seq.replace('U', 'A')
    str_seq = str_seq.replace('G', 'V')
    str_seq = str_seq.replace('C', 'G')
    str_seq = str_seq.replace('V', 'C')
    str_seq = str_seq.replace('t', 'u')
    str_seq = str_seq.replace('a', 't')
    str_seq = str_seq.replace('u', 'a')
    str_seq = str_seq.replace('g', 'v')
    str_seq = str_seq.replace('c', 'g')
    str_seq = str_seq.replace('v', 'c')
        
    return str_seq

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
def as_seq(species):
    as_coor_dic = {}
    as_fna_dic = {}

    if species == "hs":
        species = "human"
    f_coor = "original_data/{0:s}_nm_rna_region.coord".format(species)
    f_fna = "original_data/{0:s}_nm_rna.fna".format(species)
 
    # Open and read the file
    with open(f_coor, 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            data_list = re.split('\s+', line[0])
            cds_list = data_list[3].split('..')
            
            as_coor_dic[data_list[0]] = int(cds_list[0]) - 1
    fh.close()
    
    as_indic = ""

    # Open and read the file
    with open(f_fna, 'r') as fh:
        data = csv.reader(fh)
        
        for line in data:
            if line[0][0] == ">":
                temp = line[0][1:]
                as_fna_dic[line[0][1:]] = ""
                if as_indic != '':
                    as_fna_dic[as_indic] = as_fna_dic[as_indic][0:as_coor_dic[as_indic]]
            else:
                as_fna_dic[temp] += line[0]
                as_indic = temp
    as_fna_dic[as_indic] = as_fna_dic[as_indic][0:as_coor_dic[as_indic]]
    fh.close()

    return as_fna_dic 

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./main_uorf.py mouse

    as_fna_dic = as_seq('mouse')
    sys.exit(0)