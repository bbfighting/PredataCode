#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import re, sys, csv, sqlcommon


# ---------------- ---------------- ---------------- ---------------- #
# Global Variables
# ---------------- ---------------- ---------------- ---------------- #

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
def pos_add_gene(species, coor_dic):  #IREZone_size : at least IRES    
    filename = "IRESdataset/{0:s}_five/csv/pos_combine_sort_{0:s}.csv".format(species)
    fwname = "IRESdataset/{0:s}_five/csv/pos_combine_gene_{0:s}.csv".format(species)
    fw = open(fwname, 'w')
    write = csv.writer(fw)

    with open(filename, 'r') as fh:
        data = csv.reader(fh)

        for line in data:
            temp_list = [coor_dic[line[0]][0]]
            temp_list.extend(line)
            write.writerows([temp_list])
    fh.close()            
    fw.close()

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./sql_main.py mouse

    species = sys.argv[1]  
    coor_dic = sqlcommon.as_coor(species)

    pos_add_gene(species, coor_dic)
    sys.exit(0)