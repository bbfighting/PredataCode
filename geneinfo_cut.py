#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import re, sys, csv

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
def cut(species):  #IREZone_size : at least IRES    
    if species == 'hs':
        tax_id = '9606'
    elif species == 'mouse':
        tax_id = '10090'

    filename = "original_data/All_Data.gene_info"

    fw = open("IRESdataset/{0:s}_five/csv/geneinfo_{0:s}.csv".format(species), 'w')
    write = csv.writer(fw)

    # Open and read the file
    with open(filename, 'r') as fh:
        fh.next()

        for line in fh:
            detail = re.split('\t', line)
            if detail[0] == tax_id:
                write.writerows([detail])
    fh.close()
    fw.close()

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./geneinfo_cut.py mouse

    species = sys.argv[1]
    cut(species)
    
    sys.exit(0)
