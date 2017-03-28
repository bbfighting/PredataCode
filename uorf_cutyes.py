#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import sys, csv
import uorfcommon

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
def main(species):  #IREZone_size : at least IRES    
    as_fna_dic = uorfcommon.as_seq(species)
    
    if species == 'hs':
        file_species = "human"
    else:
        file_species = "mouse"

    fw = open("IRESdataset/{0:s}_five/csv/uorf_yes_{0:s}.csv".format(species), 'w')
    write = csv.writer(fw)

    with open("original_data/{0:s}_TISdb_data_1.0.csv".format(file_species), 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            if line[1] in as_fna_dic and line[-1] != "No":
                    write.writerows([line])
    fh.close()
    fw.close()

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./uorf_cutyes.py mouse
    species = sys.argv[1]

    main(species)
    sys.exit(0)