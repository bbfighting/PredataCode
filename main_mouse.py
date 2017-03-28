#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import sys
import exec_pat_mouse, parse_xml

# ---------------- ---------------- ---------------- ---------------- #
# Global Variables
# ---------------- ---------------- ---------------- ---------------- #
info_file = "original_data/" + sys.argv[1]
seq_file = "original_data/" + sys.argv[2]
species = sys.argv[3]
cds35_method = sys.argv[4]

diff = int(sys.argv[5])
IRES_up = int(sys.argv[6])

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
def main():
    gene_region = exec_pat_mouse.cds35(info_file)
    exec_pat_mouse.pre_content(seq_file, gene_region, species, cds35_method)

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./main_150803.py mouse_nm_rna_region.coord sub.fna predict five 25 3
#./main_150803.py mouse_nm_rna_region.coord mouse_nm_rna.fna mouse five 25 3

    main()
    sys.exit(0)