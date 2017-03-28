#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import sys, uorfcommon

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
def sub_seq(as_fna_dic, acc, pos1, pos2, species):
    temp_seq = "{0:s}{1:s}".format(as_fna_dic[acc], "ATG")
    return temp_seq[int(pos1)-1:int(pos2)]

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
    acc = sys.argv[1]
    pos1 = sys.argv[2]
    pos2 = sys.argv[3]
    species = sys.argv[4]
    as_fna_dic = uorfcommon.as_seq(species)
    print sub_seq(as_fna_dic, acc, pos1, pos2, species)

    sys.exit(0)