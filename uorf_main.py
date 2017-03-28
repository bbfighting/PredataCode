#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import sys
import uorf_parse_mouse

# ---------------- ---------------- ---------------- ---------------- #
# Global Variables
# ---------------- ---------------- ---------------- ---------------- #
species = sys.argv[1]


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
    uorf_parse_mouse.parse_pos(species)

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./main_uorf.py mouse

    main()
    sys.exit(0)