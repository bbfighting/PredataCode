#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import sys, csv
import uorfcommon, uorf_parse_mouse, subseq

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
def main(method, species):  #IREZone_size : at least IRES    
    as_fna_dic = uorfcommon.as_seq(species)

    fw = open("IRESdataset/{1:s}_five/csv/pos_{0:s}_{1:s}.csv".format(method, species), 'w')
    write = csv.writer(fw)
    nomatch_up1000_dic = {}

    with open("IRESdataset/{1:s}_five/csv/nomatch_{0:s}_rel_{1:s}.csv".format(method, species), 'r') as fh:
        data = csv.reader(fh)

        for line in data:
            temp_seq = subseq.sub_seq(as_fna_dic, line[1], line[3], line[5], species).lower()
            temp_start_c = temp_seq[0:3]
            temp_stop_c = temp_seq[-3::]

            if temp_start_c != line[4].lower() or temp_stop_c != line[6].lower():
                print line
                is_match = uorf_parse_mouse.compare_sub_rel_seq(line[1], line[2], line[3], line[5], line[4].lower(), line[6].lower(), "{0:s}{1:s}".format(as_fna_dic[line[1]].lower(), "atg"), write)

                if not is_match:
                    if not line[0] in nomatch_up1000_dic:
                        nomatch_up1000_dic[line[0]] = [line]
                    else:
                        nomatch_up1000_dic[line[0]].extend([line])
            else:
                write.writerows([[line[1], line[3], line[5]]])
    fh.close()
    fw.close()

    fw = open("IRESdataset/{1:s}_five/csv/nomatch_{0:s}_{1:s}.csv".format(method, species), 'w')
    write = csv.writer(fw)

    key = nomatch_up1000_dic.keys()
    key.sort()

    for i in key:
        for j in nomatch_up1000_dic[i]:
            temp_list = []
            temp_list.extend(j)
            write.writerows([temp_list])
    fw.close()

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./uorf_rel_main.py parse mouse
#./uorf_rel_main.py up1000 mouse
    method = sys.argv[1]
    species = sys.argv[2]

    main(method, species)
    sys.exit(0)
