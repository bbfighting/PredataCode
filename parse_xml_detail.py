#!/usr/bin/env python

"""
This script is used to parse the stock information
difference from parse_xmlis is csv format 
"""

import os, sys, csv
import xml.etree.ElementTree as ET
from time import strftime

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
def IREZoneCsv(species, diff_const, IREZone_size):  #IREZone_size : at least IRES
    
    xml_dir = "IRESdataset/{0:s}_five/xmlfile/".format(species)
    csv_path = "IRESdataset/{0:s}_five/csv/{0:s}_five_IREZone.csv".format(species)
    dirs = os.listdir(xml_dir)

    with open(csv_path, 'w') as fh:
        fwriter = csv.writer(fh)
        fwriter.writerows([["accession_number", "seq_len", "IRES", "IREZone"]])

        for filename in dirs:
            xml_path = xml_dir + filename

            tree = ET.parse(xml_path)
            root = tree.getroot()
            deviration_start = 0
            IREZone_count = 0
            IRES_count = len(root[1])
            IRES = root[1]
            sequence = root[0][1].text
            seq_len = len(sequence)
                
            gene = root[0][0].text
            IRES_arr = []
            IREZone_arr = []

            for i in range(IRES_count):
                if IRES[i].attrib['StartCoden'] == 1:
                    sc = "+1"
                else:
                    sc = IRES[i].attrib['StartCoden']
        
                temp_start = int(IRES[i].attrib['start'])
                temp_end = int(IRES[i].attrib['end'])

                IRES_arr.append(str(temp_start) + ":" + str(temp_end))

                if (deviration_start <= temp_end):
                    deviration = (deviration_start - temp_start) + (deviration_end - temp_end)

                    deviration_start = min(temp_start, deviration_start)
                
                    if deviration_start >= temp_start and deviration > diff_const:
                        deviration_end = temp_end    
                        IREZone_count += 1

                    if i == IRES_count - 1 and IREZone_count >= IREZone_size:
                        if IREZone_end == -1:
                            seq_end = None
                        else:
                            seq_end = IREZone_end + 1
                        IREZone_seq = sequence[deviration_start:seq_end]
                        IREZone_len = len(IREZone_seq)
                        IREZone_arr.append(str(deviration_start) + ":" + str(IREZone_end))
                else:
                    if IREZone_count >= IREZone_size:
                        if IREZone_end == -1:
                            seq_end = None
                        else:
                            seq_end = IREZone_end + 1
                        IREZone_seq = sequence[deviration_start:seq_end]
                        IREZone_len = len(IREZone_seq)
                        IREZone_arr.append(str(deviration_start) + ":" + str(IREZone_end))

                    IREZone_count = 1
                    deviration_start = temp_start
                    deviration_end = temp_end
                    IREZone_end = temp_end
            fwriter.writerows([[os.path.splitext(filename)[0], seq_len, ";".join(IRES_arr), ";".join(IREZone_arr)]])
    fh.close()

if __name__ == "__main__":
    #./parse_xml_detail.py mouse 25 3

    spe = sys.argv[1]
    diff = int(sys.argv[2])
    size = int(sys.argv[3])

    IREZoneCsv(spe, diff, size)
    sys.exit(0) 