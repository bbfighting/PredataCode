#!/usr/bin/env python

"""
This script is used to parse the stock information
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
def IREZoneCsv(species, diff_const, IREZone_size, method):  #IREZone_size : at least IRES
    
    xml_dir = "IRESdataset/{0:s}/xmlfile/".format(species + "_" + method)
    csv_path = "IRESdataset/{0:s}_{1:s}/csv/{0:s}_{1:s}_IREZone.csv".format(species, method)
    dirs = os.listdir(xml_dir)

    with open(csv_path, 'w') as fh:
        fwriter = csv.writer(fh)
        fwriter.writerows([["accession_number", "gene", "seq_len", "IREZone_seq", "IREZone_len", "start", "end", "Ires_count"]])

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

            for i in range(IRES_count):
                if IRES[i].attrib['StartCoden'] == 1:
                    sc = "+1"
                else:
                    sc = IRES[i].attrib['StartCoden']
        
                temp_start = int(IRES[i].attrib['start'])
                temp_end = int(IRES[i].attrib['end'])

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
                        fwriter.writerows([[os.path.splitext(filename)[0], gene, seq_len, IREZone_seq, IREZone_len, deviration_start, IREZone_end, IREZone_count]])
                else:
                    if IREZone_count >= IREZone_size:
                        if IREZone_end == -1:
                            seq_end = None
                        else:
                            seq_end = IREZone_end + 1
                        IREZone_seq = sequence[deviration_start:seq_end]
                        IREZone_len = len(IREZone_seq)
                        fwriter.writerows([[os.path.splitext(filename)[0], gene, seq_len, IREZone_seq, IREZone_len, deviration_start, IREZone_end, IREZone_count]])


                    IREZone_count = 1
                    deviration_start = temp_start
                    deviration_end = temp_end
                    IREZone_end = temp_end
    fh.close()

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
def output_IREZone(start, end, size):
    print "\nIREZone_start : " + str(start) #IREZone_start = deviration_start
    print "IREZone_end : " + str(end) 
    print "IREZone_size : " + str(size) + "\n"

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
def output_IRES(start, end, deviration_start, deviration_end, count, deviration = ""):
    print "\ndeviration : " + str(deviration)
    print "new_start : " + str(start)
    print "new_end : " + str(end)
    print "deviration_start : " + str(deviration_start)
    print "deviration_end : " + str(deviration_end)
    print "IREScount : " + str(count)

if __name__ == "__main__":
    #./parse_xml.py predict 25 3 web
    #./parse_xml.py hs 25 3 five

    spe = sys.argv[1]
    diff = int(sys.argv[2])
    size = int(sys.argv[3])
    method = sys.argv[4]

    IREZoneCsv(spe, diff, size, method)
    sys.exit(0) 