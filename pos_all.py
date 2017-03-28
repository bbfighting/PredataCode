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
    
    uorf_path = "IRESdataset/{0:s}_five/csv/pos_combine_sort_{0:s}.csv".format(species)
    xml_dir = "IRESdataset/{0:s}_five/xmlfile/".format(species)
    csv_path = "IRESdataset/{0:s}_five/csv/{0:s}_pos_all.csv".format(species)
    dirs = os.listdir(xml_dir)

    check_dic = {}
    uorf_dic = {}

    with open(uorf_path, 'r') as fh:
        data = csv.reader(fh)

        for line in data:
            if not line[0] in uorf_dic:
                uorf_dic[line[0]] = ["~".join(line[1:])]
                check_dic[line[0]] = [(int)(line[1])]
            else:
                uorf_dic[line[0]].append("~".join(line[1:]))
                
                if (int)(line[1]) < check_dic[line[0]][-1]:
                    check_dic[line[0]].append((int)(line[1]))
                else:
                    temp_uorf = uorf_dic[line[0]]
                    temp_check = check_dic[line[0]]
                    check_dic[line[0]] = []
                    uorf_dic[line[0]] = []
                    is_first = True
                    for i in range(len(temp_check)):
                        if (int)(line[1]) > temp_check[i]:
                            if is_first:
                                check_dic[line[0]].append((int)(line[1]))
                                uorf_dic[line[0]].append("~".join(line[1:]))
                                is_first = False
                        check_dic[line[0]].append(temp_check[i])
                        uorf_dic[line[0]].append(temp_uorf[i])
    fh.close()

    with open(csv_path, 'w') as fw:
        fwriter = csv.writer(fw)
        fwriter.writerows([["accession_number", "IRES", "IREZone", "uORF"]])

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
            NM_name = os.path.splitext(filename)[0]
                
            gene = root[0][0].text
            IRES_arr = []
            IREZone_arr = []
            check_IRES_list = []
            
            for i in range(IRES_count):
                if IRES[i].attrib['StartCoden'] == 1:
                    sc = "+1"
                else:
                    sc = IRES[i].attrib['StartCoden']
        
                temp_start = int(IRES[i].attrib['start'])
                temp_end = int(IRES[i].attrib['end'])
    
                write_start = temp_start + 1 + seq_len
                write_end = temp_end + 1 + seq_len
    
    
                if check_IRES_list != []:
                    if write_start > check_IRES_list[-1]:
                        temp_IRES_list = IRES_arr
                        temp_check_IRES_list = check_IRES_list
                        IRES_arr = []
                        check_IRES_list = []
                        is_first = True            
        
                        for j in range(len(temp_IRES_list)):
                            if write_start > temp_check_IRES_list[j]:
                                if is_first:
                                    check_IRES_list.append(write_start)
                                    IRES_arr.append(str(write_start) + "~" + str(write_end))
                                    is_first = False
                            check_IRES_list.append(temp_check_IRES_list[j])
                            IRES_arr.append(temp_IRES_list[j])
                    else:
                        IRES_arr.append(str(write_start) + "~" + str(write_end))
                        check_IRES_list.append(write_start)
                else:
                    IRES_arr.append(str(write_start) + "~" + str(write_end))
                    check_IRES_list.append(write_start)
    
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
    
                        write_IREZone_start = deviration_start + 1 + seq_len
                        write_IREZone_end = IREZone_end + 1 + seq_len
                        
    
                        IREZone_arr.append(str(write_IREZone_start) + "~" + str(write_IREZone_end))
                else:
                    if IREZone_count >= IREZone_size:
                        if IREZone_end == -1:
                            seq_end = None
                        else:
                            seq_end = IREZone_end + 1
                        IREZone_seq = sequence[deviration_start:seq_end]
                        IREZone_len = len(IREZone_seq)
    
    
                        write_IREZone_start = deviration_start + 1 + seq_len
                        write_IREZone_end = IREZone_end + 1 + seq_len
                        IREZone_arr.append(str(write_IREZone_start) + "~" + str(write_IREZone_end))
                    IREZone_count = 1
                    deviration_start = temp_start
                    deviration_end = temp_end
                    IREZone_end = temp_end

            uorf_key = uorf_dic.keys()

            if NM_name in uorf_key:         
                write_uorf = ";".join(uorf_dic[NM_name][::-1])
            else:
                write_uorf = ""
            fwriter.writerows([[NM_name, ";".join(IRES_arr[::-1]), ";".join(IREZone_arr[::-1]), write_uorf]])
    fw.close()

if __name__ == "__main__":
    #./pos_all.py mouse 25 3

    spe = sys.argv[1]
    diff = int(sys.argv[2])
    size = int(sys.argv[3])

    IREZoneCsv(spe, diff, size)
    sys.exit(0) 