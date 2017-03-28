#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import os, sys, csv, re
import xml.etree.ElementTree as ET
import urllib2
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
def parse_pos(species):  #IREZone_size : at least IRES    
    as_fna_dic = uorfcommon.as_seq(species)

    fw = open("IRESdataset/{0:s}_five/csv/pos_uorf_{0:s}.csv".format(species), 'w')
    write = csv.writer(fw)
    up1000_dic = {}
    nomatch_dic = {}

    if species == "hs":
        file_species = "human"
        dataset = "hg19"
    else:
        file_species = species
        dataset = "mm10"

    with open("original_data/{0:s}_TISdb_data_1.0.csv".format(file_species), 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            if line[1] in as_fna_dic and line[-1] != "No":
                is_match = False
                if abs(int(line[6]) - int(line[4])) < 1000:
                    if line[3] == "+":
                        url = "http://genome.ucsc.edu/cgi-bin/das/{0:s}/dna?segment={1:s}:{2:s},{3:s}".format(dataset, line[2], line[4], str(int(line[6]) + 2))
                        start_c = line[5].lower()
                        stop_c = line[7].lower()
                    else:
                        url = "http://genome.ucsc.edu/cgi-bin/das/{0:s}/dna?segment={1:s}:{2:s},{3:s}".format(dataset, line[2], str(int(line[6]) - 1), str(int(line[4])))
                        start_c = uorfcommon.replace_seq(line[7]).lower()[::-1]
                        stop_c =  uorfcommon.replace_seq(line[5]).lower()[::-1]

                    temp_seq = url_load(url) 
                    
                    try:
                        temp_start_c = temp_seq[0:3]
                        temp_stop_c = temp_seq[-3::]
                        
                        if line[3] == '-':
                            temp_seq = uorfcommon.replace_seq(temp_seq[::-1])
                    except:
                        print line[1]
                        print url
                        print "aaa"

                    if temp_start_c == start_c and temp_stop_c == stop_c:
                        is_match = compare_all_seq(line[1], as_fna_dic[line[1]].lower(), temp_seq, write) #compare_all_seq(acc, symble, pos_abs_start, pos_abs_end, coden_start, coden_end, seq_5utr, seq_uorf)
                        
                        if not is_match:
                            is_match = compare_sub_seq(line[1], line[3], line[4], line[6], line[5].lower(), line[7].lower(), as_fna_dic[line[1]].lower(), write)
                    else:
                        is_match = compare_sub_seq(line[1], line[3], line[4], line[6], line[5].lower(), line[7].lower(), as_fna_dic[line[1]].lower(), write)
                    if not is_match:
                        if not line[0] in nomatch_dic:
                            nomatch_dic[line[0]] = [[line[1], line[3], line[4], line[5], line[6], line[7]]]
                        else:
                            nomatch_dic[line[0]].extend([[line[1], line[3], line[4], line[5], line[6], line[7]]])
                else:
                    if not line[0] in up1000_dic:
                        up1000_dic[line[0]] = [[line[1], line[3], line[4], line[5], line[6], line[7]]]
                    else:
                        up1000_dic[line[0]].extend([[line[1], line[3], line[4], line[5], line[6], line[7]]])
    fh.close()
    fw.close()

    fw = open("IRESdataset/{0:s}_five/csv/nomatch_up1000_abs_{0:s}.csv".format(species), 'w')
    write = csv.writer(fw)

    key = up1000_dic.keys()
    key.sort()

    for i in key:
        for j in up1000_dic[i]:
            temp_list = [i]
            temp_list.extend(j)
            write.writerows([temp_list])
    fw.close()

    fw = open("IRESdataset/{0:s}_five/csv/nomatch_parse_abs_{0:s}.csv".format(species), 'w')
    write = csv.writer(fw)

    key = nomatch_dic.keys()
    key.sort()

    for i in key:
        for j in nomatch_dic[i]:
            temp_list = [i]
            temp_list.extend(j)
            write.writerows([temp_list])
    fw.close()
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
def url_load(url):  #IREZone_size : at least IRES    
    try:
        url_xml = urllib2.urlopen(url).read()
        root = ET.fromstring(url_xml)

        temp_seq = re.sub('\s+', '', root[0][0].text)
        if temp_seq is None:
            url_load(url)
        elif (len(temp_seq) < 3):
            url_load(url)
        else:
            return temp_seq
    except:
        url_load(url)

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
def compare_all_seq(acc, seq_5utr, seq_uorf, write):
    len_uorf = len(seq_uorf)

    reg = "(?={0:s})".format(seq_uorf)
    index_seq_list = [m.start() for m in re.finditer(reg, seq_5utr)]
                
    # print index_seq_list           
    if (len(index_seq_list) == 1):
        rel_pos_start = index_seq_list[0] + 1
        rel_pos_end = index_seq_list[0] + len_uorf
        write.writerows([[acc, rel_pos_start, rel_pos_end]])            
        return True
    else:
        return False

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
def compare_sub_seq(acc, symble, pos_abs_start, pos_abs_end, coden_start, coden_end, seq_5utr, write):
    if symble == '+':
        pos_start = int(pos_abs_start)
        pos_end = int(pos_abs_end) + 2
    else:
        pos_start = int(pos_abs_end) - 1
        pos_end = int(pos_abs_start)

    len_uorf = pos_end - pos_start + 1 
    num_dot = pos_end - pos_start - 6 + 1

    if num_dot == -3:
        re_reg = coden_start
    else:
        re_reg = coden_start + ".{" + str(num_dot) + "}" + coden_end

    index_seq_list = [m.start() for m in re.finditer(re_reg, seq_5utr)]

    if (len(index_seq_list) == 1):
        rel_pos_start = index_seq_list[0] + 1
        rel_pos_end = index_seq_list[0] + len_uorf
        write.writerows([[acc, rel_pos_start, rel_pos_end]])
        return True
    else:
        return False

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
def compare_sub_rel_seq(acc, symble, pos_rel_start, pos_rel_end, coden_start, coden_end, seq_5utr, write):
    pos_start = int(pos_rel_start)
    pos_end = int(pos_rel_end)

    len_uorf = pos_end - pos_start + 1 
    num_dot = pos_end - pos_start - 6 + 1

    if num_dot == -3:
        re_reg = coden_start
    else:
        re_reg = coden_start + ".{" + str(num_dot) + "}" + coden_end

    index_seq_list = [m.start() for m in re.finditer(re_reg, seq_5utr)]
    if (len(index_seq_list) == 1):
        rel_pos_start = index_seq_list[0] + 1
        rel_pos_end = index_seq_list[0] + len_uorf
        write.writerows([[acc, rel_pos_start, rel_pos_end]])
        return True
    else:
        return False