#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import re, os, sys, csv
import sqlcommon
import MySQLdb

# ---------------- ---------------- ---------------- ---------------- #
# Global Variables
# ---------------- ---------------- ---------------- ---------------- #
db = MySQLdb.connect("localhost", "testuser", "testuser", "temp_hs")
cursor = db.cursor()

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
def ribosome_cut(species):
    coor_dic = sqlcommon.as_coor(species)
    fna_dic = sqlcommon.as_fna(coor_dic, species)

    if species == "hs":
        file_array = ['SRR970538.fdp', 'SRR970490.fdp', 'SRR970565.fdp', 'SRR970561.fdp', 'SRR970587.fdp', 'SRR970588.fdp', 'H_Rep1.fdp', 'H_Rep2.fdp', 'N_Rep1.fdp', 'N_Rep2.fdp']
    else:
        file_array = ['SRR3208406.fdp']
    f_ri_depth = "original_data/"
    
    ri_depth_dic = {}
    f_count = 0

    # Open and read the file
    for i in file_array:
        with open(f_ri_depth + i, 'r') as fh:
            data = csv.reader(fh)
    
            for line in data:
                print i
                print line[0]
                temp_array = line[0].split('\t')
                temp_depth_array = [temp_array[1]]

                if temp_array[1] != 'NA':
                    temp_ri_len = len(fna_dic[temp_array[0]]) + 1
                    for depth in line[1:temp_ri_len]:
                        temp_depth_array.append(depth)
                    temp_depth_str = ';'.join(temp_depth_array)
                    
                    count_depth = len(temp_depth_str) / 5000;
                    
                    for k in range(count_depth):
                        temp_key = temp_array[0] + "_" + str(k)
                        temp_split_str = temp_depth_str[k * 5000 : (k + 1) * 5000]
                        if not temp_key in ri_depth_dic:
                            ri_depth_dic[temp_key] = []
                        while f_count != len(ri_depth_dic[temp_key]):
                            ri_depth_dic[temp_key].append('')
                        ri_depth_dic[temp_key].append(temp_split_str)

                    temp_key = temp_array[0] + "_" + str(count_depth)
                    temp_split_str = temp_depth_str[count_depth * 5000::]
                    if not temp_key in ri_depth_dic:
                        ri_depth_dic[temp_key] = []

                    while f_count != len(ri_depth_dic[temp_key]):
                        ri_depth_dic[temp_key].append('')

                    ri_depth_dic[temp_key].append(temp_split_str)

                    if len(line) != len(fna_dic[temp_array[0]]):
                        print line
            f_count += 1
    fh.close()

    return ri_depth_dic

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
def sql_ri(species, ri_depth_dic):
    finalkey = ri_depth_dic.keys()
    finalkey.sort()
    acc_array = []
    accession_id = 0
    
    if species == 'hs':
        col_count = 11
    else:
        col_count = 2

    for j in finalkey:
        temp_dep_array = [j]
        temp_dep_array.extend(ri_depth_dic[j])
        while len(temp_dep_array) != col_count:
            temp_dep_array.append('')
        print len(temp_dep_array)
    
        temp_acc = temp_dep_array[0].split('_')[1]
        temp_acc = 'NM_' + temp_acc        

        if not temp_acc in acc_array:
            acc_array.append(temp_acc)
            accession_id += 1
            row_count = 1
            sql = """INSERT INTO saa_acc_id_%s(id, acc) VALUES (%d, '%s')""" % \
                    (species, accession_id, temp_acc)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()
        else:
            row_count += 1


        if species == 'hs':
            sql = """INSERT INTO saa_ri_dep_%s(acc_id, row_id, G1_1, G1_2, S1, S2, M1, M2, HRep1, HRep2, NRep1, NRep2) VALUES (%d, %d, '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')""" % \
                    (species, accession_id, row_count, temp_dep_array[1], temp_dep_array[2], temp_dep_array[3], temp_dep_array[4], temp_dep_array[5], temp_dep_array[6], temp_dep_array[7], temp_dep_array[8], temp_dep_array[9], temp_dep_array[10])
        else:
            sql = """INSERT INTO saa_ri_dep_%s(acc_id, row_id, SRR3208406) VALUES (%d, %d, '%s')""" % \
                    (species, accession_id, row_count, temp_dep_array[1])

            cursor.execute(sql)
            db.commit()
        except:
            db.rollback()
    db.close()

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
    #./exec_pat.py Homo\ sapiens\ jun\ 5UTR\ seq.txt gua
    species = sys.argv[1]

    ri_dep_dic = ribosome_cut(species)
    print "over"
    sql_ri(species, ri_dep_dic)        

    sys.exit(0)

