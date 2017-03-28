#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import re, sys, csv, sqlcommon, MySQLdb


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
def data_IRES_IREZone(species, coor_dic):  #IREZone_size : at least IRES    
    filename = "IRESdataset/{0:s}_five/csv/{0:s}_five_IREZone.csv".format(species)
    IREZone_dic = {}
    gene_IREZone_list = []
    nm_IREZone_list = []

    with open(filename, 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            if line[3] != "":
                #len_IREZone = len(line[3].split(';'))
                IREZone_dic[line[0]] = line[3]
                if not line[0] in nm_IREZone_list:
                    nm_IREZone_list.append(line[0])
                if not coor_dic[line[0]][0] in gene_IREZone_list:
                    gene_IREZone_list.append(coor_dic[line[0]][0])
            else:
                IREZone_dic[line[0]] = ""
    fh.close()
   
    return IREZone_dic

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
def data_uorf_sort(coor_dic, species):  #IREZone_size : at least IRES    
    f_uorf = "IRESdataset/{0:s}_five/csv/pos_combine_sort.csv".format(species)

    uorf_dic = {}
    uorf_str_dic = {}
    uorf_start_dic = {}
    uorf_end_dic = {}

    # Open and read the file
    with open(f_uorf, 'r') as fh:
        data = csv.reader(fh)

        for line in data: #len:260 3 ~ 258, start = -(260 - 3 + 1) = -258; end = -(260 - 258 + 1)  = -3
            temp_start = -(coor_dic[line[0]][1] - int(line[1]) + 1)
            temp_end = -(coor_dic[line[0]][1] - int(line[2]) + 1)

            if not line[0] in uorf_end_dic:
                uorf_end_dic[line[0]] = [temp_end]
            else:
                uorf_end_dic[line[0]].append(temp_end)

            temp_start_key = line[0] + "_" + str(temp_end)

            if not temp_start_key in uorf_end_dic:
                uorf_start_dic[temp_start_key] = [temp_start]
            else:
                uorf_start_dic[temp_start_key].append(temp_start)            
    fh.close()

    for i in uorf_end_dic.keys():
        end_uniq_list = list(set(uorf_end_dic[i]))
        end_uniq_list = sorted(end_uniq_list, reverse=True)

        for j in end_uniq_list:
            temp_key = str(i) + "_" + str(j)
            uorf_start_list = sorted(uorf_start_dic[temp_key], reverse=True)

            for k in uorf_start_list:
                if not i in uorf_dic:
                    uorf_dic[i] = [str(k) + "~" + str(j)]
                else:
                    uorf_dic[i].append(str(k) + "~" + str(j))            

    nm_uorf_list = []
    gene_uorf_list = []
        
    for q in uorf_dic.keys():
        uorf_str_dic[q] = ';'.join(uorf_dic[q])
                
        if not q in nm_uorf_list:
            nm_uorf_list.append(q)
        if not coor_dic[q][0] in gene_uorf_list:
            gene_uorf_list.append(coor_dic[q][0])

    return uorf_dic        


# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./sql_main.py mouse

    species = sys.argv[1]
    gene_list = []    

    coor_dic = sqlcommon.as_coor(species)
    nm_list = []
    
    for nm in coor_dic.keys():
        if coor_dic[nm][1] == 0:
            del coor_dic[nm]
        else:
            if not coor_dic[nm][0] in gene_list:
                gene_list.append(coor_dic[nm][0])
            if not nm in nm_list:
                nm_list.append(coor_dic[nm][1])

    IREZone_dic = data_IRES_IREZone(species, coor_dic)
    uorf_dic = data_uorf_sort(coor_dic, species)
    all_data_dic = {}

    gene_IREZone_uorf_list = []
    nm_IREZone_uorf_list = []
    gene_IREZone_list = []
    uorf_all_count = 0
    uorf_overlap_count = 0

    for i in IREZone_dic.keys():
        is_uorf = False
        is_gene_IREZone = False
        is_nm_IREZone = False
        sub_IREZone = IREZone_dic[i].split(';') 
        temp_uorf = []
        count = 0
        if IREZone_dic[i] != "":
            all_data_dic[i] = [len(sub_IREZone)]
            is_nm_IREZone = True
            if not coor_dic[i][0] in gene_IREZone_list:
                gene_IREZone_list.append(coor_dic[i][0])
                is_gene_IREZone = True
        else:
            all_data_dic[i] = [0]
        for j in sub_IREZone:
            if i in uorf_dic:
                if is_nm_IREZone:
                    nm_IREZone_uorf_list.append(i)
                if is_gene_IREZone:
                    gene_IREZone_uorf_list.append(coor_dic[i][0])
                is_uorf = True
                for k in uorf_dic[i]:
                    if not k in temp_uorf:
                        if j != "":                        
                            IREZone_start = int(j.split(':')[0])
                            IREZone_end = int(j.split(':')[1])

                            uorf_start = int(k.split('~')[0])
                            uorf_end = int(k.split('~')[1])

                            if ((uorf_end < IREZone_end) and (uorf_end > IREZone_start)) or ((uorf_start < IREZone_end) and (uorf_start > IREZone_start)):
                                    
                                temp_uorf.append(k)
                                count += 1
        if is_uorf:
            all_data_dic[i].append(len(uorf_dic[i]))
            uorf_all_count += len(uorf_dic[i])
        else:
            all_data_dic[i].append(0)
        if is_uorf and IREZone_dic[i] != "":
            all_data_dic[i].append(count)
            all_data_dic[i].append(len(uorf_dic[i]) - count)
            all_data_dic[i].append(float(count) / float(len(uorf_dic[i])) * 100)
            all_data_dic[i].append((float(len(uorf_dic[i]) - count) / float(len(uorf_dic[i]))) * 100)
            uorf_overlap_count += (count)
        else:         
            all_data_dic[i].append(0)
            all_data_dic[i].append(0)
            all_data_dic[i].append(0)
            all_data_dic[i].append(0)
            all_data_dic[i].append(0)

    filename = "IRESdataset/{0:s}_five/csv/IREZone_uorf_num_{0:s}.csv".format(species)

    with open(filename, 'w') as fw:
        for k in all_data_dic.keys():
            write = csv.writer(fw)
            temp_list = [k] 
            temp_list.extend(all_data_dic[k])
                
            write.writerows([temp_list])
    fw.close()
    db.close()

    print "uorf_all_count"
    print uorf_all_count
    print "uorf_overlap_count"
    print uorf_overlap_count
                            
    sys.exit(0)