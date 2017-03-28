#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import csv, sys, re
import IRES, parse_xml

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
def cds35(info_name):
    gene_dic = {}
    
    # Open and read the file
    with open(info_name, 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            seperate_list = [1]
            data_list = re.split('\s+', line[0])
            seq_list = data_list[2].split('..')
            cds_list = data_list[3].split('..')
            seperate_list.extend([int(cds_list[0]) - 1, int(cds_list[0]), int(cds_list[1]), int(cds_list[1]) + 1, int(seq_list[1])])
            
            gene_dic[data_list[0]] = [seperate_list]
            gene_dic[data_list[0]].append(data_list[1])
    fh.close()
    return gene_dic
    

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
def pre_content(filename, gene_dic, species, method, sc = ""):
    #print gene_dic

    # Open and read the file
    fh = open(filename, 'r')
    fcontent = fh.read()
    sub_content_list = fcontent.split('>')
    sub_content_list.pop(0)

    if method == "cds":
        index_p = 3
    elif method == "three":
        index_p = 5

    for item in sub_content_list:
        temp_item = re.sub('(\s+)$', '', item)
        temp_item = temp_item.replace('\n', '  ')
        gene_seq_list = re.split('\s{2,}', temp_item)
        gene = gene_seq_list[0]
        seq = "".join(gene_seq_list[1:])

        if gene in gene_dic:
            genename = gene_dic[gene][1]
        else:
            genename = gene

        if method == "web":
            gene = re.sub('\s.+', '', gene)
            IRES.scan(gene, seq, species, sc, genename, method).main()
        elif method == "five":
            UTR5(gene, seq, gene_dic[gene][0][1], species, genename, method)
        else: # method == "cds" || method == "three"
            CDS_UTR3(gene, seq, gene_dic[gene][0], index_p, species, genename, method)

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
def UTR5(gene, sequence, position, species, genename, method):
    if position != 0:
        UTR5_seq = sequence[0:position]
        temp_sc = sequence[position:position + 3][::-1].lower()
        sc = temp_sc.replace("t", "u")
        IRES.scan(gene, UTR5_seq, species, sc, genename, method).main()

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
def CDS_UTR3(gene, sequence, position, index_p, species, genename, method):
    if position[index_p] > position[index_p - 1]: 
        temp_seq = sequence[position[index_p - 2]:position[index_p]]
        print gene
        print temp_seq + "\n"
        print_sequence = temp_seq #print
        last_ATG = temp_seq.rfind("ATG")
        sc = "gua"

        while not (check(temp_seq[last_ATG + 3:], gene) or last_ATG == -1): 
            temp_seq = temp_seq[0:last_ATG]
            last_ATG = temp_seq.rfind("ATG")
        print "\nnew_seq-----------------------\n", temp_seq[0:last_ATG]

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
def check(seq, gene):
    print "\n", seq
    len_seq = len(seq)

    if not len_seq < 45:
        for i in range(15):
            temp_seq =  seq[i * 3:(i + 1) * 3]
            print temp_seq            

            if temp_seq == 'TAA' or temp_seq == 'TAG' or temp_seq == 'TGA': 
                print "break"
                return False
        
        return True
    else:
        return False

# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
    #./exec_pat.py Homo\ sapiens\ jun\ 5UTR\ seq.txt gua

    upload_file = "original_data/" + sys.argv[1] #ex : original_data/Homo sapiens jun 5UTR seq.txt
    sc = sys.argv[2] #ex : gua
    pre_content(upload_file, {}, "predict", "web", sc)

    sys.exit(0)

