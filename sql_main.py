#!/usr/bin/env python

"""
This script is used to parse the stock information
"""

import re, sys, csv, MySQLdb
import sqlcommon
import urllib2
from bs4 import BeautifulSoup

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
def sql_as(fna_dic, coor_dic, g_list, species):  #IREZone_size : at least IRES    

    key = coor_dic.keys()
    key.sort()

    for ac in key:
        if len(fna_dic[ac]) != 0:
            sql = """INSERT INTO saa_as_%s(accession_number, gene, seq_len) VALUES ('%s', '%s', '%s')""" % \
                    (species, ac, coor_dic[ac][0], len(fna_dic[ac]))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()
    
            sql = """INSERT INTO saa_IRES_IREZone_uorf_%s(accession_number, seq_len) VALUES ('%s', %d)""" % \
                    (species, ac, len(fna_dic[ac]))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()
    
            sql = """INSERT INTO saa_5utr_%s(ass, seq_5utr) VALUES ('%s', '%s')""" % \
                    (species, ac, fna_dic[ac])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()

    for gene in g_list:
        sql = """INSERT INTO saa_gene_%s(official_symbol) VALUES ('%s')""" % \
                (species, gene)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            db.rollback() 

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
def sql_entrezNCBI(gene_list, species):  #IREZone_size : at least IRES    
    if species == 'hs':
        tax_id = '9606'
    elif species == 'mouse':
        tax_id = '10090'

    #filename = "original_data/human_NCBI_160112.tsv"
    filename = "original_data/All_Data.gene_info"
    geneid_list = []

    # Open and read the file
    with open(filename, 'r') as fh:
        fh.next()

        for line in fh:
            detail = re.split('\t', line)
            if detail[0] == tax_id:
                if detail[2] in gene_list:
                    detail[5] = detail[5].replace("|", ", ")
                    detail[5] = detail[5].replace('\'', '\\\'')
                    detail[4] = detail[4].replace("|", ", ")
                    detail[4] = detail[4].replace('\'', '\\\'')
                    detail[8] = detail[8].replace('\'', '\\\'')
                    geneid_list.append(detail[1])
                    
                    sql = """UPDATE saa_gene_%s SET gene_id = '%s', fullname = '%s', related = '%s', known = '%s' WHERE official_symbol = '%s' """ % \
                        (species, detail[1], detail[8], detail[5], detail[4], detail[2])
                    try:
                        cursor.execute(sql)
                        db.commit()
                    except:
                        db.rollback()

    fh.close()
    return geneid_list

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
def sql_IRES_IREZone(species):  #IREZone_size : at least IRES    
    filename = "IRESdataset/{0:s}_five/csv/{0:s}_five_IREZone.csv".format(species)

    with open(filename, 'r') as fh:
        data = csv.reader(fh)
        data.next()

        for line in data:
            if line[2] != "":
                len_IRES = len(line[2].split(';'))
            else:
                len_IRES = 0

            if line[3] != "":
                len_IREZone = len(line[3].split(';'))
            else:
                len_IREZone = 0
            
            sql = """UPDATE saa_IRES_IREZone_uorf_%s SET IRES = '%s', IREZone = '%s' WHERE accession_number = '%s' """ % \
                    (species, line[2], line[3], line[0])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()

            sql = """UPDATE saa_as_%s SET numIRES = %d, numIREZone = %d WHERE accession_number = '%s' """ % \
                    (species, len_IRES, len_IREZone, line[0])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                db.rollback()
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
def sql_uorf_sort(coor_dic, species):  #IREZone_size : at least IRES    
    f_uorf = "IRESdataset/{0:s}_five/csv/pos_combine_sort_{0:s}.csv".format(species)

    uorf_dic = {}
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
        
    for q in uorf_dic.keys():

        sql = """UPDATE saa_IRES_IREZone_uorf_%s SET uorf = '%s' WHERE accession_number = '%s' """ % \
                (species, ';'.join(uorf_dic[q]), q)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            db.rollback()

        sql = """UPDATE saa_as_%s SET numuorf = %d WHERE accession_number = '%s' """ % \
                (species, len(uorf_dic[q]), q)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            db.rollback()

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
def sql_graphic(geneid_list, species, count, load_count): #*2 : "TEC-7006", "MEMO1-7795", "SFPQ-6421", "HBD-3045", "MMD2-221938, 100505381", "LRRC37A-9884"
    for gid in geneid_list[count::]:
        url = 'http://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=ShowDetailView&TermToSearch={0:s}'.format(gid)
        try:
            soup = BeautifulSoup(urllib2.urlopen(url).read())
            title_search = soup.find(title="Nucleotide Graphic report")

            if (title_search):
                graph_href = title_search['href']
            else:
                graph_href = "no"
    
            if soup.find(id="summaryDl").find_all('dt')[-2].string == 'Summary':
                summary = soup.find(id="summaryDl").find_all('dd')[-2].string.replace('\'', '\\\'')
            else:
                summary = ""

            sql = """UPDATE saa_gene_%s SET summary = '%s', graphic = '%s' WHERE gene_id = '%s' """ % \
                    (species, summary, graph_href, gid)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                print "sql_wrong"
                print gid
            count += 1
            load_count = 0
            print gid
            print "success"
        except:
            if (load_count < 3):
                load_count += 1
                sql_graphic(geneid_list, species, count, load_count)
            else:
                print "load_wrong"
                print gid

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
def sql_getgeneid(gene_list, species):  #IREZone_size : at least IRES    
    if species == 'hs':
        tax_id = '9606'
    elif species == 'mouse':
        tax_id = '10090'

    filename = "original_data/All_Data.gene_info"
    geneid_list = []

    # Open and read the file
    with open(filename, 'r') as fh:
        fh.next()

        for line in fh:
            detail = re.split('\t', line)
            if detail[0] == tax_id:
                if detail[2] in gene_list:
                    detail[5] = detail[5].replace("|", ", ")
                    detail[4] = detail[4].replace("|", ", ")
                    geneid_list.append(detail[1])

    fh.close()
    return geneid_list
# ---------------- ---------------- ---------------- ---------------- #
# Main Program
# ---------------- ---------------- ---------------- ---------------- #
if __name__ == "__main__":
#./sql_main.py mouse

    species = sys.argv[1]

    filename = "IRESdataset/{0:s}_five/csv/saa_gene_{0:s}.csv".format(species)
    geneid_list = []

    # Open and read the file
    with open(filename, 'r') as fh:
        #fh.next()
        data = csv.reader(fh)

        for line in data:
            geneid_list.append((str(line[0])))
    print len(geneid_list)
    print "geneid_list finish"
    sql_graphic(geneid_list, species, 0, 0)
    print "sql_graphic finish"    

    db.close()
    sys.exit(0)
