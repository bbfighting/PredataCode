#!/user/bin/env python

"""
This script is used to parse the stock information
"""

import re, os
#import subprocess

# ---------------- ---------------- ---------------- ---------------- #
# Global Variables
# ---------------- ---------------- ---------------- ---------------- #

class scan:

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
    def __init__(self, xmlname, sequence, species, startcoden, genename, method):
        self.gene = xmlname
        self.seq = sequence
        self.sc = startcoden
        self.spe = species
        self.method = method
        if genename == "":
            self.genename = xmlname
        else:
            self.genename = genename
        self.pat_num = 25

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
    def main(self):
#        path_scan = "/home/testuser/public_html/temp_web/other/tool_lib/scan_for_matches/scan_for_matches"
        path_scan = "/home/testuser/public_html/IRES_test/SingleIRES/library/scan_for_matches"

        xml_name = "IRESdataset/" + self.spe + "_" + self.method + "/xmlfile/" + self.gene + ".xml"

        total_len = self.scan_input()
        pat_count = total_len / self.pat_num + 1
        jpg_count = 0

        #create xml outfile
        with open(xml_name, 'w') as fh:
            fh.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
            fh.write("<IRESearch>\n")
            fh.write("\t<InputInfo>\n")
            fh.write("\t\t<GeneName>")
            fh.write(self.genename)
            fh.write("</GeneName>\n")
            fh.write("\t\t<Sequence>")
            fh.write(self.seq.replace('T', 'U'))
            fh.write("</Sequence>\n")
            fh.write("\t</InputInfo>\n")

            fh.write("\t<OutputInfo>\n")

            for i in range(pat_count):
                pattern = self.scan_pattern(i)
                output_str = ""

                cmd = path_scan + " -m1 pattern < input"
                output = os.popen(cmd)

                for line in output:
                    output_str += line
                
                if output_str != "":
                    IRES_write = self.IRES_info(output_str, i, pattern, jpg_count)
                    fh.write(IRES_write)
                    jpg_count += 1            

            fh.write("\t</OutputInfo>\n");
            fh.write("</IRESearch>");
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
    def scan_input(self):
        input_seq = self.sc + self.seq[::-1]

        # Open and read the file
        with open("input", 'w') as fh:
            fh.write(">" + self.gene + "\n")
            fh.write(input_seq)
        fh.close()

        return len(input_seq)

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
    def scan_pattern(self, count):        
        if count == 0:
            temp_pat = "{0:d}...{1:d}".format(3, 25 * (count + 1))
        else:
            temp_pat = "{0:d}...{1:d}".format(25 * count + 1, 25 * (count + 1))
        pat = "r1={au,ua,cg,gc,gu,ug} ^" + self.sc + " " + temp_pat + " p9=5...6 p10=3...8 r1~p9[1,0,0] 2...5 p1=5...6 0...6 p2=5...6 p8=0...5 p6=5...8 p7=3...5 r1~p6[1,0,0] p4=5...8 p5=3...8 r1~p4[1,0,0] p3=0...2 r1~p2[1,0,0] 0...6 r1~p1[1,0,0] (10...10|(9...9|(8...8|(7...7|(6...6|(5...5|(4...4|(3...3|(2...2|(1...1|0...0))))))))))"

        # Open and read the file
        with open("pattern", 'w') as fh:
            fh.write(pat)
        fh.close()
    
        return temp_pat

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
    def IRES_info(self, scan_output, count, pat, jpg_count):
        seq_len = len(self.seq)

        output_list = scan_output.split(' ')
        output_list.pop(-1)

        temp_structure = self.get_strcture(output_list[1:])

        output_start_end = re.search('(?<=\[)\d+,\d+', output_list[0])
        start_end = output_start_end.group(0).split(',')
        output_start = start_end[0]
        output_end = start_end[1]

        if count == 0:
            structure = temp_structure[0:len(temp_structure) - 3]
        else:
            structure = re.sub('(\.)+$', '..........', temp_structure)

        temp_start = seq_len - int(output_end) + 1 + 3
        temp_end = seq_len - int(output_start) + 1        

        diff_len = len(temp_structure) - len(structure)
        new_end = temp_end - temp_start + 1 - diff_len + 3

        xml_pat = pat
        xml_start = -(seq_len - temp_start + 1)
        xml_end = -(seq_len - (temp_end - diff_len + 3) + 1)
        xml_structure = structure
        xml_seq = self.seq.replace("T", "U")[temp_start - 1:temp_start + len(structure) - 1]
        xml_sc = temp_end - seq_len + 1

        xml_PSfile = "IRESdataset/{0:s}_{1:s}/RNAplot/{2:s}/{3:d}_ss.ps".format(self.spe, self.method, self.gene, jpg_count)
        xml_JPGfile = "IRESdataset/{0:s}_{1:s}/RNAplot/{2:s}/{3:d}_ss.png".format(self.spe, self.method, self.gene, jpg_count)

        temp_write = "\t\t<IRES StartCoden=\"{0:d}\" start=\"{1:d}\" end=\"{2:d}\" IRES2AUG=\"{3:s}\" PSfile=\"{4:s}\" JPGfile=\"{5:s}\">\n".format(
            xml_sc, xml_start, xml_end, xml_pat, xml_PSfile, xml_JPGfile)
        temp_write += "\t\t\t<Sequence>"          
        temp_write += xml_seq              
        temp_write += "</Sequence>\n"
        temp_write += "\t\t\t<Structure>"
        temp_write += xml_structure
        temp_write += "</Structure>\n"
        temp_write +=  "\t\t</IRES>\n"
    
        return temp_write

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
    def get_strcture(self, scan_list):
        replace_list = [".", "(", ".", "(", ".", "(", ".", ")", "(", ".", ")", ".", ")", ".", ")", ".", "(", ".", ")", ".", "."]
        temp_structure = ""

        for i in range(21):
            for j in range(len(scan_list[20 - i])):
                temp_structure += replace_list[i]
        return temp_structure