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
def Align_print(header_list, body_list, diff = 5):    
    temp_print = ""

    for i in header_list:
        if i == header_list[0]:
            temp_print += "%{0:d}s".format(len(i) + diff)
        else:
            temp_print += "%{0:d}s".format(len(i) + diff + 5)

    print temp_print % tuple(header_list) 

    len_col = len(header_list)
    len_row = len(body_list) / len_col

    for j in range(len_row):
        print temp_print % tuple(body_list[j * len_col: (j + 1) * len_col]) + "\n"