import os

file_name = 'filtered_PDB_RSCC.txt'

filterFile = None

with open(file_name) as base_file:
    for lineNo, line in enumerate(base_file):
        if lineNo != 0:
            print(lineNo)
            temp = line.split(',')
            if temp[6] == "":
                continue
            if float(temp[6]) > 0.8:
                if filterFile is not None:
                    filterFile.close()
                filterFile = open("filtered_PDB_RSCC_0.8.txt",'a')
                filterFile.write(line)
    if filterFile:
        filterFile.close()
