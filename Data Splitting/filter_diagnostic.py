import os

file_name = 'privateerPDB_RSCC.txt'

filterFile = None

with open(file_name) as base_file:
    for line in base_file:
        temp = line.split(',')
        if temp[12] == "yes\n":
            if filterFile is not None:
                filterFile.close()
            filterFile = open("filtered_PDB_RSCC.txt",'a')
            filterFile.write(line)
    if filterFile:
        filterFile.close()
