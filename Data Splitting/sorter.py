import os

file_name = 'filtered_PDB_RSCC_0.8.txt'
header = "PDB,Sugar,Rsln,Q,Phi,Theta,RSCC,Detected type,Cnf,<mFo>,<Bfac>,Ctx,Ok?"
count = 0

with open(file_name) as file:
    for line in file:
        count += 1

numberOfLines = count - 1
fileNumber = int(input("How many files to output? "))
linesPerFile = round(numberOfLines / fileNumber)
splitFile = None
fileNo = 0

with open(file_name) as base_file:
    for lineNumber, line in enumerate(base_file):
        if lineNumber % linesPerFile == 0:
            fileNo += 1
            if splitFile is not None:
                splitFile.close()
            splitFileName = 'PDB_split_no_{}.txt'.format(fileNo)
            splitFile = open(splitFileName, "w")
            if lineNumber != 0:
                splitFile.write(header + "\n")
        splitFile.write(line)
    if splitFile:
        splitFile.close()
