import os

file_name = input("What is the name of the text file (enter in format; NAME.txt) ")

sugar_name = ""
sugar_chain = ""
sugar_number = ""
pdb_code = ""

with open(file_name, "r") as myFile:
    for lineNo, line in enumerate(myFile):
        if lineNo != 0:
            split_line = line.split(',')
            if split_line[0] == pdb_code:
                sugar_chain_position = split_line[1].split('-')
                sugar_chain_position[2].replace(" ", "")

                #if sugar_chain_position[0] == sugar_name:
                if sugar_chain_position[1] == sugar_chain:
                    if int(sugar_chain_position[2]) == (int(sugar_number) + 1):
                        print(f"Potential Joined Chains found \n PDB CODE: {pdb_code}, Sugars: {sugar_name} and {sugar_chain_position[0]}, Chain: {sugar_chain_position[1]}, Numbers:{sugar_number} and {sugar_chain_position[2]}")

                sugar_name = sugar_chain_position[0]
                sugar_chain = sugar_chain_position[1]
                sugar_number = sugar_chain_position[2].strip()
            
            pdb_code = split_line[0] 
