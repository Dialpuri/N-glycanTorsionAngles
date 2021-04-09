import privateer
import Bio
from Bio.PDB import PDBList
import json
import csv 

pdbl = PDBList()
output_data = []
residue_data_list = []

#CHOOSE OUTPUT FORMAT.
output_format = 'csv'
#output_format = 'json'

def sugar_sugar_read(): 
    #Accepts input from console for file name.
    file_name = raw_input("What is the name of the text file (enter in format; NAME.txt) ")
    sug_sug_output_file = "1.1_output_{}.txt".format(output_format) #If you wanted to change the name of the output file do it here.

    #INITAL VARIABLES - NO NEED TO EDIT.
    sugar_name = ""
    sugar_chain = ""
    sugar_number = ""
    pdb_code = ""
    rscc_1 = ""

    #Opens your input file - NO NEED TO EDIT ANY OF THIS.
    with open(file_name, "r") as myFile:
        for lineNo, line in enumerate(myFile):
            if lineNo != 0:
                split_line = line.split(',')
                if split_line[0] == pdb_code:
                    sugar_chain_position = split_line[1].split('-')
                    sugar_chain_position[2].replace(" ", "")
                    res = split_line[2]
                    q = split_line[3]
                    phi = split_line[4]
                    theta = split_line[5]
                    rscc = split_line[6]
                    detectedType = split_line[7]
                    cnf = split_line[8]
                    mFo = split_line[9]
                    Bfac = split_line[10]

                    if sugar_chain_position[1] == sugar_chain:
                        if int(sugar_chain_position[2]) == (int(sugar_number) + 1):
                            print("PDB CODE: {}, Sugars: {} and {}, Chain: {}, Numbers:{} and {}".format(pdb_code, sugar_name, sugar_chain_position[0], sugar_chain_position[1], sugar_number, sugar_chain_position[2]))
                            pdbl.retrieve_pdb_file(pdb_code ,pdir='.',file_format="pdb")
                            temp_pdb_code = "pdb" + pdb_code + ".ent"

                            x = privateer.print_glycosidic_torsions(temp_pdb_code, sugar_name, sugar_chain_position[0])
                            
                            for entry in x: 
                                split_sugar_data_1 = entry[0].split('/')
                                split_sugar_data_2 = entry[1].split('/')
                                sugar_chain_1 = split_sugar_data_1[1]
                                sugar_chain_2 = split_sugar_data_2[1]
                                position_sugar_1 =  split_sugar_data_1[2].split('(')
                                position_sugar_2 = split_sugar_data_2[2].split('(')
                                position_1 = position_sugar_1[0]
                                position_2 = position_sugar_2[0]
                                sugar_name_1 = position_sugar_1[1]
                                sugar_name_2 = position_sugar_2[1]
                                sugar_name_1 = sugar_name_1[:-1]
                                sugar_name_2 = sugar_name_2[:-1]

                                if sugar_chain_2 == sugar_chain_position[1] and sugar_chain_1 == sugar_chain: 
                                    if sugar_name == sugar_name_1 and sugar_chain_position[0] == sugar_name_2:
                                        if int(position_1) == int(sugar_number) and int(sugar_chain_position[2]) == int(position_2):
                                            if output_format == 'json':
                                                sugar_data_json = {
                                                    "PDB": pdb_code,
                                                    "Sugar 1": sugar_name_1,
                                                    "Sugar 2": sugar_name_2,
                                                    "Chain" : sugar_chain,
                                                    "Position 1": position_1,
                                                    "Position 2": position_2,
                                                    "Phi": entry[2],
                                                    "Psi": entry[3],
                                                    "Resolution": res,
                                                    "Cremer-Pope": {
                                                        "Q" : q,
                                                        "Phi": phi,
                                                        "Theta": theta
                                                    },
                                                    "Sugar 1 RSCC": rscc_1,
                                                    "Sugar 2 RSCC": rscc,
                                                    "Detected Type": detectedType,
                                                    "CNF": cnf,
                                                    "mFo": mFo,
                                                    "Bfac": Bfac
                                                }
                                                output_data.append(sugar_data_json)
                                            elif output_format == 'csv':
                                                output_data.append([pdb_code,sugar_name_1,sugar_name_2,sugar_chain,position_1,position_2,entry[2],entry[3],res,q,phi,theta,rscc_1,rscc,detectedType,cnf,mFo,Bfac])
                                            else:
                                                print("The output format has not been detected correctly, check that the output_file = csv or output_file = json")

                    sugar_name = sugar_chain_position[0]
                    sugar_chain = sugar_chain_position[1]
                    sugar_number = sugar_chain_position[2].strip()
                    rscc_1 = rscc

                pdb_code = split_line[0] 

            elif lineNo == 0:
                if output_format == 'csv':
                    output_data.append(['PDB','Sugar 1','Sugar 2','Chain','Position 1','Position 2','Phi','Psi','Resolution','Q','Phi','Theta','Sugar 1 RSCC','Sugar 2 RSCC','Detected Type','CNF','mFo','Bfac'])

    with open(sug_sug_output_file, "wb") as output:
        if output_format == 'json':
            json.dump(output_data, output, indent=4)
        elif output_format == 'csv':
            writer = csv.writer(output, delimiter=',')
            writer.writerows(output_data)

def residue_sugar_read():

    #Accepts input from console for file name.
    file_name = raw_input("What is the name of the text file (enter in format; NAME.txt) ")
    res_sug_output_file = "Residue And Sugar Output.txt" #If you wanted to change the name of the output file do it here
    pdb_code = ""
    
    #Opens your input file - NO NEED TO EDIT ANY OF THIS.
    with open(file_name, "r") as residue_read:
        for lineNo, line in enumerate(residue_read):
            if lineNo != 0:
                split_line = line.split(',')
                sugar_chain_position = split_line[1].split('-')
                sugar_chain_position[2].replace(" ", "")
                pdb_code = split_line[0]
                res = split_line[2]
                q = split_line[3]
                phi = split_line[4]
                theta = split_line[5]
                rscc = split_line[6]
                detectedType = split_line[7]
                cnf = split_line[8]
                mFo = split_line[9]
                Bfac = split_line[10]

                pdbl.retrieve_pdb_file(pdb_code ,pdir='.',file_format="pdb")
                temp_pdb_code = "pdb" + pdb_code + ".ent"

                y = privateer.print_glycosidic_torsions(temp_pdb_code)

                for entry in y:
                    residue = entry[0].split('/')
                    sugar_data = entry[1].split('/')
                    residue_chain = residue[1]
                    sugar_chain = sugar_data[1]
                    position_AA = residue[2].split('(')
                    AA_position = position_AA[0]
                    AA = position_AA[1]
                    AA = AA[:-1]

                    position_sugar = sugar_data[2].split('(')
                    sugar_position = position_sugar[0]
                    sugar = position_sugar[1]
                    sugar = sugar[:-1]

                    if sugar == sugar_chain_position[0] and sugar_chain == sugar_chain_position[1] and int(sugar_position) == int(sugar_chain_position[2]):
                        if output_format == 'json':
                            residue_data_json = {
                                "PDB": pdb_code,                            
                                "Phi": entry[2],
                                "Psi": entry[3],
                                "Resolution": res,
                                "Cremer-Pope": {
                                    "Q" : q,
                                    "Phi": phi,
                                    "Theta": theta
                                },
                                "Residue Data": {
                                    "Residue Chain": residue_chain,
                                    "Residue Position": AA_position,
                                    "Amino Acid": AA
                                },
                                "Sugar Data": {
                                    "Sugar Chain": sugar_chain,
                                    "Sugar Position": sugar_position, 
                                    "Sugar": sugar                               
                                },
                                "Sugar RSCC": rscc,
                                "Detected Type": detectedType,
                                "CNF": cnf,
                                "mFo": mFo,
                                "Bfac": Bfac
                            }
                            residue_data_list.append(residue_data_json)          
                        elif output_format =='csv':
                            residue_data_list.append([pdb_code,residue_chain, AA_position, AA, sugar_chain, sugar_position, sugar, entry[2], entry[3], rscc, q, phi, theta, detectedType, cnf, mFo, Bfac])
                        else:
                            print("The output format has not been detected correctly, check that the output_file = csv or output_file = json")
            elif lineNo == 0:
                if output_format == 'csv':
                    residue_data_list.append(['PDB','Residue Chain','Residue Position','Amino Acid','Sugar Chain','Sugar Position','Sugar','Phi','Psi','Sugar RSCC','Q','Phi','Theta','Detected Type','CNF','mFo','Bfac'])

    with open(res_sug_output_file, "wb") as output:
        if output_format == 'json':
            json.dump(residue_data_list, output, indent=4)
        elif output_format == 'csv':
            writer = csv.writer(output, delimiter=',')
            writer.writerows(residue_data_list)


#CHOOSE WHICH OUTPUT YOU WANT
#residue_sugar_read()
sugar_sugar_read()


