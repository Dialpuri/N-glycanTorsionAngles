import privateer
import Bio
from Bio.PDB import PDBList
import json
import csv 
import pandas as pd

pdbl = PDBList()
output_data = []
residue_data_list = []

#CHOOSE OUTPUT FORMAT.
output_format = 'csv'
#output_format = 'json'

#IF YOU WANT DON'T WANT ALL THE EXTRA DATA THAT WAS IN THE INPUT FILE SET reduced_output = True
reduced_output = False 

def sugar_sugar_read(): 
    #Accepts input from console for file name.
    file_name = raw_input("What is the name of the text file (enter in format; NAME.txt) ")
    if reduced_output == False:
        sug_sug_output_file = "ALL_sug_sug_output_{}.txt".format(output_format) #If you wanted to change the name of the output file do it here.
    else:
        sug_sug_output_file = "reduced_sug_sug_output_{}.txt".format(output_format)

    #INITAL VARIABLES - NO NEED TO EDIT.
    pdb_code = ""
    rscc_1 = ""

    data = pd.read_csv(file_name,delimiter=',')
    data['Sugar'] = data['Sugar'].apply(lambda x: x.strip())
    grouped_data = data.groupby('PDB')['Sugar'].apply(list)
    
    pdb_checklist = []
    pdb_completed = []

    for name, group in grouped_data.items():
        pdb = {"PDB":name, "Sugars": group}
        pdb_checklist.append(pdb)
      
    #Opens your input file - NO NEED TO EDIT ANY OF THIS.
    with open(file_name, "r") as myFile:
        for lineNo, line in enumerate(myFile):
            if lineNo != 0:
                split_line = line.split(',')
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
                pdb_code = split_line[0]
                if pdb_code not in pdb_completed:

                    pdbl.retrieve_pdb_file(pdb_code, pdir='.',file_format = "pdb")
                    temp_pdb_code = "pdb" + pdb_code + ".ent"

                    branched_glycans = privateer.get_branched_glycans(temp_pdb_code)

                    for entry in branched_glycans:
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

                        concated_sugar_1 = sugar_name_1 + "-" + sugar_chain_1 + "-" + position_1
                        concated_sugar_2 = sugar_name_2 + "-" + sugar_chain_2 + "-" + position_2

                        for value in pdb_checklist:
                            if value['PDB'] == pdb_code:
                                if concated_sugar_1 in value['Sugars']:
                                    if concated_sugar_2 in value['Sugars']:
                                        sug_1 = data.loc[(data['PDB']==pdb_code) & (data['Sugar'] == concated_sugar_1)]
                                        sug_2 = data.loc[(data['PDB']==pdb_code) & (data['Sugar'] == concated_sugar_2)]
                                        if output_format == 'json':
                                            if reduced_output == False:
                                                sugar_data_json = {
                                                    "PDB": pdb_code,
                                                    "Sugar 1": sugar_name_1,
                                                    "Sugar 2": sugar_name_2,
                                                    "Chain" : sugar_chain_1,
                                                    "Position 1": position_1,
                                                    "Position 2": position_2,
                                                    "Phi": entry[2],
                                                    "Psi": entry[3],
                                                    "Sug 1 Rsln": float(sug_1['Rsln'].values),
                                                    "Sug 2 Rsln": float(sug_2['Rsln'].values),
                                                    "Sug 1 Q":float(sug_1['Q'].values),
                                                    "Sug 2 Q":float(sug_2['Q'].values),
                                                    "Sug 1 Phi":float(sug_1['Phi'].values),
                                                    "Sug 2 Phi":float(sug_2['Phi'].values),
                                                    "Sug 1 Theta":float(sug_1['Theta'].values),
                                                    "Sug 2 Theta":float(sug_2['Theta'].values),
                                                    "Sug 1 RSCC":float(sug_1['RSCC'].values),
                                                    "Sug 2 RSCC":float(sug_2['RSCC'].values),
                                                    "Sug 1 Detected Type":str(sug_1['Detected type'].item()),
                                                    "Sug 2 Detected Type":str(sug_2['Detected type'].item()),
                                                    "Sug 1 Cnf":str(sug_1['Cnf'].item()),
                                                    "Sug 2 Cnf":str(sug_2['Cnf'].item()),
                                                    "Sug 2 mFo":float(sug_1['<mFo>'].values),
                                                    "Sug 2 mFo":float(sug_2['<mFo>'].values),
                                                    "Sug 1 Bfac":float(sug_1['<Bfac>'].values),
                                                    "Sug 2 Bfac":float(sug_2['<Bfac>'].values),
                                                    "Sug 1 Ctx":str(sug_1['Ctx'].item()),
                                                    "Sug 2 Ctx":str(sug_2['Ctx'].item()),
                                                    "Sug 1 Ok":str(sug_1['Ok?'].item()),
                                                    "Sug 2 Ok":str(sug_2['Ok?'].item())
                                                }
                                                output_data.append(sugar_data_json)
                                            else:
                                                sugar_data_json = {
                                                    "PDB": pdb_code,
                                                    "Sugar 1": sugar_name_1,
                                                    "Sugar 2": sugar_name_2,
                                                    "Chain" : sugar_chain_1,
                                                    "Position 1": position_1,
                                                    "Position 2": position_2,
                                                    "Phi": entry[2],
                                                    "Psi": entry[3]
                                                }
                                                output_data.append(sugar_data_json)

                                        elif output_format == 'csv':
                                            if reduced_output == False:
                                                output_data.append([pdb_code,sugar_name_1,sugar_name_2,sugar_chain_1,position_1,position_2,entry[2],entry[3],float(sug_1['Rsln'].values),float(sug_2['Rsln'].values),float(sug_1['Q'].values),float(sug_2['Q'].values),float(sug_1['Phi'].values),float(sug_2['Phi'].values),float(sug_1['Theta'].values),float(sug_2['Theta'].values),float(sug_1['RSCC'].values),float(sug_2['RSCC'].values),str(sug_1['Detected type'].item()),str(sug_2['Detected type'].item()),str(sug_1['Cnf'].item()),str(sug_2['Cnf'].item()),float(sug_1['<mFo>'].values),float(sug_2['<mFo>'].values),float(sug_1['<Bfac>'].values),float(sug_2['<Bfac>'].values),str(sug_1['Ctx'].item()),str(sug_2['Ctx'].item()),str(sug_1['Ok?'].item()),str(sug_2['Ok?'].item())])
                                            else:
                                                output_data.append([pdb_code,sugar_name_1,sugar_name_2,sugar_chain_1,position_1,position_2,entry[2],entry[3]])
                                        else:
                                            print("The output format has not been detected correctly, check that the output_file = csv or output_file = json")
                 
                    rscc_1 = rscc
                    pdb_completed.append(pdb_code)

            elif lineNo == 0:
                if output_format == 'csv':
                    if reduced_output == False:
                        output_data.append(['PDB','Sugar 1','Sugar 2','Chain','Position 1','Position 2','Phi','Psi','Sugar 1 Resolution', 'Sugar 2 Resolution','Sugar 1 Q','Sugar 2 Q','Sugar 1 Phi','Sugar 2 Phi','Sugar 1 Theta','Sugar 2 Theta','Sugar 1 RSCC','Sugar 2 RSCC','Sugar 1 Detected Type','Sugar 2 Detected Type','Sugar 1 CNF','Sugar 2 CNF','Sugar 1 mFo','Sugar 2 mFo','Sugar 1 Bfac','Sugar 2 Bfac',"Sugar 1 Ctx", "Sugar 2 Ctx", "Sugar 1 Ok?", "Sugar 2 Ok?"])
                    else:
                        output_data.append(['PDB','Sugar 1','Sugar 2','Chain','Position 1','Position 2','Phi','Psi'])

    with open(sug_sug_output_file, "wb") as output:
        if output_format == 'json':
            json.dump(output_data, output, indent=4)
        elif output_format == 'csv':
            writer = csv.writer(output, delimiter=',')
            writer.writerows(output_data)

def residue_sugar_read():

    #Accepts input from console for file name.
    file_name = raw_input("What is the name of the text file (enter in format; NAME.txt) ")
    
    if reduced_output == False:
        res_sug_output_file = "res_sug_output_{}.txt".format(output_format) #If you wanted to change the name of the output file do it here
    else:
        res_sug_output_file = "reduced_res_sug_output_{}.txt".format(output_format)

    pdb_code = ""
    
    data = pd.read_csv(file_name,delimiter=',')
    data['Sugar'] = data['Sugar'].apply(lambda x: x.strip())
    grouped_data = data.groupby('PDB')['Sugar'].apply(list)
    pdb_checklist = []
    pdb_completed = []

    for name, group in grouped_data.items():
        stripped_group = [sug.strip() for sug in group]
        pdb = {"PDB":name, "Sugars": stripped_group}
        pdb_checklist.append(pdb)

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

                if pdb_code not in pdb_completed:
                    y = privateer.get_residue_sugar(temp_pdb_code)

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

                        concated_sugar = sugar + "-" + sugar_chain + "-" + sugar_position

                        for value in pdb_checklist:
                            if value['PDB'] == pdb_code:
                                if concated_sugar in value['Sugars']:
                                    sug = data.loc[(data['PDB']==pdb_code) & (data['Sugar'] == concated_sugar)]
                                    if output_format == 'json':
                                        if reduced_output == False:
                                            residue_data_json = {
                                                "PDB": pdb_code,                            
                                                "Phi": entry[2],
                                                "Psi": entry[3],
                                                "Resolution": float(sug['Rsln'].values),
                                                "Cremer-Pope": {
                                                    "Q" : float(sug['Q'].values),
                                                    "Phi": float(sug['Phi'].values),
                                                    "Theta": float(sug['Theta'].values)
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
                                                "Sugar RSCC": float(sug['RSCC'].values),
                                                "Detected Type": str(sug['Detected type'].item()),
                                                "CNF": str(sug['Cnf'].item()),
                                                "mFo": str(sug['<mFo>'].item()),
                                                "Bfac": str(sug['<Bfac>'].item())
                                            }
                                            residue_data_list.append(residue_data_json)  
                                        else:
                                            residue_data_json = {
                                                "PDB": pdb_code,                            
                                                "Phi": entry[2],
                                                "Psi": entry[3],
                                                "Residue Data": {
                                                    "Residue Chain": residue_chain,
                                                    "Residue Position": AA_position,
                                                    "Amino Acid": AA
                                                },
                                                "Sugar Data": {
                                                    "Sugar Chain": sugar_chain,
                                                    "Sugar Position": sugar_position, 
                                                    "Sugar": sugar                               
                                                }
                                            }
                                            residue_data_list.append(residue_data_json)         
                                    elif output_format =='csv':
                                        if reduced_output == False:
                                            residue_data_list.append([pdb_code,residue_chain, AA_position, AA, sugar_chain, sugar_position, sugar, entry[2], entry[3],float(sug['RSCC'].values),float(sug['Q'].values),float(sug['Phi'].values),float(sug['Theta'].values),float(sug['Rsln'].values),str(sug['Detected type'].item()),str(sug['Cnf'].item()),float(sug['<mFo>'].values),float(sug['<Bfac>'].values),str(sug['Ctx'].item()),str(sug['Ok?'].item())])
                                        else:
                                            residue_data_list.append([pdb_code,residue_chain, AA_position, AA, sugar_chain, sugar_position, sugar, entry[2], entry[3]])
                                        
                                    else:
                                        print("The output format has not been detected correctly, check that the output_file = csv or output_file = json")  
                    
                    pdb_completed.append(pdb_code)
                        
            elif lineNo == 0:
                if output_format == 'csv':
                    if reduced_output == False:
                        residue_data_list.append(['PDB','Residue Chain','Residue Position','Amino Acid','Sugar Chain','Sugar Position','Sugar','Phi','Psi','Sugar RSCC','Q','Phi','Theta','Rsln','Detected Type','CNF','mFo','Bfac','Ctx','Ok?'])
                    else:
                        residue_data_list.append(['PDB','Residue Chain','Residue Position','Amino Acid','Sugar Chain','Sugar Position','Sugar','Phi','Psi'])
    with open(res_sug_output_file, "wb") as output:
        if output_format == 'json':
            json.dump(residue_data_list, output, indent=4)
        elif output_format == 'csv':
            writer = csv.writer(output, delimiter=',')
            writer.writerows(residue_data_list)

#CHOOSE WHICH OUTPUT YOU WANT
residue_sugar_read()
#sugar_sugar_read()



