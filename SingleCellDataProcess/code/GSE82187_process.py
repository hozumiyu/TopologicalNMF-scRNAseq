# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 20:34:29 2023

@author: yutah
"""

import wget, shutil, tarfile, gzip, sys
from auxilary import makeFolder, writeCellGeneCSV, writeCellLabelsCSV
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join, exists



data = 'GSE82187'   #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)



temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE82187_cast_all_forGEO.csv.gz'):   #download raw data
    wget.download('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82187/suppl/GSE82187_cast_all_forGEO.csv.gz', out=temp_folder)   #download data
else:
    print('RAW data already in temporary folder')

    
with gzip.open(temp_folder + 'GSE82187_cast_all_forGEO.csv.gz') as f:
    raw_data  = pd.read_csv(f)


meta = pd.read_csv(data_process_path + 'code/meta_data/%s_meta.csv'%( data))   #load meta data, processed SRA file
GEO_list = list(meta['GEO_Accession (exp)'])   #get the GEO accession name
aux = pd.read_csv(data_process_path + 'code/auxilary_file/' + 'GSE82187_sample.csv')
Exp_Name_list = list(raw_data.keys())[1:]
SRR_list = []


Exp_Name_list = []
Cell_type_list = []
for idx_p, protocol in enumerate(raw_data['protocol']):
    if protocol == 'Mic-scRNA-Seq':
        Exp_Name_list.append( raw_data['cell.name'][idx_p])
        Cell_type_list.append(raw_data['type'][idx_p])

Sample_Name_list = []
for exp_name in Exp_Name_list:
    found = False
    for idx_s, sample in enumerate(aux['Accession']):
        if exp_name in aux['Title'][idx_s]:
            found = True
            break
    if found == False:
        print('Exp name in raw data and aux file does not match')
    Sample_Name_list.append(sample)

no_dup = 0
for idx_s, sample in enumerate(Sample_Name_list):
    Exp_Name_list[idx_s] = Exp_Name_list[idx_s].split()[0]
    found = False
    for idx, geo in enumerate(meta['GEO_Accession (exp)']):
        if sample == geo:
            if found == False:
                found = True
                index = idx
            else:
                no_dup += 1
                if meta['major_cell_type'][idx] != meta['major_cell_type'][index]:
                    print('cell type of duplicate does not match!')
                    break
    if found == False:
        print('Sample name not found in meta data')
    else:
        SRR_list.append(meta['Run'][index])



unique = list(set(Cell_type_list))
unique.sort()

N = len(Sample_Name_list)
print('Number of cells:', N)


#Do a quick check to make sure that raw data and labels match

cell_name_check = list(raw_data['cell.name'])
cell_type_check = list(raw_data['type'])
#make sure that the cell names match
if cell_name_check == Exp_Name_list:
    print('Cell name of raw data and meta matches!!')
if cell_type_check == Cell_type_list:
    print('Cell type of raw data and meta matches!!')


#construct the gene list
gene_list = list(raw_data.keys())[5:]   #the first row contains the gene info
for idx_gene, gene in enumerate(gene_list):
    gene_list[idx_gene] = gene

MATRIX = raw_data.values[:705, 5:].astype(float)
MATRIX = MATRIX.T

if MATRIX.shape[0] == len(gene_list) and MATRIX.shape[1] == N:
    print('Dimension of the gene and cell matches!')
else:
    print('WARNING: DIMENSION OF THE MATRIX DOES NOT MATCH')



#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = unique

#write the dictionary
cell_type_dict = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}  #reverse dict from cell type to index

cell_count =  {cell_type_unique[i]:0 for i in range(len(cell_type_unique))}
Sample_Name_list = Sample_Name_list
Cell_label_list = []   #map cell type to number
for idx in range(len(Cell_type_list)): 
    cell_type = Cell_type_list[idx]     #get the current cell type
    try:   #see if the cell type exist in the dictionary
        cell_label = cell_type_dict[cell_type]
        cell_count[cell_type] += 1
        Cell_label_list.append(cell_label)
        found = True
    except:
        found = False    #if it is not found, print message
        print(Sample_Name_list[idx], cell_type, 'Wrong cell type')
print('Cell Count')
for k in cell_count.keys():
    print('%s:'%k, cell_count[k])

outfile = outpath + '%s_data.csv'%(data)
writeCellGeneCSV(outfile, Sample_Name_list, gene_list, MATRIX)
 
outfile = outpath + '%s_labels.csv'%(data)
writeCellLabelsCSV(outfile, Sample_Name_list, SRR_list, Cell_type_list, Cell_label_list)
shutil.rmtree(temp_temp_folder) 
