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



data = 'GSE57249'   #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)

temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)


if not exists(temp_folder + 'GSE57249_fpkm.txt.gz'):   #download raw data
    wget.download('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57249/suppl/GSE57249_fpkm.txt.gz', out=temp_folder)   #download data
else:
    print('RAW data already in temporary folder')

    
with gzip.open(temp_folder + 'GSE57249_fpkm.txt.gz') as f:
    lines = f.readlines()


meta = pd.read_csv(data_process_path + 'code/meta_data/%s_meta.csv'%( data))   #load meta data, processed SRA file
GEO_list = list(meta['GEO_Accession (exp)'])   #get the GEO accession name
Exp_Name_list = lines[0].split()[1:]   #the first row of the fpkm file is the sample name
SRR_list = []

unique = [ 'zygote',  'Two-cell Embryo Blastomere', 'Four-cell Embryo Blastomere']
Sample_Name_list = []
Cell_type_list = [] 

exp_index = []   #used to extract column from the raw file
for idx_e, exp_name in enumerate(Exp_Name_list):
    found = False
    exp_name = exp_name.decode("utf-8")  
    for idx_s, sample in enumerate(GEO_list):  #run through GEO_list to find matching GSM name
        if exp_name == sample:  
            found = True
            break
    if found == False:
        print('Cannot find matching sample names in the meta data')
    
    if meta['Developmental_stage'][idx_s] in unique:
        Sample_Name_list.append(sample)
        Cell_type_list.append(meta['Developmental_stage'][idx_s])
        SRR_list.append(meta['Run'][idx_s])
        exp_index.append(idx_e)
exp_index = np.array(exp_index)

    
    
        

#check to see if there is any sample name duplicates
if len(list(set(Sample_Name_list)))!= len(Exp_Name_list):
    print('Number of exp and sample does not match')

#construct the gene list. First column of the data contains gene info
gene_list = []
for idx, line in enumerate(lines):
    if idx != 0:
        line = line.split()
        gene_list.append(line[0].decode("utf-8")  )




N = len(Sample_Name_list)
M = len(gene_list)
print('Number of cells:', N)
print('Number of Genes:', M)


if N != 56:
    print('There should be 44 samples. You have %d samples'%(N))

MATRIX = np.zeros([M, N])

for idx, line in enumerate(lines):
    if idx > 0:
        line = line.split()
        r = np.array(line[1:]).astype(float)
        r = r[exp_index]
        MATRIX[idx-1, :] = r

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
