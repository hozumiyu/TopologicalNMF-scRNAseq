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



data = 'GSE45719', #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)

temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE45719_RAW.tar'):   #download raw data
    wget.download('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45719&format=file', out=temp_folder)   #download data
else:
    print('RAW data already in temporary folder')

with tarfile.open(temp_folder + 'GSE45719_RAW.tar') as tf:
    tf.extractall(path=temp_temp_folder)
    
    
meta = pd.read_csv(data_process_path + 'code/meta_data/%s_meta.csv'%( data))   #load meta data, processed SRA file
source_name_list = meta['source_name']   #source name contains the cell type for GSE45719
GEO_list = list(meta['GEO_Accession (exp)'])   #get the GEO accession name


Cell_type_list = []   #cell type list
Sample_Name_list = []    #sample name
SRR_list = []
count = 0
for idx, source in enumerate(source_name_list):   #Only extract the GEO corresponding to the cell types
    if '2-cell' in source:
        Cell_type_list.append('2-cell')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif '4-cell' in source:
        Cell_type_list.append('4-cell')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif '8-cell' in source:
        Cell_type_list.append('8-cell')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif '16-cell' in source:
        Cell_type_list.append('16-cell')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif 'Early blastocyst' in source:
        Cell_type_list.append('Early blastocyst')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif 'Mid blastocyst' in source:
        Cell_type_list.append('Mid blastocyst')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif 'Late blastocyst' in source:
        Cell_type_list.append('Late blastocyst')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    elif 'fibroblast' in source:
        Cell_type_list.append('fibroblast')
        Sample_Name_list.append(GEO_list[idx])
        SRR_list.append(meta['Run'][idx])
    else:
        print(source)
        count += 1

unique_cell_type = list(set(list(Cell_type_list)))   #get the list of the unique cell types



N = len(Sample_Name_list)
print('Number of cells:', N)

onlyfiles = [f for f in listdir(temp_temp_folder) if isfile(join(temp_temp_folder, f))]   #the path of all the raw files

#construct the gene list
gene_list = []
sample = Sample_Name_list[0]
found = False
for file in onlyfiles:
    if sample in file:
        found = True
        break
if not found:
    print(sample, 'not found')

raw_data_file = gzip.open(temp_temp_folder + file, 'r')
raw_data_lines = raw_data_file.readlines()
raw_data_file.close()
gene_list = []
dup_gene = []
for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first line is header
        raw_data = raw_data.split()[0]  #first column is the gene name
        if raw_data.decode('utf-8') not in gene_list:
            gene_list.append(raw_data.decode('utf-8'))
        else:
            dup_gene.append(raw_data.decode('utf-8') )

print('Revmoved', len(dup_gene), 'due to duplicate')

#Construct the data matrix
M = len(gene_list)
MATRIX = np.zeros([M, N])
for idx_cell in range(N):
    sample = Sample_Name_list[idx_cell]
    for file in onlyfiles:
        if sample in file:
            found = True
            break
    if not found:
        print(sample, 'not found')
    raw_data_file = gzip.open(temp_temp_folder + file, 'r')
    raw_data_lines = raw_data_file.readlines()
    raw_data_file.close()
    idx_gene = 0
    for idx_raw, raw_data in enumerate(raw_data_lines):
        if idx_raw > 0: #first line is 
            raw_data = raw_data.split()
            gene = raw_data[0].decode('utf-8')
            if gene == gene_list[idx_gene]:
                MATRIX[idx_gene, idx_cell] = float(raw_data[2].decode('utf-8'))
                idx_gene += 1
    if idx_gene != M:
        print('ERROR: Some of the gene index is wrong')


#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = ['2-cell', 'Late blastocyst', 'fibroblast', '8-cell', '16-cell', '4-cell', 'Early blastocyst', 'Mid blastocyst']

#write the dictionary
cell_type_dict = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}  #reverse dict from cell type to index

cell_count =  {cell_type_unique[i]:0 for i in range(len(cell_type_unique))}
Sample_Name_list = Sample_Name_list
Cell_label_list = []   #map cell type to number
for idx in range(len(Cell_type_list)):  #iterate over all the snd(raw_data.decode('utf-8').upper())
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