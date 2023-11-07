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



data = 'GSE75748time'   #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)



temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE75748_sc_time_course_ec.csv.gz'):   #download raw data
    wget.download('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75748/suppl/GSE75748_sc_time_course_ec.csv.gz', out=temp_folder)   #download data
else:
    print('RAW data already in temporary folder')

    
with gzip.open(temp_folder + 'GSE75748_sc_time_course_ec.csv.gz') as f:
    raw_data  = pd.read_csv(f)


meta = pd.read_csv(data_process_path + 'code/meta_data/GSE75748_meta.csv')   #load meta data, processed SRA file
GEO_list = list(meta['GEO_Accession (exp)'])   #get the GEO accession name
aux = pd.read_csv(data_process_path + '/code/auxilary_file/' + 'GSE75748_sample.csv')
Exp_Name_list = list(raw_data.keys())[1:]
SRR_list = []


Sample_Name_list = []
Cell_type_list = [] 
for exp_name in Exp_Name_list:
    found = False
    for idx_s, sample in enumerate(aux['Accession']):
        title = aux['Title'][idx_s]
        title = title.split()[-1]
        title = title[1:-1]
        if exp_name == title:
            found = True
            break
    if found == False:
        print('Exp name in raw data and aux file does not match')
        break
    Sample_Name_list.append(sample)
    cell_type = exp_name.split('.')[1]
    cell_type = cell_type[:3]
    Cell_type_list.append(cell_type)

for sample in Sample_Name_list:
    for idx_geo, geo in enumerate(GEO_list):
        if sample == geo:
            found = True
            break
    if found == False:
        print('GEO Accession and Sample does not match')
        break
    SRR_list.append(meta['Run'][idx_geo])
    
        

#check to see if there is any sample name duplicates
if len(list(set(Sample_Name_list)))!= len(Exp_Name_list):
    print('Number of exp and sample does not match')

#construct the gene list
gene_list = list(raw_data['Unnamed: 0'])
for idx_gene, gene in enumerate(gene_list):
    gene_list[idx_gene] = gene


unique = list(set(Cell_type_list))
unique.sort()

N = len(Sample_Name_list)
M = len(gene_list)
print('Number of cells:', N)
print('Number of Genes:', M)

if len(raw_data) != M:
    print('Number of gene does not match raw data')



MATRIX = raw_data.values[:, 1:].astype(float)

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
