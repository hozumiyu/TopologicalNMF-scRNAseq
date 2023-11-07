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
import os 
import shutil
#this data does not have a standard meta data. Series list is used. SRA run slector is not available

data = 'GSE94820'   #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)



temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz'):   #download raw data
    wget.download('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94820&format=file&file=GSE94820%5Fraw%2EexpMatrix%5FDCnMono%2Ediscovery%2Eset%2Esubmission%2Etxt%2Egz', out=temp_folder)   #download data
    #os.rename(temp_folder + 'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz', temp_folder + 'GSE94820_raw.txt.gz')
else:
    print('RAW data already in temporary folder')


with gzip.open(temp_folder + 'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz', 'rb') as f:
    with open(temp_temp_folder + 'GSE94820.txt', 'wb') as f_out:
        shutil.copyfileobj(f, f_out)
    raw_data_file = open(temp_temp_folder + 'GSE94820.txt')
    raw_data_lines = raw_data_file.readlines()
    raw_data_file.close()



aux = pd.read_csv(data_process_path + 'code/auxilary_file/' + 'GSE94820_sample.csv')
#ge the exp name from the raw data, and the cell type and sample name from aux
Exp_Name_list = raw_data_lines[0].split()
Cell_type_list = []
Sample_Name_list = []
#map the exp name to class names
for idx, exp in enumerate(Exp_Name_list):
    for idx_title, title in enumerate(aux['Title']):
        if exp == title:
            found = True
            break
    cell_type =  title.split('_')
    cell_type = cell_type[0]
    
    Cell_type_list.append(cell_type)
    Sample_Name_list.append(aux['Accession'][idx_title])
    if found == False:
        print('WARNING', exp, 'not found in aux!')
        break

if len(Sample_Name_list) != len(list(set(Sample_Name_list))):
    print('WARNING: There is dupolicates in the sample name. consider reprocessing')
    
#construct the gene list
gene_list = []   #the first col contains the gene info
for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first row is header
        gene = raw_data.split()[0]
        gene_list.append(gene.upper())

unique = list(set(Cell_type_list))
unique.sort()

M = len(gene_list)
N = len(Sample_Name_list)
print('Number of genes:', M)
print('Number of cells:', N)
print('Number of cell type:', len(unique))  #should be 4

if M != len(raw_data_lines) - 1:
    print('Gene length and raw data does not match')

MATRIX = np.zeros([M,N])

for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first line is the header
        idx_row = idx_raw -1 
        raw_data = raw_data.split('\t')[1:]
        
        rowdata = np.array(raw_data).astype(float)
        MATRIX[idx_row, : ] = rowdata




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
writeCellLabelsCSV(outfile, Sample_Name_list, Sample_Name_list, Cell_type_list, Cell_label_list)
shutil.rmtree(temp_temp_folder) 
