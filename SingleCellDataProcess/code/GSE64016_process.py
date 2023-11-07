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

data = 'GSE64016'   #data name
outpath = sys.argv[1]  + '%s/'%(data)
data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)
makeFolder(outpath)


temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE64016_H1andFUCCI_normalized_EC.csv.gz'):   #download raw data
    wget.download('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64016&format=file&file=GSE64016%5FH1andFUCCI%5Fnormalized%5FEC%2Ecsv%2Egz', out=temp_folder)   #download data
    #os.rename(temp_folder + 'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz', temp_folder + 'GSE94820_raw.txt.gz')
else:
    print('RAW data already in temporary folder')


with gzip.open(temp_folder + 'GSE64016_H1andFUCCI_normalized_EC.csv.gz') as f:
    raw_record = pd.read_csv(f)


meta = pd.read_csv(data_process_path + 'code/meta_data/' + 'GSE64016_meta.csv')
aux = pd.read_csv(data_process_path + 'code/auxilary_file/' + 'GSE64016_sample.csv')
#ge the exp name from the raw data, and the cell type and sample name from aux
Exp_Name_list = list(raw_record.keys())[1:]
Cell_type_list = []
Sample_Name_list = []
SRR_list = []
#map the exp name to class names
for idx, exp in enumerate(Exp_Name_list):
    found = False
    for idx_title, title in enumerate(aux['Title']):
        if exp == title:
            found = True
            break
    if found == False:
        print('WARNING', exp, 'not found in aux!')
        break
    cell_type =  exp.split('_')
    cell_type = cell_type[0]
    
    Cell_type_list.append(cell_type)
    sample = aux['Accession'][idx_title]
    Sample_Name_list.append(sample)
    
    found = False
    for idx_geo, geo in enumerate(meta['GEO_Accession (exp)']):
        if geo == sample:
            found = True
            break
    if found == False:
        print('WARNING', 'GEO', 'not found in meta!')
        break
    SRR_list.append(meta['Run'][idx_geo])

if len(Sample_Name_list) != len(list(set(Sample_Name_list))):
    print('WARNING: There is dupolicates in the sample name. consider reprocessing')

print(raw_record.keys())
#construct the gene list
gene_list = raw_record['Unnamed: 0']   #the first col contains the gene info

unique = list(set(Cell_type_list))
unique.sort()

M = len(gene_list)
N = len(Sample_Name_list)
print('Number of genes:', M)
print('Number of cells:', N)
print('Number of cell type:', len(unique))  #should be 4


MATRIX = raw_record.values[:, 1:].astype(float)

if MATRIX.shape[1] != N:
    print('sample numbers do not match')




#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = unique
print(cell_type_unique)
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
