# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 20:34:29 2023

@author: yutah
"""

import wget, shutil, tarfile, gzip, sys
from auxilary import makeFolder, writeCellGeneCSV, writeCellLabelsCSV
import pandas as pd
import numpy as np
from os import listdir, makedirs
from os.path import isfile, join, exists



data = 'GSE84133'   #data name

data_process_path = sys.argv[2] 

temp_folder = data_process_path + '/temporary/'; makeFolder(temp_folder)



temp_temp_folder = temp_folder + '%s_temporary/'%(data); makeFolder(temp_temp_folder)

if not exists(temp_folder + 'GSE84133_RAW.tar'):   #download raw data
    wget.download('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file', out=temp_folder)   #download data
else:
    print('RAW data already in temporary folder')

with tarfile.open(temp_folder + 'GSE84133_RAW.tar') as tf:
    tf.extractall(path=temp_temp_folder)
    
    
meta = pd.read_csv(data_process_path + 'code/meta_data/GSE84133_meta.csv')   #load meta data, processed SRA file
source_name_list = meta['source_name']   #source name contains the cell type for GSE45719
GEO_list = list(meta['GEO_Accession (exp)'])   #get the GEO accession name

exp_list = ['human1', 'human2', 'human3', 'human4', 'mouse1', 'mouse2']
GSM_list = [   'GSM2230757', 'GSM2230758', 'GSM2230759', 'GSM2230760', 'GSM2230761', 'GSM2230762']
onlyfiles = [f for f in listdir(temp_temp_folder) if isfile(join(temp_temp_folder, f))]
#meta data
for idx_e, exp in enumerate(exp_list):
    GSM = GSM_list[idx_e]
    for file in onlyfiles:
        if file[-2:] == 'gz':
            if exp in file:
                outpath2 = sys.argv[1]  + '%s%s/'%(data, exp); makeFolder(outpath2)
                raw_data = pd.read_csv(temp_temp_folder + file)
                #headers are the genes
                gene_list = list(raw_data.keys())[3:]
                for idx_gene, gene in enumerate(gene_list):
                    gene_list[idx_gene] = gene
                if len(gene_list) != len(list(raw_data.keys())) - 3:
                    print('Gene list size does not match raw data')
                Sample_Name_list = list(raw_data['Unnamed: 0'])
                Cell_type_list = list(raw_data['assigned_cluster'])
                unique = list(set(Cell_type_list))
                unique.sort()
                

                N = len(Sample_Name_list)
                M = len(gene_list)
                
                SRR_list = N*[GSM]
                if len(raw_data) != N:
                    print('Number of cells do not match')
                
                print('Current exp:', exp)
                print('Number of cells:', N)
                print('Number of Genes:', M)
                print('Number of cell type:', len(unique))
                
                MATRIX = raw_data.values[:, 3:].astype(float)
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
                Cell_label_list = []
                for idx in range(len(Cell_type_list)):  #iterate over all the samples
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

                outfile = outpath2 + '%s%s_data.csv'%(data, exp)
                writeCellGeneCSV(outfile, Sample_Name_list, gene_list, MATRIX)
                 
                outfile = outpath2 + '%s%s_labels.csv'%(data, exp)
                writeCellLabelsCSV(outfile, Sample_Name_list, SRR_list, Cell_type_list, Cell_label_list)
shutil.rmtree(temp_temp_folder) 