# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 20:35:38 2023

@author: yutah
"""

import os, csv



def makeFolder(outpath):
    try:
        os.makedirs(outpath)
    except:
        print('outfolder already exists')
        
    return
    

def writeCellGeneCSV(outfile, Sample_Name_list, gene_list, MATRIX):
    with open(outfile, "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow( [None] + Sample_Name_list)
        for idx_row in range(len(gene_list)):
            row = list(MATRIX[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
            
            
def writeCellLabelsCSV(outfile, Sample_Name_list, SRR_list, Cell_type_list, Cell_label_list):
    with open(outfile, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sample Name', 'SRA run', 'Cell type', 'Label'])
        for idx in range(len(Cell_type_list)):
            writer.writerow([Sample_Name_list[idx], SRR_list[idx], Cell_type_list[idx] , Cell_label_list[idx]])