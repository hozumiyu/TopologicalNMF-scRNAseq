#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 14:38:48 2022

@author: yutaho
"""

import numpy as np
import csv
import pandas as pd
from os import listdir
from os.path import isfile, join
import gzip
import sys


def outpath(data):
    import os
    try:
        os.makedirs('../%s/'%(data))
    except:
        return
    return
def constructMetaLine(line, column_length):
    line = line.split(',')
    newline = []
    contain = False
    
    for idx in range(len(line)):
        if '"' in line[idx]:
            if contain == False:
                new_word = line[idx]
                contain = True
            else:
                contain = False
                new_word = new_word +  '_' + line[idx]
                newline.append(new_word)
        else:
            if contain == True:
                new_word = new_word +  '_' + line[idx]
            else:
                new_word = line[idx].strip()
                if new_word == '' and column_length > 0:
                    newline.append('-')
                else:
                    newline.append(new_word)
    return newline


def constructMetaData():
    file = open('SraRunTable.txt', 'r')
    sra = file.readlines()
    file.close()
    length = len(sra)
    with open('meta.csv', "w", newline = '') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for idx_length in range(length):
            if idx_length == 0:
                newline = constructMetaLine(sra[idx_length], column_length = 0)
                #line = sra[idx_length]#.decode('utf8')
                #line = line.split(',')
                #column_length = len(line)
                column_length = len(newline)
                writer.writerow( newline)
            
            else:
                newline = constructMetaLine(sra[idx_length], column_length)
                if len(newline) > column_length:
                    print('Too many entry at idx', idx_length)
                writer.writerow( newline)
    '''
    df = pd.read_csv('meta.csv'%(data,data))
    df.drop_duplicates()
    df.to_csv('../%s/%s_meta.csv'%(data,data), index = False)
    '''
if __name__ == "__main__":
    constructMetaData()