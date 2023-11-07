# -*- coding: utf-8 -*-
"""
Created on Sun May 14 13:41:21 2023

@author: Yuta
"""

import sys, os


#option = int(sys.argv[1])
#random_state = int(sys.argv[2])

data_process_path = './SingleCellDataProcess/'   #change this line if you have a different path for the data processing
data_path = './data/'  #location of where the processed raw data will be saved


data_vec = ['GSE45719', 'GSE67835', 'GSE75748cell', 'GSE75748time', 'GSE75748cell', 'GSE82187', 'GSE84133human4', 'GSE84133mouse1', 'GSE94820', 'GSE57249']


for data in data_vec:
    os.system('python main.py --data ' + data)

