# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 20:42:44 2023

@author: yutah
"""

import os, argparse




parser = argparse.ArgumentParser(description='scRNA-seq data processing')
parser.add_argument('--data', type=str, default = 'GSE64016',
                    help='data for preprocessing')
parser.add_argument('--process_directory', type=str,  default= './',
                    help='input the SingleCellDataProcess directory')
parser.add_argument('--outpath', type=str,  default= './data/',
                    help='input the output directory')

args = parser.parse_args()

data = args.data
data_process_path = args.process_directory
outpath = args.outpath


if data[:8] in 'GSE84133':
    data = 'GSE84133'

os.system('python %s/code/%s_process.py %s %s'%(data_process_path, data, outpath, data_process_path))