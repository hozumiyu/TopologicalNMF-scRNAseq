# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 20:46:18 2022

@author: yutah
"""

import numpy as np
import pandas as pd
import csv
import os
from sklearn.cluster import KMeans

def makeFolder(outpath):
    try:
        os.makedirs(outpath)
    except:
        return
    return

def load_X(data, inpath, data_process_path):
    if not os.path.exists(inpath + data + '/%s_data.csv'%(data)):
        os.system('python %s/main.py --data %s --process_directory %s --outpath %s'%(data_process_path, data, data_process_path, inpath) )
    inpath = inpath + data + '/'
    X = pd.read_csv(inpath + '%s_data.csv'%(data))
    X = X.values[:, 1:].astype(float)
    return X




def load_y(data, inpath, data_process_path):
    if not os.path.exists(inpath + data + '/%s_labels.csv'%(data)):
        os.system('python %s/main.py --data %s --process_directory %s --outpath %s'%(data_process_path, data, data_process_path, inpath) )
    inpath = inpath + data + '/'
    y = pd.read_csv(inpath + '%s_labels.csv'%(data))
    y = np.array(list(y['Label'])).astype(int)
    return y



def drop_sample(X, y):
    original = X.shape[1]
    labels = np.unique(y)
    good_index = []
    for l in labels:
        index = np.where(y == l)[0]
        if index.shape[0] > 15:
            good_index.append(index)
        else:
            print('label %d removed'%(l))
    good_index = np.concatenate(good_index)
    good_index.sort()
    new = good_index.shape[0]
    print(original - new, 'samples removed')
    return X[:, good_index], y[good_index]


def variance_cutoff(X, cutoff = 0.8):
    X = np.log2(1+X)
    index = np.arange(X.shape[0])
    
    variance = np.var(X, axis = 1)
    v_idx = np.where(variance > 1e-6)[0]
    new_index  = index[v_idx]
    
    variance = variance[new_index]
    
    M = variance.shape[0]
    cut = M - int(np.ceil(M * cutoff))
    print('Removing %d genes'%cut)
    new_larget_variance_index = np.argsort(variance)   #order the nonzero variance from smallest to largest
    print(variance[new_larget_variance_index])
    new_index = new_index[new_larget_variance_index[cut:]]
    new_index.sort()
    
    return X[new_index], new_index