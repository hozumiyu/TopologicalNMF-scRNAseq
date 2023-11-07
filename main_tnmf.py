# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:37:34 2023

@author: yutah
"""

import numpy as np
from algorithm.auxilary import load_X, load_y, drop_sample, makeFolder
from algorithm.clustering import computeClusteringScore
from algorithm.gnmf import GNMF
from algorithm.rgnmf import rGNMF
from algorithm.initialization import initialize
from algorith.pl_graph import computePersistentLaplacian
from algorith.pl_graph_cutoff import computePersistentLaplacian_cutoff
import sys, os


data_path = './SingleCellDataProcess/data/' 
data_process_path = './SingleCellDataProcess/'
init_path = './initialization/'

data = sys.argv[1]
method = sys.argv[2]
n_neighbors = sys.argv[3]
weights = sys.argv[4]
l = float(sys.argv[5])
weights = weights.split(',')
weights = np.array(weights).astype(float)


X = load_X(data, data_path, data_process_path)
y = load_y(data, data_path, data_process_path)
if data != 'GSE57249':
    X, y = drop_sample(X, y)
else:
    X = np.log10(1+X).T
n_clusters = np.unique(y).shape[0]

M, N = X.shape
    
X = np.log10(1+X)
X = X/ np.linalg.norm(X, axis = 0) [None , :]    

W, H = initialize(init_path, data, X)

if method == 'kTNMF':
    L, A, D = computePersistentLaplacian(X, n_neighbors, weights)
    myNMF = GNMF(n_components = W.shape[1], l = l)
    W, H = myNMF.fit_transform(X, W, H, L, A, D)
elif method == 'TNMF':
    L, A, D = computePersistentLaplacian_cutoff(X, n_neighbors, weights)
    myNMF = GNMF(n_components = W.shape[1], l = l)
    W, H = myNMF.fit_transform(X, W, H, L, A, D)
elif method == 'krTNMF':
    L, A, D = computePersistentLaplacian(X, n_neighbors, weights)
    myNMF = rGNMF(n_components = W.shape[1], l = l)
    W, H = myNMF.fit_transform(X, W, H, L, A, D)
elif method == 'rTNMF':
    L, A, D = computePersistentLaplacian_cutoff(X, n_neighbors, weights)
    myNMF = rGNMF(n_components = W.shape[1], l = l)
    W, H = myNMF.fit_transform(X, W, H, L, A, D)

Ht = H.T
ari_nmf, nmi_nmf, purity_nmf, acc_nmf, LABELS = computeClusteringScore(Ht, y, max_state = 10)
    
  

print('Scores for %s \n'%data )

print('Method \t ARI \t NMI \t Purity \t ACC \n')
print(method , ari_nmf.round(4), nmi_nmf.round(4), purity_nmf.round(4), acc_nmf.round(4))