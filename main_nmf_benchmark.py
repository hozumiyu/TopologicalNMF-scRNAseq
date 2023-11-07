# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:37:34 2023

@author: yutah
"""

import numpy as np
from algorithm.auxilary import load_X, load_y, drop_sample, makeFolder
from algorithm.clustering import computeClusteringScore
from algorithm.nmf import NMF
from algorithm.rnmf import rNMF
from algorithm.graphs import computeGraphLaplacian
from algorithm.gnmf import GNMF
from algorithm.rgnmf import rGNMF
from algorithm.initialization import initialize
import sys, os

def write_results(outfile, param, ari, nmi, purity, acc):
    file = open(outfile, 'w')
    file.writelines('Param,ARI,NMI' + '\n')
    file.writelines('%s,%.4f,%.4f,%.4f,%.4f'%(param, np.mean(ari), np.mean(nmi), np.mean(purity), np.mean(acc)))
    file.close()
    return

data_path = './SingleCellDataProcess/data/' 
data_process_path = './SingleCellDataProcess/'
init_path = './initialization/'

data = sys.argv[1]



X = load_X(data, data_path, data_process_path)
y = load_y(data, data_path, data_process_path)
if data != 'GSE57249':
    X, y = drop_sample(X, y)
n_clusters = np.unique(y).shape[0]

M, N = X.shape
    
X = np.log10(1+X)
X = X/ np.linalg.norm(X, axis = 0) [None , :]    

W, H = initialize(init_path, data, X)



myNMF = NMF(n_components = W.shape[1])
W, H = myNMF.fit_transform(X, W, H)
Ht = H.T
ari_nmf, nmi_nmf, purity_nmf, acc_nmf, LABELS = computeClusteringScore(Ht, y, max_state = 10)
    
    
W, H = initialize(init_path, data, X)
myrNMF = rNMF(n_components = W.shape[1])
W, H = myrNMF.fit_transform(X, W, H)
Ht = H.T
ari_rnmf, nmi_rnmf, purity_rnmf, acc_rnmf, LABELS = computeClusteringScore(Ht, y, max_state = 10)


L, A, D = computeGraphLaplacian(X.T, n_neighbors = 8)

W, H = initialize(init_path, data, X)
myGNMF = GNMF(n_components = W.shape[1], l = 1)
W, H = myGNMF.fit_transform(X, W, H, L, A, D)
ari_gnmf, nmi_gnmf, purity_gnmf, acc_gnmf, LABELS = computeClusteringScore(Ht, y, max_state = 10)


W, H = initialize(init_path, data, X)
myrGNMF = rGNMF(n_components = W.shape[1], l = 1)
W, H = myrGNMF.fit_transform(X, W, H, L, A, D)
ari_rgnmf, nmi_rgnmf, purity_rgnmf, acc_rgnmf, LABELS = computeClusteringScore(Ht, y, max_state = 10)

print('Scores for %s \n'%data )

print('Method \t ARI \t NMI \t Purity \t ACC \n')
print('NMF', ari_nmf.round(4), nmi_nmf.round(4), purity_nmf.round(4), acc_nmf.round(4))
print('rNMF', ari_rnmf.round(4), nmi_rnmf.round(4), purity_rnmf.round(4), acc_rnmf.round(4))
print('GNMF', ari_gnmf.round(4), nmi_gnmf.round(4), purity_gnmf.round(4), acc_gnmf.round(4))
print('rGNMF', ari_rgnmf.round(4), nmi_rgnmf.round(4), purity_rgnmf.round(4), acc_rgnmf.round(4))