# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:23:38 2023

@author: Yuta
"""

import numpy as np
from sklearn.utils.extmath import squared_norm #,randomized_svd
from math import sqrt
import os

'''
This is nndsvd implementation taken from sklearn
'''

def makeFolder(outpath):
    try:
        os.makedirs(outpath)
    except:
        return
    return
    
def norm(x):
    """Dot product-based Euclidean norm implementation.
    See: http://fa.bianp.net/blog/2011/computing-the-vector-norm/
    Parameters
    ----------
    x : array-like
        Vector for which to compute the norm.
    """
    return sqrt(squared_norm(x))

def nndsvd(X, n_components, eps = 1e-8):
    '''
        This is modification from the nndsvd initialization from sklearn
        Instead of using randomised, actual svd is used
    '''
    #U, S, V = randomized_svd(X, n_components = n_components, random_state=random_state)
    u,s,vh = np.linalg.svd(X, full_matrices= False)
    U = u[:,:n_components]
    V = vh[:n_components,:]
    S = s[:n_components]
    W = np.zeros_like(U)
    H = np.zeros_like(V)

    # The leading singular triplet is non-negative
    # so it can be used as is for initialization.
    W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
    H[0, :] = np.sqrt(S[0]) * np.abs(V[0, :])

    for j in range(1, n_components):
        x, y = U[:, j], V[j, :]

        # extract positive and negative parts of column vectors
        x_p, y_p = np.maximum(x, 0), np.maximum(y, 0)
        x_n, y_n = np.abs(np.minimum(x, 0)), np.abs(np.minimum(y, 0))

        # and their norms
        x_p_nrm, y_p_nrm = norm(x_p), norm(y_p)
        x_n_nrm, y_n_nrm = norm(x_n), norm(y_n)

        m_p, m_n = x_p_nrm * y_p_nrm, x_n_nrm * y_n_nrm

        # choose update
        if m_p > m_n:
            u = x_p / x_p_nrm
            v = y_p / y_p_nrm
            sigma = m_p
        else:
            u = x_n / x_n_nrm
            v = y_n / y_n_nrm
            sigma = m_n

        lbd = np.sqrt(S[j] * sigma)
        W[:, j] = lbd * u
        H[j, :] = lbd * v
    avg = np.mean(X)
    W[W < eps] = avg
    H[H < eps] = avg
    return W, H

def initialize(init_path, data, X):
    path = init_path + data + '_init/'
    file = data 
    if not os.path.exists(path + file + '_Winit.npy') or not os.path.exists(path + file + '_Hinit.npy'):
        makeFolder(path)
        N = X.shape[1]
        v = int(np.ceil(np.sqrt(N)))
        W, H = nndsvd(X, n_components = v, eps = 1e-8)
        np.save(path + file + '_Winit.npy', W)
        np.save(path + file + '_Hinit.npy', H)
    else:
        W = np.load(path + file + '_Winit.npy')
        H = np.load(path + file + '_Hinit.npy')
    return W, H