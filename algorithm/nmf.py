# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:04:41 2023

@author: yutah
"""

import numpy as np



class NMF():
    def __init__(self, n_components, max_iter = 1000, tol = 1e-3):
        '''

        Parameters
        ----------
        n_components : int
            Number of components.
        max_iter : int, optional
            maximum number of iteration before the algorithm stops. The default is 1000.
        tol : float, optional
            relative error tolerance of the loss function. The default is 1e-5.
        random_state : int, optional
            used for random initialization. The default is 1.

        Returns
        -------
        None.

        '''
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        return        
    
    def compute_loss(self, W, H):
        loss = np.linalg.norm( self.X - np.dot(W, H), ord = 'fro')
        return loss
    
    def updateW(self, W, H):
        W = W * ( (np.dot(self.X, H.T)) / ( np.dot(np.dot(W, H), H.T) + 1e-6)  )
        return W
    
    
    def updateH(self, W, H):
        H = H * (  (np.dot(W.T,self.X))  /   (  np.dot(W.T,np.dot(W, H)) + 1e-6)  )
        return H
    
        
    def multiplicative_update(self, W, H):
        W = self.updateW(W,H)
        H = self.updateH(W,H)
        scale = np.linalg.norm(W, axis = 0)
        W = W / scale[None, :]
        H = H * scale[:, None]
        return W, H
    
    
    def fit_transform(self, X, W, H):
        self.X = X
        iteration = 1
        loss_prev = 1e9; curr_loss = self.compute_loss( W, H)
        
        while iteration < self.max_iter and np.abs(loss_prev - curr_loss) / loss_prev > self.tol:
            loss_prev = curr_loss
            W, H = self.multiplicative_update(W, H)
            curr_loss = self.compute_loss( W, H)
            print('Iteration:', iteration, 'Current loss:', curr_loss, 'dloss:',np.abs(loss_prev - curr_loss) / loss_prev)
            iteration+=1
        
        if iteration == self.max_iter:
            print('maximum iteration reached')
        return W, H