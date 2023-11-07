# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:04:41 2023

@author: yutah
"""

import numpy as np



class rGNMF():
    def __init__(self, n_components, l, max_iter = 1000, tol = 1e-3, random_state = 1):
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
        self.random_state = 1
        self.max_iter = max_iter
        self.tol = tol
        self.l = l
        return        
    
    def adjacency_diagonal(self):
        D = np.diagonal(self.L)#.copy()
        A = self.L.copy()
        np.fill_diagonal(A, 0)
        mask = np.zeros(A.shape)
        mask[A != 0] = 1
        A = -A
        A[mask == 0]  = 0
        self.D = D
        self.A = A
        return 
    
    def compute_loss(self, W, H):
        loss = np.sum(np.linalg.norm( self.X - W@H, axis = 0)) + self.l*np.trace( H@self.L@H.T)
        return loss
    
    def updateW(self, W, H, Q):
        QH = Q  @ H.T
        W = W * (( self.X@QH) / ( W@ H @ QH   + 1e-9) )
        return W
    
    
    def updateH(self, W, H, Q):
        H = H * (  (W.T @ self.X@ Q + 2*self.l*H@self.A )   /   (  W.T@W@ H@ Q + 2*self.l *H@self.D+ 1e-9)  )
        return H
    
    def updateQ(self, W, H):
        Q = np.diag(1/np.linalg.norm(self.X - np.dot(W, H), axis = 0))
        return Q
    
    
    def multiplicative_update(self, W, H):
        Q = self.updateQ(W, H)
        W = self.updateW(W,H, Q)
        H = self.updateH(W,H, Q)
        scale = np.linalg.norm(W, axis = 0)
        W = W / scale[None, :]
        H = H * scale[:, None]
        return W, H
    
    
    def fit_transform(self, X, W, H, L, A, D):
        self.L = L
        self.A = A
        self.D = D
        #self.adjacency_diagonal()
        self.X = X
        iteration = 1
        loss_prev = 1e9; curr_loss = self.compute_loss( W, H)
        
        while iteration < self.max_iter and np.abs(loss_prev - curr_loss) / loss_prev > self.tol:
            loss_prev = curr_loss
            W, H = self.multiplicative_update(W, H)
            curr_loss = self.compute_loss( W, H)
            print('Iteration:', iteration, 'Current loss:', curr_loss)
            iteration+=1
        
        if iteration == self.max_iter:
            print('maximum iteration reached')
        return W, H