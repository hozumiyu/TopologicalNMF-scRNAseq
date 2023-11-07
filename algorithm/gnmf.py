# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:04:41 2023

@author: yutah
"""

import numpy as np



class GNMF():
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
    
    def compute_loss(self, W, H):
        loss = np.linalg.norm( self.X -W@ H, ord = 'fro')**2 + self.l*np.trace( H@self.L@H.T)
        return loss
    
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
        
        
    
    def updateW(self, W, H):
        W = W * ( ( self.X @H.T) / ( W@H@H.T + 1e-9)  )
        return W
    
    
    def updateH(self, W, H):
        numerator = W.T @self.X + 2*self.l*H@self.A
        denominator = W.T@W@ H + 2*self.l *H@self.D +  1e-9
        H = H * numerator / denominator
        return H
    
        
    def multiplicative_update(self, W, H):
        W = self.updateW(W,H)
        H = self.updateH(W,H)
        scale = np.linalg.norm(W, axis = 0)
        W = W / scale[None, :]
        H = H * scale[:, None]
        return W, H
    
    
    def fit_transform(self, X, W, H, L, A ,D):
        self.X = X
        self.L = L
        self.A = A
        self.D  =D 
        #self.adjacency_diagonal()
        iteration = 1
        loss_prev = 1e9; curr_loss = self.compute_loss( W, H)
        print(self.L)
        while iteration < self.max_iter and np.abs(loss_prev - curr_loss) / loss_prev > self.tol:
            loss_prev = curr_loss
            H_prev = H.copy()
            W, H = self.multiplicative_update(W, H)
            curr_loss = self.compute_loss( W, H)
            print('Iteration:', iteration, 'Current loss:', curr_loss, np.abs(loss_prev - curr_loss) / loss_prev)
            iteration+=1
        
        if iteration == self.max_iter:
            print('maximum iteration reached')
        return W, H