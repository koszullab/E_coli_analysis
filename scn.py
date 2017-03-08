# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 16:58:59 2014

@author: axel
"""
import numpy as np
import matplotlib.pyplot as plt
from pylab import *  
    
# Normalisation function
# A is the matrix to normalise  
# the treshold to remove bins from the normalization but not from the matrice  
    
def scn_func(A,threshold):    
    n1 = A.shape[0];
    n_iterations=1;
    
    B = np.zeros((n1, n1));
    keep = np.zeros((n1, 1));

    b=A.sum(axis=0)
    hist, bins = np.histogram(b,bins=100)
    width = 0.7 * (bins[1] - bins[0]);
    center = (bins[:-1] + bins[1:]) / 2;
    #plt.bar(center, hist, align='center', width=width)
    #plt.show();
        
    for i in range(0,n1) :
        if np.sum(A[i,]) > threshold:
            keep[i] = 1
        else :
            keep[i] = 0
      
    for n in range(0,n_iterations) :
        print(n)
        for i in range(0,n1) :
            for j in range(0,n1) :
                if keep[i] == 1:
                    B[i,j]=A[i,j]/ np.sum(A[i,]) 
                else :
                    B[i,j]=0
                    
        print(n)        
        for i in range(0,n1) :
           for j in range(0,n1) :
               if keep[j] == 1 :
                   A[i,j] = B[i,j] / np.sum(B[:,j])
               else :
                   A[i,j]=0            
            
    return A
    