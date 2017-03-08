# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:11:54 2015
@author: Axel KournaK
"""

# MIRNY type nomalisation:
# mij = mij * (total reads) / (total reads in bin i * total reads in bin j)
# formula used in the science paper of Le et al. 

import numpy as np
import matplotlib.pyplot as plt  
from pylab import *  
    
# A is the matrix to normalise  
# july 2016: adding of a threshold filtering 
    
def ice_func(A,threshold): 
    n = A.shape[0];
    n_iterations=51;
    
    B = np.zeros((n, n))
    keep = np.zeros((n, 1))

    b=A.sum(axis=0)
    hist, bins = np.histogram(b,bins=100)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    #plt.bar(center, hist, align='center', width=width)
    #plt.show();
     
    total_reads= A.sum();
    print(total_reads)
    
    for i in range(0,n_iterations) :
        S1= np.multiply( np.ones((n,n)),  np.sum(A,axis=0) );
        S2=( np.multiply( np.ones((n,n)),  np.sum(A,axis=1) ) ).T;
    
        B = np.divide(A , np.multiply(S1,S2) );   
        where_are_NaNs = np.isnan(B)
        B[where_are_NaNs] = 0
        # adding a threshold filtering:
        B[b<threshold,:]=0
        B[:,b<threshold]=0        
        A=B 
    
    return A