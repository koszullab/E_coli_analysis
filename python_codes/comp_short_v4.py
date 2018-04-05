# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:10:16 2015
@author: Axel KournaK 
v2: Other way of Computation taking 3 identical matrices 
"""

import numpy as np 
#  COMP FUNCTION :
    
def comp_func(A, nw):
    n1 = A.shape[0];
    B = np.concatenate((A,A,A), axis=1);

    somme_short = np.zeros((n1, 1));  
    n_int = np.zeros((n1, 1)); 
    
    for i in range(0,n1) :
        for k in (n1+i-nw,n1+i+nw) :
            somme_short[i] = somme_short[i] + B[i,k];
            n_int[i] = n_int[i] +1;
                             
    return somme_short
       
                    

