# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:58:56 2015
To compute a distance law plot from a Matrix.
@author: Axel
"""

import numpy as np

#  DIST FUNCTION :

def dist_law(A):
    
    n1 = A.shape[0];
    th= 0;  # a threshold for a minimal distance above with calculate the signal 
    print("Size of the matrix entetered for the Compaction signal :")
    print(n1)
    
    dist = np.zeros((n1/2 +1, 1));
    n_int =   np.zeros((n1, 1));
  
    for nw in range(0,n1/2 +1) :   # scales
        somme = [];         
        for i in range(0,n1) :    #  over the whole matrix
                kp= i - nw;
                lp= i + nw;
                if kp < 0 :
                    kp = n1 +kp ;
                if lp < 0 :
                    lp = n1 +lp ; 
                if kp >= n1 :
                    kp = kp - n1;
                if lp >= n1:
                    lp = lp - n1;
                    
                somme.append(A[i,kp]);
                somme.append(A[i,lp]);
                
        dist[nw] = np.mean(somme)                      
    return dist;
           
