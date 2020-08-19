# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:58:56 2015
To compute a distance law plot from a Matrice for a specific group of the genome (ex: highly expressed genes)
@author: Axel KournaK
"""

import numpy as np


def dist_law(A,indices):   
    n1 = A.shape[0];
    n2 = A.shape[1];
    print "Number of positions of the matrix entetered:"
    print n1;
    print "Number of scales of the matrix entetered:"
    print n2;
    
    dist = np.zeros((n2/2 +1, 1));
    n_int =   np.zeros((n1, 1));
  
    for nw in range(0,n2/2 +1) :   # scales
        somme = [];         
        for j in range(0,n1) :
                i = indices[j];
                kp= i -nw;
                lp= i + nw;
                if kp < 0 :
                    kp = n2 +kp ;
                if lp < 0 :
                    lp = n2 +lp ; 
                if kp >= n2 :
                    kp = kp - n2;
                if lp >= n2:
                    lp = lp - n2;
                    
                somme.append(A[j,kp]);
                somme.append(A[j,lp]);
                
        dist[nw] = np.mean(somme) ;                         
    return dist;
           
