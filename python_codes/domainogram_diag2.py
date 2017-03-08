# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:38:48 2014
@author: axel
This prog calculates domainograms for matrices 
2: add of a distance threshold to remove the effects of the inner squares
"""

import numpy as np 
#  DOMAINOGRAM FUNCTION : take a matrice and calculate the mean interaction 
#  inside a square along the main diagonal
#  Input parameters :
#  A : a matrix of HiC (better if it is a correlation matrix)
#  nw :  number of bins to make the sum in the square    
#  diag : only take the sum of the tranverse diag  (can be better and more rapid)
    
    
def dom_diag(A, nw):
    n1 = A.shape[0];
    th= 0;  # a threshold for a minimal distance above with calculate the signal 
    
    "Size of the matrix entetered for the domainogram diag :"
    print(n1)
    
    somme =   np.zeros((n1, 1));
    signal1 = np.zeros((n1, 1));
    n_int =   np.zeros((n1, 1));
    
    for i in range(0,n1) :    
        for k in range(-nw,nw+1) :
               kp =i-k; 
               lp =i+k;
               # Circularity conditions: 
               if kp < 0 :
                    kp = n1 +kp ;
               if lp < 0 :
                    lp = n1 +lp ;
               if kp >= n1 :
                   kp = kp - n1;
               if lp >= n1:
                   lp = lp - n1;
               # If we want to put a threshold above with, the computation is done:     
               if (kp >= (i + th) or kp <= (i - th)  and (lp >= (i + th) or lp <= (i - th) ) ):
                   somme[i] = somme[i] + A[kp,lp];
                   n_int[i] = n_int[i] + 1;
               
    signal1 = somme / n_int;
    return signal1