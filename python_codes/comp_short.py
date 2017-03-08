# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:10:16 2015

@author: Axel
"""
import numpy as np 

#  COMP FUNCTION :
    
def comp_func(A, nw):    
    n1 = A.shape[0];
    th= 0;  # a threshold for a minimal distance above with calculate the signal 
    #print "Size of the matrix entetered for the Compaction signal :"
    #print n1;
    
    somme_short = np.zeros((n1, 1));
    somme_long = np.zeros((n1, 1));    
    signal1 = np.zeros((n1, 1));
    n_int =   np.zeros((n1, 1));
    mat_test = np.zeros((n1, n1));
    
    for i in range(0,n1) :  
        if i<=nw:
            p1=n1 + i-nw;
            p2=i + nw;
            for k in range(0,n1) :
                if k<=p2 or k>=p1:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;

        elif  (n1 - nw) <= i:
            p1=i- nw;
            p2=i+ nw-n1;
            for k in range(0,n1) :
                if k<=p2 or k>=p1:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;

        else :
            p1= i - nw;
            p2= i + nw;
            for k in range(0,n1) :
                if p1<=k and k<=p2:
                    somme_short[i] = somme_short[i] + A[i,k];
                    n_int[i] = n_int[i] +1;
      
    signal1 = somme_short ;                       
    return signal1;
       
                    