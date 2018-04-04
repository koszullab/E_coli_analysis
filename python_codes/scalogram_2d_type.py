# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:58:14 2016
@author: Axel KournaK
return the mean contacts score at a certain distance along the genome ("type 2 scalograms")
"""
import numpy as np
import comp_short
import comp_short_v4

def scalo2(A,scales,BIN): 
    comp_scale1 = np.zeros(  (len(scales), A.shape[0] ) );
    i=0;
    c_a=0;
    for nw in scales:
        print(nw)
        c=comp_short_v4.comp_func(A,nw);       
        #c=comp_short.comp_func(A,nw);               
        comp_scale1[i,] =  c.T;
        i=i+1; 
    return comp_scale1;
