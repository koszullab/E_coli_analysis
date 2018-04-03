# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:38:48 2014
@author: Axel KournaK
This function convert matrices into triangular matrices 
to facilitate comparison with other signals.
"""

import numpy as np 

def triangular_diag(A):

    n1 = A.shape[0];
    print("Size of the matrix entetered for the triangular representation:")
    print(n1)
    scales = range(0,n1/2);
    TRI = np.zeros((len(scales), n1));
    
    for i in range(0,n1) :
        for k in scales :
               kp =i-k;
               lp =i+k;
               # Circularity conditions: 
               if kp < 0 :
                    kp = n1 +kp ; TRI[k,i] = 'NaN';
               elif lp < 0 :
                    lp = n1 +lp ; TRI[k,i] = 'NaN';
               elif kp >= n1 :
                   kp = kp - n1;  TRI[k,i] = 'NaN';
               elif lp >= n1:
                   lp = lp - n1;  TRI[k,i] = 'NaN';
               else :
                   TRI[k,i] = A[kp,lp];           
    return TRI

# For plotting then:
# scales = range(0,500)
#TRI1=triangular_diag(matscn1);
#imshow(TRI1[range( (len(scales)) -1,-1,-1)]**0.2, interpolation="none")
#imshow(TRI1[range( (len(scales)) -1,-1,-1)]**0.15, cmap=cm,interpolation="none");
