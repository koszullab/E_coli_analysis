# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 17:40:37 2015
@author: Axel KournaK
Process of Boccard data to confront 3C and Recombination essays 
"""
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

from scipy.integrate import simps, trapz
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

import scn


#---------------------------------------------------------------------------------
data1 = loadtxt("/home/axel/Bureau/data/vick/2014_05_HiSeq_Boccard/matrices_Vicky/ecoli_MM_async/BC164_GTGT_fused7banks_Ecoli_MM_async/matrices_10kb/mat_temp.dat")
matscn1 = scn.scn_func(data1,0)

length_genome=4639675;
bin_mat = 10000;

bins = [1055417, 3250851, 331520, 3841019, 4024867, 806540, 1379810, 2383627, 2470410];

i=0;
for x_position in bins:
    i=i+1;
    print x_position;
     
    file = "/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd1/ECOLI_PROJECT/CMP_3C_Recombination/boccard_DATA/"+str(x_position)+"_36_20.dat";
    d1 = loadtxt(file);
    d2=d1;

    x_binned = x_position / bin_mat;
    y_binned = d1[:,0]  / bin_mat;
    y_binned = [ int(x) for x in y_binned ];
    d2[:,1] = matscn1[ [x_binned] * len(y_binned) ,y_binned ];
    
    dist = abs(d1[:,0] - x_position); 
    d2[:,0] = dist;
    
    file = "/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd1/ECOLI_PROJECT/CMP_3C_Recombination/boccard_DATA/"+str(x_position)+"_36_20.dat";
    d1 = loadtxt(file);
    d1[:,0] = dist;
       
    if i==1:
        v=np.array(d1);
        v2=np.array(d2);
    else:
        v= np.concatenate((v,np.array(d1)) , axis=0);
        v2= np.concatenate((v2,np.array(d2)) ,axis=0);


dist=v[:,0];
dist[dist > (length_genome / 2) ] = length_genome - dist[dist > (length_genome / 2) ];


#-------------------------------------------------------------------------------------
# Removal of points of recombination < 3%

c2=np.concatenate((  np.vstack( np.array(v[:,1]) ), np.vstack( np.array(v2[:,1] ) )  )  ,axis=1);
c2 =np.array(c2);
c2 = c2[ (c2[:,0] > 3)  ]

plt.yscale('log', nonposy='clip');
plt.xscale('log', nonposy='clip');

plot(c2[:,0],c2[:,1] ,'o');

# LOWESS Fit
#z = lowess(c2[:,1], c2[:,0]);
#plot(z[:,0],z[:,1],lw=3,color="yellow",label="LOWESS Fit")

# linear fit:
slope, intercept = np.polyfit(c2[:,0],c2[:,1] , 1);
plot(c2[:,0], c2[:,0]*slope + intercept, 'r',lw=6,label="linear fit");

plt.xlabel('% of Recombination',fontsize=15);
plt.ylabel('3C score',fontsize=15);

plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)} , loc=1)
plt.title('Comparions between 3C and Recombination methods',fontsize=15);

xticks(fontsize=15)
yticks(fontsize=15)

PCC = pearsonr(c2[:,0],c2[:,1] )

text(28, 0.00065, "PCC = "+"0.80", fontsize=15)
text(28, 0.00050, "p-value = "+"6.32e-27", fontsize=15)
grid()


