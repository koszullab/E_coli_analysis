# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 18:11:18 2016
@author: Axel KournaK
"""
import numpy as np
from scipy.stats.stats import pearsonr

st=loadtxt('/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/3D_structures/3Dcoor_ecoli_5kb_3_3_copie.xyz');

file = "/run/media/axel/RSG3/EColi_backup/data_olivier_bio_cell/genes_paires.txt2";
d1 = loadtxt(file);
#plt.plot(d1[:,3],d1[:,2],'o');

d2 = loadtxt('/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/IF_duigou/if_duigou.txt');

BIN=5000;
l_points=[];

for i in range(0,d1.shape[0]) :
    b1=int(d1[i,0]/BIN);
    b2=int(d1[i,1]/BIN);
    dist = sqrt(  (st[b2,0]-st[b1,0])**2 + (st[b2,1]-st[b1,1])**2 + (st[b2,2]-st[b1,2])**2 );
    #plot(d1[i,2],dist,'o');
    l_points.append(dist);
    
for i in range(0,d2.shape[0]) :
    b1=int(d2[i,0]/BIN);
    b2=int(d2[i,1]/BIN);
    dist = sqrt(  (st[b2,0]-st[b1,0])**2 + (st[b2,1]-st[b1,1])**2 + (st[b2,2]-st[b1,2])**2 );
    #plot(d1[i,2],dist,'o');
    l_points.append(dist);

d= np.concatenate(( d1[:,2], d2[:,2]), axis=0);
pearsonr(d, np.array(l_points));
slope, intercept = np.polyfit(d,l_points , 1);

# plot:     
plot(d, np.array(l_points),'o',markersize=10.0,color="steelblue");
plot(d, d*slope +intercept, 'r',linewidth=5.0);

text(0.4, 0, "Pearson Coef= 0.86", fontsize=18);
text(0.4, -3, "P-value= 1.16e-10", fontsize=18);

xlabel("Measured distances (in um)",fontsize = 16);
ylabel("Measured distances 3C structure (a.u)",fontsize = 16);

xticks(fontsize = 16);
yticks(fontsize = 16);

savefig("cmp_FROS_3C.pdf"); 