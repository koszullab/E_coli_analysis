# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:08:22 2016
@author: Axel KournaK
To confront measures of mobilities of Marco or Espeli with cumulative contacts signals.
"""
import numpy as np
import scipy
import domainogram_diag2
import directional_indice
import scn
import ice_mirny3
import comp_short
from pylab import *
import gc_content
import matplotlib.gridspec as gridspec
import pandas as pd

# Extraction from Marco data:
dm=loadtxt("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/data_marco/fig3_MSD_time_condition.txt")

#dm=loadtxt("/home/axel/Bureau/fig3_MSD_time_condition_resent.txt")
#dm[np.isnan(dm)] = 0.0

dm0=dm[0,]
dm2=np.reshape(dm0, (10,27))

# Computation of 3C cumulative signal:

data1=loadtxt("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/matrices_NextSeqs_june_October2015/mat_temp_5000_fused4banks_30C_MM_ppmm.dat");
#data1=loadtxt("run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/matrices_NextSeqs_june_October2015/mat_temp_1000_fused4banks_30C_MM.dat");

times_laps=loadtxt("/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/data_marco/fig3_time_lags.txt")

BIN=5000;
matscn1=ice_mirny3.ice_func(data1,100);
matscn1= scn.scn_func(matscn1,0);
n1 = matscn1.shape[0];

ii=0;
scales=range(0,400);
comp_scale1 = np.zeros(  (len(scales), n1) );
for nw in scales :
    print(nw);
    sc = nw*BIN
    c=comp_short.comp_func(matscn1,nw);
    comp_scale1[ii,] =  c.T;
    ii=ii+1

imshow(comp_scale1,interpolation="none")

# MSD measures:
msd=loadtxt('/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/data_marco/fig1b_MSD_10s.txt')
df=pd.read_table('/run/media/axel/RSG3/ECOLI_PROJECT_IN_PROGRESS/data_marco/names_positions_probes.txt',header=None)
v=np.array(df[1]);

# look for the points in the scalogram cumulative:
scale_chosen=200000;
ii=0;
#comp_local = zeros((msd.shape[0],0));
comp_local = [];
indices = [];
names=[];
for i in msd[:,0] *1000/5 :
    ind = int(i+0.5); 
    comp_local.append(comp_scale1[scale_chosen/BIN,ind]);
    indices.append(ind);
    #names.append();
    ii=ii+1;

# Correlation between both signals:
c = scipy.stats.spearmanr(comp_local,msd[:,1]);
#plot(comp_local,msd[:,1],'o');
 
plot(comp_local,msd[:,1],'o', markersize=8.0);

#ylabel("MSD with tau = 10 sec",fontsize=16);
#xlabel("Cumulative 3C Contacts signal below 200kb",fontsize=16);
text(0.6,0.014,"Spearman Coef = -0.85");
text(0.6,0.012,"p-value = 1.77e-08");

# Linear regression:
x = comp_local;
y = msd[:,1];
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y);
x=np.array(x);
plot(x, x*slope + intercept, 'r',linewidth = 4.0);
xticks(fontsize = 16);
yticks(fontsize = 16);

na=["Ter1","NSR2","NSR1","NSR4","OriC","Ori1","Ter6","Ori2","Ori3","Ori4","NSR3","NSR5","Right1","Right2","Right4","Ter2","Ter5","Ter7","Left1","NSL3","NSL2","NSL1","Left3","Ter8","Left2","Ter3","Ter4"]
for i in range(0,len(comp_local) ) :
    text(comp_local[i]+0.007,msd[i,1], na[i],fontsize = 14) 



#-----------------------------------------------------------------------------------------
# Computation of matrices of scales
# Spatial scales versus temporal scales 

#scales = range(0,400);
#scales = (100,200,300)
#for i in range(321,322) :
#    for j in range(6,10)  :  

#for i in (50,50) :
#    for j in (8,8)  : 


MAT_SCALES = np.zeros(  (len(scales) ,  10)  );
for i in scales:
     for j in range(0,10):
        scale_spatial = i;
        scale_temporal = j;
        
        comp_local_i = comp_scale1[scale_spatial, indices ];
        #msd_j = msd[:,1]
        msd_j = msd[ scale_temporal  , range(1,msd.shape[1])  ] * 4 * msd[ scale_temporal ,0 ];
        
        z_comp_local_i = scipy.stats.mstats.zscore(comp_local_i);        
        z_msd_j        = scipy.stats.mstats.zscore( msd_j);
        
        x= comp_scale1[scale_spatial, sort(indices) ]
        y= dm2[scale_temporal,:]
        
        #plot(x, y,'o');
        #xlabel("Local Compaction from 3C");
        #ylabel("MSDXY / sum ");        
        
        c= scipy.stats.pearsonr( x, y );
        #c= scipy.stats.spearmanr( x, y );
        #c_rounded=c[0];
        #c_rounded=round(c_rounded,2);
        #print i,c[0],c[1];

        #slope, intercept, r_value, p_value, std_err = stats.linregress(x,y);
        #plot(x, x*slope + intercept, 'r',linewidth = 2.0);
        #text(0.9*max(x), 0.9*max(y), "Spearman Coefficient = "+ str(c_rounded)  );
        #text(0.9*max(x), 0.8*max(y), "p-value = "+ str(c[1])  );        
        #print i,j,abs(c);
        MAT_SCALES[i,j] = - c[0]
#        if MAT_SCALES[i,j]  > 0.75 :
#            print i,j,MAT_SCALES[i,j];
  
imshow(MAT_SCALES,interpolation="none",extent= (0, 100, 400, 0))


imshow(MAT_SCALES[:,times_laps.argsort()],interpolation="none",extent= (0, 100, 400, 0))
ylabel("Spatial scales")
xlabel(r'Temporal scales $\tau$ in sec')
cbar = colorbar()
cbar.set_label('Spearman Correlation Coefficient', rotation=270, fontsize = 14)
plt.xticks(rotation=70)

# Plot with one temporal scale:

plot(MAT_SCALES[:,0],linewidth=3.0)
ylabel("Spearman Correlation Coefficient"+"\n" +"between MSD measured at 10.00 sec (Javer et al. data)" +"\n"+" and 3C cumulative signal (MM, 30($^\circ$C)",fontsize=14)
xlabel("Cumulative 3C contact signal below a specific spatial scale (in kb)",fontsize=14)
plt.axvline(40,linewidth=3.,color="red")
grid()

      
# PLOT of matrix of scales:               
gs = gridspec.GridSpec(1,1);
ax1 = plt.subplot(gs[0]);
extent=(0,1100,0,400);        
ax1.imshow(MAT_SCALES, interpolation='none', extent=extent);
#contourf(MAT_SCALES, interpolation='none', extent=extent);
im = imshow(MAT_SCALES, interpolation='none', extent=extent);
#imshow(MAT_SCALES, interpolation='none');
ax1.set_ylabel("Spatial scales used in 3C (in kb)");
ax1.set_xlabel("Temporal scales used in microscopy (in sec)")

tick_locs = list( np.array(range(0,msd.shape[0] ) ) *100 +50 ) ;
tick_lbls = msd[range(1,msd.shape[0] ),0];
ax1.set_xticks(tick_locs);
ax1.set_xticklabels(tick_lbls, fontdict=None, minor=False);

tick_locs = range(390,0,-20);
tick_locs = range(10,390,20);
tick_lbls = ((array(tick_locs))*5).tolist();
tick_locs = range(390,0,-20);
ax1.set_yticks(tick_locs);
ax1.set_yticklabels(tick_lbls, fontdict=None, minor=False);
plt.colorbar(im,shrink = 0.2,orientation="vertical");

