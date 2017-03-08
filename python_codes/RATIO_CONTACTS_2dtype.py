# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:26:01 2016
@author: Axel KournaK
To plot ratio of scalograms 2d type to observe the effects of a certain mutant. 
"""

import numpy as np 
from scipy import stats
import scn
import ice_mirny3
import comp_short
import domainogram_diag
import sliding_genome

from pylab import *
import sys
import os
import shutil
import matplotlib.pyplot as plt
from scalogram_2d_type import *
import scipy.ndimage

data1 = loadtxt(sys.argv[1]);
data2 = loadtxt(sys.argv[2]);
name = sys.argv[3];

dir = name;
if not os.path.exists(dir):
    os.makedirs(dir);
else:
    shutil.rmtree(dir);          
    os.makedirs(dir);

print name;
print "Directory created."

#data1 = loadtxt('/data/vick/exp5_dec2014_NextSeq/fastq/matrices_5kb_10kb_NextSeq_dec2014/mat_temp_5000_WT.dat'); #  wt 
#data2 = loadtxt('/data/vick/exp5_dec2014_NextSeq/fastq/matrices_5kb_10kb_NextSeq_dec2014/mat_temp_5000_MatP.dat');   #  mutant

#data1 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC172_CGGT_MG1655_LB_30C/tmp/mat_LB30C_1000.dat'); #  wt 
#data2 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC78_ACTT_HNS_LB_30C/tmp/mat_HNS_1000.dat');   #  mutant

#data1 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC172_CGGT_MG1655_LB_30C/tmp/mat_ppmm_1000_BC172_CGGT_MG1655_LB_30C.dat');
#data2 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC78_ACTT_HNS_LB_30C/tmp/mat_ppmm_1000_BC78_ACTT_HNS_LB_30C.dat');

#data1 = loadtxt('/data/vick/exp7_juin2015_NextSeq/seqs/BC172_CGGT_MG1655_LB_30C/matrices/mat_temp_MG1655_LB_30C_5kb.dat
#data2 = loadtxt('/data/vick/exp7_juin2015_NextSeq/seqs/BC78_ACTT_HNS_LB_30C/matrices/mat_temp_HNS_5kb.dat');

#data1 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC70_TACT_MG1655_22C/tmp/mat_temp_BC70_TACT_MG1655_22C_resequencing.dat')
#data2 = loadtxt('/home/axel/Bureau/data/vick/exp7_octobre2015_NextSeq-resequencing/BC72_AGCT_MukB_22C/tmp/mat_temp_BC72_AGCT_MukB_22C_resequencing')

matscn1=ice_mirny3.ice_func(data1,100)  #  wt
matscn1= scn.scn_func(matscn1,0)

matscn2=ice_mirny3.ice_func(data2,100); #  mutant
matscn2 = scn.scn_func(matscn2,0)

BIN= 5000

# Choice of a scales:
#scales=range(0,40);

#v=10**(np.array(range(200,640)) *0.01) / BIN;
#v=v.astype(int);
#scales=sorted(set(v));

scales=range(0,200)

extent=(0, 928, 0, len(scales));

s_w=scalo2(matscn1,scales,BIN);  #  wt
s_m=scalo2(matscn2,scales,BIN);  #  mutant

#plots:
s_nat=s_m/s_w;
s=log(s_m/s_w);
s=s[range(s.shape[0]-1,-1,-1),:];
#imshow((s_m/s_w), interpolation='none', extent=extent, cmap="bwr", vmin=0.0, vmax=2.0);
imshow( s, interpolation='none', extent=extent, cmap="bwr", vmin=-1.0, vmax=1.0);

plt.xlabel('Position along the genome (in kb)');
plt.ylabel('Spatial scales (in kb)');
cbar = plt.colorbar(shrink = 0.17,orientation="horizontal", ticks=None);
cbar.ax.tick_params(labelsize=5);

tick_locs = range(len(scales),0,-1);
#v2 = range(0,400,10);
v2=scales[::-1];
tick_lbls = ((array(v2)+1)*5).tolist();

tick_locs = np.array(tick_locs);
tick_locs = tick_locs[range(0,len(tick_locs),8)];
tick_lbls = np.array(tick_lbls);
tick_lbls = tick_lbls [range(0,len(tick_lbls),8)]; 

plt.yticks(tick_locs, tick_lbls,fontsize=5);

tick_locs = range(0,801,200);
tick_lbls = (array( range(0,801,200)) * 5 ).tolist();
plt.xticks(tick_locs, tick_lbls,fontsize=5);
plt.title(name)
n=name+"/"+name+'_scalogram_RATIO_2D_LOG'+'.png';
savefig(n);
close('all');

n1=name+"/"+name+'_Ratio_scalogram2D_LOG'+'.txt';
np.savetxt( n1, s );

n2=name+"/"+name+'_Ratio_scalogram2D_NAT'+'.txt';
np.savetxt( n2, s_nat );

#  with gaussian filter:
s=loadtxt(n1); 
s[np.isnan(s)] = 0.0;
s[np.isinf(s)] = 0.0;

extent = (0, 928, 0, s.shape[0])
img_gaus = scipy.ndimage.filters.gaussian_filter(s, 2, mode='wrap');
imshow(img_gaus, interpolation='none', extent=extent, cmap="bwr", vmin=-1.0, vmax=1.0)
colorbar(shrink = 0.2);
tick_locs = range(len(scales),0,-1);
#v2 = range(0,400,10);
v2=scales[::-1];
tick_lbls = ((array(v2)+1)*5).tolist();

tick_locs = np.array(tick_locs);
tick_locs = tick_locs[range(0,len(tick_locs),8)];
tick_lbls = np.array(tick_lbls);
tick_lbls = tick_lbls [range(0,len(tick_lbls),8)]; 

plt.yticks(tick_locs, tick_lbls,fontsize=5);

tick_locs = range(0,801,200);
tick_lbls = (array( range(0,801,200)) * 5 ).tolist();
plt.xticks(tick_locs, tick_lbls,fontsize=5);
plt.title(name)
n=name+"/"+name+'_scalogram_RATIO_2D_GAUSSIAN_2_LOG'+'.png';
savefig(n);
close('all');

# with Seismic colormap: 
imshow(img_gaus, interpolation='none', extent=extent, cmap="seismic", vmin=-1.0, vmax=1.0);
colorbar(shrink = 0.2);
tick_locs = range(len(scales),0,-1);
v2=scales[::-1];
tick_lbls = ((array(v2)+1)*5).tolist();

tick_locs = np.array(tick_locs);
tick_locs = tick_locs[range(0,len(tick_locs),8)];
tick_lbls = np.array(tick_lbls);
tick_lbls = tick_lbls [range(0,len(tick_lbls),8)]; 

plt.yticks(tick_locs, tick_lbls,fontsize=5);

tick_locs = range(0,801,200);
tick_lbls = (array( range(0,801,200)) * 5 ).tolist();
plt.xticks(tick_locs, tick_lbls,fontsize=5);
n=name+"/"+name+'_scalogram_RATIO_2D_GAUSSIAN_2_LOG_SEISMIC'+'.png';
plt.title(name)
savefig(n);
close('all');




