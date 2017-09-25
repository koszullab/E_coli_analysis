# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:32:44 2017
@author: Axel Kournak
To plot SUM SUB MAT of a group of positions like matS sites or TSS.
"""
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Restriction import *
import random
import time
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
import matplotlib.gridspec as gridspec
import scn 
import ice_mirny3
import scipy.ndimage
import scipy.io as sio
import distance_law_human

C2 = sio.loadmat('/home/axel/Bureau/bacteries_project/z-python_scripts/MyColormaps2.mat');
C = C2['mycmap2']
cm = mpl.colors.ListedColormap(C)

# Hiseq data:
MATRICE1 = loadtxt("/home/axel/Bureau/data/vick/2014_05_HiSeq_Boccard/matrices_Vicky/ecoli_MM_async/BC164_GTGT_fused7banks_Ecoli_MM_async/matrices_2kb/mat_temp.dat")
MATRICE2 = loadtxt("/home/axel/Bureau/data/vick/2014_05_HiSeq_Boccard/matrices_Vicky/ecoli_MM_deltaMatP/BC172_CGGT_fused7banks_ecoli_MM_deltaMatP/matrices_DeltaMatP/matrices_2kb/mat_temp.dat")      

# Next seq data:
MATRICE1 = loadtxt("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/vick_data_backup/exp5_dec2014_NextSeq/fastq/BC70_TACT_WT/matrices_BC70_TACT_WT/mat_temp.dat")
MATRICE2 = loadtxt("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/vick_data_backup/exp5_dec2014_NextSeq/fastq/BC72_AGCT_MatP/matrices_BC72_AGCT_MatP/mat_temp.dat")      

# Normalisations: 
matscn1=ice_mirny3.ice_func(MATRICE1,100)
matscn1 = scn.scn_func(matscn1,0)

matscn2=ice_mirny3.ice_func(MATRICE2,0)
matscn2 = scn.scn_func(matscn2,0)

# Plots:
imshow(matscn1**0.2,interpolation="none",cmap=cm,vmin=0.0,vmax=matscn1.max()**0.2)
colorbar()
plot(mats/BIN,mats/BIN,"o",color="orange")

# Ratio:
r=matscn2/matscn1
r = log(r)
r[np.isnan(r)] = 0.0
r[np.isinf(r)] = 0.0

imshow(r, interpolation='none', cmap="seismic", vmin=-1.0, vmax= 1.0)
xlabel("Position along the genome (2 kb bins)")
ylabel("Position along the genome (2 kb bins)")
colorbar()

img_gaus = scipy.ndimage.filters.gaussian_filter(r, 2, mode='wrap')
imshow(img_gaus, interpolation='none', cmap="seismic", vmin=-1.0, vmax=1.0)

df=pd.read_table('/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/vick_data_backup/matS_cell.dat2',header=None, delimiter=" ")
tss = loadtxt("/home/axel/Bureau/TSS_ECOLI/coli/genes_K12_UCSC.dat.tss")
df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/higly_epxressed_genes5.txt',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/SRP_genes.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/SECB_genes.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/poorly_expressed_genes.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/SRP_genes_without_10pc_higly_expressed.txt',header=None, delimiter=" ")

df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/heEPOD.txt',header=None, delimiter=" ")
df=pd.read_table('/home/axel/Bureau/TSS_ECOLI/tsEPOD.txt',header=None, delimiter=" ")

# Along the diagonale :
BIN = 2000
n1 = shape(matscn1)[0]
mats = df[0]
area = 20  #  number of bins around the mats sites
MAT_SUM = np.zeros(   (area*2+1,area*2+1)  )
ns =0
for site in mats :
    print(site)
    site = int(site/BIN)
    pi =0
    ns +=1
    for i in range(site-area,site+area+1) :
        pj =0
        for j in range(site-area,site+area+1) :
            if i >=0 and j>=0 and i<n1 and j<n1 :
                MAT_SUM[pi,pj] = MAT_SUM[pi,pj] + matscn1[i,j]
            else :
                MAT_SUM[pi,pj] = MAT_SUM[pi,pj] + 0
            pj +=1
        pi +=1    
print(ns)

imshow( (MAT_SUM/ns)**0.2,interpolation="none")
imshow( (MAT_SUM/ns)**0.2,interpolation="none",cmap=cm,vmin=1.5,vmax=2.6)

imshow( (MAT_SUM/ns)**0.2,interpolation="none",vmin=1.8,vmax=3.)
#imshow( (MAT_SUM/ns)**0.2,interpolation="none",vmin=1.8,vmax=3.)
tick_locs=(0,20,40)
tick_lbls=('-40 kb','tsEPOD (from Vora et al.)','+ 40kb')
colorbar()
xticks(tick_locs, tick_lbls,fontsize=15)

imshow( (MAT_SUM/ns)**0.2,interpolation="none")
tick_locs=(0,10,20)
tick_lbls=('-20 kb','HEG','+ 20kb')
colorbar()
xticks(tick_locs, tick_lbls,fontsize=15)
yticks(tick_locs, tick_lbls,fontsize=15)

# Between sites (removing the ones on the diagonale):
random_area = 10000
BIN = 2000
mats1 = df[0] + random_area
mats2 = df[0] + random_area
area = 10  #  number of bins around the mats sites
MAT_SUM = np.zeros(   (area*2+1,area*2+1)  )
ns =0
for site1 in mats1 :
    site1 = int(site1/BIN)
    for site2 in mats2 :
        site2 = int(site2/BIN)
        if site1 != site2:
            pi =0
            ns +=1
            for i in range(site1-area,site1+area+1) :
                pj =0
                for j in range(site2-area,site2+area+1) :
                    if i >=0 and j>=0 and i<n1 and j<n1 :
                        MAT_SUM[pi,pj] = MAT_SUM[pi,pj] + matscn2[i,j]
                    else :
                        MAT_SUM[pi,pj] = MAT_SUM[pi,pj] + 0
                    pj +=1
                pi +=1    
print(ns)

imshow( (MAT_SUM/ns)**0.2,interpolation="none",vmin=1,vmax=3.)

tick_locs=(0,10,20)
tick_lbls=('-20 kb','Random genes (HEG 480kb shifted)','+ 20kb')
colorbar()
xticks(tick_locs, tick_lbls,fontsize=15)
yticks(tick_locs, tick_lbls,fontsize=15)


# Plots: 
imshow( (MAT_SUM/ns)**0.2,interpolation="none",cmap=cm,vmin=0.,vmax=0.72)
tick_locs=(0,10,20)
tick_lbls=('-20 kb','matS','+ 20kb')
xticks(tick_locs, tick_lbls,fontsize=15)
yticks(tick_locs, tick_lbls,fontsize=15)
title("MatP NextSeq 2kb bins (between matS sites)")
colorbar()

# ratio plot: 
imshow(log(MAT_SUM2/MAT_SUM1),interpolation="none",vmin=-1,vmax=1.,cmap="bwr")
tick_locs=(0,10,20)
tick_lbls=('-20 kb','random site','+ 20kb')
xticks(tick_locs, tick_lbls,fontsize=15)
yticks(tick_locs, tick_lbls,fontsize=15)
colorbar()
title("MatP/WT NextSeq 2kb bins (between sites)")




# Detrentage to the genomic distance effect:
d=distance_law_human.dist_law(matscn1)
n1=matscn1.shape[0]
## Computation of genomic distance law matrice:
MAT_DIST =  np.zeros((n1, n1));
for i in range(0,n1) :
    for j in range(0,n1) :
        MAT_DIST[i,j] =  d[abs(j-i)]  
        
    ##imshow(MAT_DIST**0.2)
    ## detrendage:
MAT2=matscn1/MAT_DIST
MAT2[np.isnan(MAT2)] = 1.0



r=log(matscn2/matscn1)
imshow(log(r),interpolation="none",cmap="bwr",vmin=-1.0,vmax=1.0)

