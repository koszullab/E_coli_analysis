# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:31:33 2017
@author: Axel KournaK
To look for correlation between transcription and 3C contacts. 
"""
from pylab import *
import ice_mirny3
import scn
from itertools import islice
from scipy.stats import gaussian_kde
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

# smooth function:
def smoothListGaussian(list,strippedXs=False,degree=5):  
    window=degree*2-1  
    weight=np.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(np.exp((4*(frac))**2))  
        weightGauss.append(gauss)  
    weight=np.array(weightGauss)*weight  
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  
    return smoothed


data1=loadtxt("/home/axel/Bureau/data/vick/exp5_dec2014_NextSeq/fastq/matrices_5kb_10kb_NextSeq_dec2014/mat_temp_5000_WT.dat")

matscn1=ice_mirny3.ice_func(data1,100)
matscn1= scn.scn_func(matscn1,0)

# plots:
b=np.diagonal(matscn1)
plot(b)
plot(smoothListGaussian(b) )

# transcription data (Rnaseq from Olivier Espeli lab):
chip=loadtxt("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/vick_data_backup/espeli_data/fastq/EV-4_TGACCA_L002_R1_001.fastq.sam.MQ0.hist5000")
chip=loadtxt("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd/vick_data_backup/espeli_data/fastq/hist_EV-4_TGACCA_L002_R1_001.MQ30.hist5000")
plot(log(chip))

# Correlation:
c=log(chip[:,1])

plot( (b-mean(b))/std(b))
plot( (c-mean(c))/std(c))

pearsonr(b,c)

# with smoothing the signals:
pearsonr(smoothListGaussian( (b-mean(b))/std(b)  ), smoothListGaussian( (c-mean(c)/std(c) )  ))
spearmanr(smoothListGaussian( (b-mean(b))/std(b)  ), smoothListGaussian( (c-mean(c)/std(c) )  ))

#  plot of both signal:
plot(smoothListGaussian( (b-mean(b))/std(b) ))
plot(smoothListGaussian( (c-mean(c))/std(c) ))






