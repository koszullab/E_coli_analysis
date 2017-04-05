# -*- coding: utf-8 -*-
"""
@author: Axel KournaK
"""
import numpy as np
from pylab import *
import matplotlib.gridspec as gridspec
import scipy
from operator import truediv
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

# mock signal:
m=loadtxt("/home/axel/Bureau/data/HN-S_FIS_ChIP_rpo/mockIP_recovered.fq.sam.MQ0.hist5000")

# chip signals:
chip2 = loadtxt("/home/axel/Bureau/data/HN-S_FIS_ChIP_rpo/hist_chipseq_HNS_WT_coli_5kb.txt");
chip2=chip2[:,1]

chip2 = loadtxt("/home/axel/Bureau/data/HN-S_FIS_ChIP_rpo/hist_chipseq_FIS_WT_EE_5kb.txt");
chip2=chip2

chip2 = loadtxt("/home/axel/Bureau/data/MukB_ChIPseq_unknown_team/ALL_MukB_MQ0.txt.hist5000");
chip2=chip2

# ratio of contacts:
#scalo=loadtxt("/home/axel/Bureau/bacteries_project/z-python_scripts/HNS_LB_30C_rep2_5/HNS_LB_30C_rep2_5_Ratio_scalogram2D_NAT.txt");
#scalo=loadtxt("/home/axel/Bureau/bacteries_project/z-python_scripts/HNS_LB_30C_rep1/HNS_LB_30C_rep1_Ratio_scalogram2D_NAT.txt");

#scalo=loadtxt("/home/axel/Bureau/bacteries_project/z-python_scripts/Fis_4/Fis_4_Ratio_scalogram2D_NAT.txt");
scalo=loadtxt("/home/axel/Bureau/bacteries_project/z-python_scripts/MukBrep2_WTrep2_5/MukBrep2_WTrep2_5_Ratio_scalogram2D_NAT.txt");

#  look for best spatial scale
scalo[np.isnan(scalo)] = 1.0;
scalo[np.isinf(scalo)] = 1.0;

c2=scalo
p_vector = np.zeros(c2.shape[0])
for i in range(0,c2.shape[0]):
    p= spearmanr(chip2, c2[i,:].T)
    p_vector[i] = p[0]

plot(p_vector,linewidth=3.,color="darkred")
plot(p_vector,'o',color="darkred")

# Density plot: 
x=chip2
y=c2[30,:]

pearsonr(x,y)
spearmanr(x,y)

# normal plot
plot( (x-mean(x) )/x.std() ,label="ChIP",linewidth=2.)
plot( (y -mean(y) )/ std(y) ,label="Contacts signal between 0 and 5kb",linewidth=2.);
xlabel("Position along the genome (in bins of 5kb)", fontsize=18);
ylabel(" Z score of the signal", fontsize=18);
xticks(fontsize = 16);
yticks(fontsize = 16);
legend(fontsize=18);
grid();


# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

fig, ax = plt.subplots()
cax = ax.scatter(x, y, c=z, s=50, edgecolor='')
fig.colorbar(cax)
plt.show()

legend()
xlabel("ChIP Fis mutant", fontsize=18);
ylabel("Contacts signal 1 kb", fontsize=18);
xticks(fontsize = 10);
yticks(fontsize = 10);
text(150000,0.23,"Spearman =0.36",fontsize=18)
text(150000,0.13,"pvalue = 2.3e-30",fontsize=18)
grid()

z = lowess(y, x,frac= 0.6);
plot(z[:,0],z[:,1],lw=3,label="LOWESS FIT",color="darkred");
legend(fontsize=18)
