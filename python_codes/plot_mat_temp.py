# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 16:34:40 2016
@author: Axel KournaK
"""
from pylab import *
import scn
import ice_mirny3
import triangular_diag
import scalogram_2d_type
import matplotlib.gridspec as gridspec
import distance_law
from scipy.stats.stats import pearsonr
import triangular_diag
import domainogram_diag2
import directional_indice
import gc_content
import scipy.ndimage
import scipy.io as sio
import sys

MAT_BIN = loadtxt(sys.argv[1])
name = sys.argv[2]

# Normalisation:
matscn1=ice_mirny3.ice_func(MAT_BIN,0)
matscn1= scn.scn_func(matscn1,0)
name = sys.argv[2]

#  PLOTTING MATRIX:
matplotlib.rcParams.update({'font.size': 8});
subplots_adjust(hspace= 0.5);
plt.figure(num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
#gs = gridspec.GridSpec(1, 1,);

#ax0 = plt.subplot(gs[0]);
C2 = sio.loadmat('/home/axel/Bureau/bacteries_project/z-python_scripts/MyColormaps2.mat');
C = C2['mycmap2'];
cm = mpl.colors.ListedColormap(C);

exponent = 0.2
bmax=0.757858283255199
imshow( matscn1**exponent, cmap=cm,interpolation="none", extent=[0,4640,4640,0], aspect=1,vmin=0,vmax=bmax);
bounds = [0, bmax]
bounds2 = [0, 0.25]
cb1 = colorbar(shrink = 0.2, orientation="horizontal",ticks=bounds, spacing='proportional');
cb1.ax.set_xticklabels(bounds2)
#cb1 = colorbar(shrink = 0.2, orientation="horizontal", spacing='proportional');
cb1.set_label('Contacts Normalised Score');
plt.xlabel('Genome position (in kb)');
plt.title(name);
plt.savefig(name+"_Natural.png",  dpi=600, format='png');

img_gaus = scipy.ndimage.filters.gaussian_filter(matscn1**exponent, 2, mode='wrap')
imshow(img_gaus, cmap=cm,interpolation="none", extent=[0,4640,4640,0], aspect=1,vmin=0,vmax=bmax)
bounds = [0, bmax]
bounds2 = [0, 0.25]
cb1 = colorbar(shrink = 0.2, orientation="horizontal",ticks=bounds, spacing='proportional');
cb1.ax.set_xticklabels(bounds2) 
#cb1 = colorbar(shrink = 0.2, orientation="horizontal", spacing='proportional');
cb1.set_label('Contacts Normalised Score');
plt.xlabel('Genome position (in kb)');
plt.title(name);
plt.savefig(name+"_Natural_Gaussian.eps",  dpi=600, format='eps')


#Vmax= 0.7;  #  Matlab take by default the maximum of the heatmap
#imshow(matscn1**0.2, cmap=cm,interpolation="none", vmin=0.0, vmax=Vmax, extent=[0,4640,4640,0], aspect=1);
ml=log10(matscn1);
ml[np.isinf(ml)] = 0 ;
ml[ml==0]=ml[ml<0].min();
#ml[np.isnan(ml)] = 0.0;
imshow( ml, cmap=cm,interpolation="none", extent=[0,4640,4640,0], aspect=1,vmin=-4,vmax=-0.5)

#TRI1=triangular_diag.triangular_diag(matscn1);
#imshow(TRI1[range( matscn1.shape[0]/2 -1,-1,-1)]**0.2, cmap=cm, interpolation="none",vmin=0.0, vmax=Vmax, aspect=1);

bounds = [-4, -1];
cb1 = colorbar(shrink = 0.2, orientation="horizontal",ticks=bounds, spacing='proportional');
#cb1 = colorbar(shrink = 0.2, orientation="horizontal", spacing='proportional');
cb1.set_label('Contacts Normalised Score log10');
plt.xlabel('Genome position (in kb)');
#plt.ylabel('Genome position (in kb)');

#tick_locs = range(0,801,200);
#tick_lbls = (array( range(0,801,200)) * 5 ).tolist();
#plt.xticks(tick_locs, tick_lbls,fontsize=10);
plt.title(name);

# saving:
plt.savefig(name+"_LOG.eps",  dpi=600, format='eps')
np.savetxt("mat_5000_"+name+"_SCN"+".txt", matscn1)
print(name+" was saved!")
close('all');