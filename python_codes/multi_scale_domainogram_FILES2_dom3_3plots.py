# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:39:30 2015
@author: Axel KournaK
2 DI + SCALOGRAM  
"""
import numpy as np
import domainogram_diag2
import directional_indice
import scn
import ice_mirny3
import comp_short
from pylab import *
import gc_content
import matplotlib.gridspec as gridspec

BIN = 5000
data1 = loadtxt(sys.argv[1]);

matscn1=ice_mirny3.ice_func(data1,100);
matscn1= scn.scn_func(matscn1,0);

#matscn1=data1;
name = sys.argv[2];

#  PLOTTING MATRIX AND MULTI-SCALES DOMAINOGRAM:
M = np.corrcoef(matscn1);
#M = matscn1;
n1 = M.shape[0];
# Choice of scales:
#v=10**(np.array(range(50,161)) *0.04) / BIN
#v=v.astype(int);    
#v2=sorted(set(v));
#scales=v2;
scales = range(0,100);
#---------------------------------------------------------------------------------
DOM = np.zeros((len(scales), n1));
DI = np.zeros((len(scales), n1));

ii=0;
for i in scales :
    DOM[ii,:] = domainogram_diag2.dom_diag(M,i).T
    DI[ii,:] = directional_indice.directional(M,i).T
    ii=ii+1;

#--------------------------------------------------------------------------------
ii=0;
comp_scale1 = np.zeros(  (len(scales), n1) );
for nw in scales:
    print nw;
    sc= nw*BIN;
    c=comp_short.comp_func(matscn1,nw);
    comp_scale1[ii,] =  c.T;
    ii=ii+1;

#------------------------------------------------------------
#  PLOTTING MATRIX AND MULTI-SCALES DOMAINOGRAM:

gs = gridspec.GridSpec(3, 1);
matplotlib.rcParams.update({'font.size': 8});
subplots_adjust(hspace= 0);
plt.figure(num=None, figsize=(6.52, 6.52), dpi=80, facecolor='w', edgecolor='k');

#ax0 = plt.subplot(gs[0]);
#ax0.imshow(DOM[range( len(scales) -1,-1,-1)],interpolation='none',cmap='PuOr',extent=[0,matscn1.shape[0],0,200]);
#ax0.contourf(DOM,cmap='PuOr',extent=[0,matscn1.shape[0],0,200]);
#ax0.set_ylabel("Scales (in kb)");

scales2 = range(0,len(scales)*2,20);
tick_locs = scales2;
tick_lbls = ((array(scales2))*5).tolist();
plt.yticks(tick_locs, tick_lbls,fontsize=10);

tick_locs = range(0,matscn1.shape[0]-100,200);
tick_lbls = (array( range(0,matscn1.shape[0]-100,200)) * 5 ).tolist();
plt.xticks(tick_locs, tick_lbls,fontsize=10);

DI80=DI[80,:];
ax1 = plt.subplot(gs[0])
b1=0;
b2=len(DI80.T);
borders = list(range(b1,b2));
borders2 = array(DI80.T);
#borders2 = array(borders2[:,0]);
ax1.set_xlim([b1, b2]);
ax1.set_ylim([-2.0, 2.0]);
ax1.fill_between( borders, 0,borders2, color='red');
ax1.fill_between( borders, 0,borders2,borders2<0 ,color='green' );
ax1.set_ylabel("DI (scale = 400kb)");


DI20=DI[20,:];

ax2 = plt.subplot(gs[1],sharex=ax1);
b1=0;
b2=len(DI20.T);
borders = list(range(b1,b2));
borders2 = array(DI20.T);
#borders2 = array(borders2[:,0]);
ax2.set_xlim([b1, b2]);
ax2.set_ylim([-2.0, 2.0]);
ax2.fill_between( borders, 0,borders2, color='red' );
ax2.fill_between( borders, 0,borders2,borders2<0 ,color='green' );
ax2.set_ylabel("DI scale = 100kb)");

ax3 = plt.subplot(gs[2],sharex=ax1);
#ax4.imshow(comp_scale1,interpolation='none', vmin=0.0, vmax= 0.9);
ax3.contourf(  comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matscn1.shape[0],0,400],cmap="rainbow");
im = contourf( comp_scale1, vmin=0.0, vmax= 1.0,levels=[0, .15, .30, .45, 0.6,0.75, 1.0],extent=[0,matscn1.shape[0],0,400],cmap="rainbow");
ax3.set_ylabel("Scales (in kb)");
#ax3.set_xlabel("Position along the genome (in kb)");

#tick_locs = range(515,00,-10);
scales2 = range(0,len(scales)*4,40);
tick_locs = scales2;
tick_lbls = ((array(scales2))*5/2).tolist();
plt.yticks(tick_locs, tick_lbls,fontsize=10);

tick_locs = range(0,matscn1.shape[0]-100,200);
tick_lbls = (array( range(0,matscn1.shape[0]-100,200)) * 5 ).tolist();
plt.xticks(tick_locs, tick_lbls,fontsize=10);
bounds = [0., 1.0];
cb1 = colorbar(im,shrink = 0.2, orientation="horizontal",ticks=bounds, spacing='proportional');


# saving:
plt.savefig(name+"_DOM"+".svg",  dpi=600, format='svg');
plt.savefig(name+"_DOM"+".eps",  dpi=600, format='eps'); 
plt.savefig(name+"_DOM"+".jpeg", dpi=600, format='jpeg');
print name+" was saved!";
close('all');
