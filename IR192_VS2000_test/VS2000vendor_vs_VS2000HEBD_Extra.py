# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:00:46 2022

@author: mikell
"""

import numpy as np
import sys
import matplotlib as mpl
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#below is to import the class jkcm_TG43_calc
sys.path.append(r"H:\src\TG43")
from jkcm_TG43_calc import *

#with the source pointing in the +z direction, 
#perform comparisons from x,y =(0,-15) cm to x,y = (15,15)

#create array of points to calculate dose to
minx_cm=0
maxx_cm=15
dx_cm = 0.05
x_epsilon = 1E-6

xArr = np.arange((maxx_cm-minx_cm)/dx_cm+1)
xArr *= dx_cm
xArr += minx_cm +x_epsilon
nx = len(xArr)

miny_cm=-15
maxy_cm=maxx_cm
dy_cm = dx_cm


yArr = np.arange((maxy_cm-miny_cm)/dy_cm+1)
yArr *= dy_cm
yArr += miny_cm
ny = len(yArr)

zArr = np.zeros(1)

posArr = np.vstack(np.meshgrid(xArr,yArr,zArr)).reshape(3,-1).T

#compare vendor VS2000 with HEBD VS2000
##################################################################
####HEBD consensus data for VS2000################################
rootDir = r'H:\src\TG43'
sourceDir = rootDir + r"\sources\Ir192_VS2000"
frthetafile = sourceDir+r"\Ir192_VS2000_frtheta.txt"
grfile = sourceDir+r"\Ir192_VS2000_gr.txt"
sourcedatafile = sourceDir+r"\Ir192_VS2000_source_data.txt"
o = jkcm_TG43_calc()

o.import_aniso_table(frthetafile)
o.import_gr_table(grfile)
o.import_source_data(sourcedatafile)
#position the source with tip pointing in the +y direction
o.setSourceCenterAndTipPos(0,0,0,0,1,0)

Sk_U = 1.
dwell_h = 1.
q = o.calc_to_points(posArr)
#make any negative values 1: these are likely in the source and result of G(r,theta) breakdown.
index = np.where(q < 0)
q[index] = 1

img_VS2000_HEBD = np.reshape(q, (ny,nx), order='C')


#compare vendor VS2000 with HEBD VS2000
##################################################################
####vendor supplied data for VS2000################################
rootDir = r'H:\src\TG43'
sourceDir = rootDir + r"\sources\Ir192_pwg"
frthetafile = sourceDir+r"\Ir192_pwg_frtheta.txt"
grfile = sourceDir+r"\Ir192_pwg_gr.txt"
sourcedatafile = sourceDir+r"\Ir192_pwg_source_data.txt"
o = jkcm_TG43_calc()

o.import_aniso_table(frthetafile)
o.import_gr_table(grfile)
o.import_source_data(sourcedatafile)
#position the source with tip pointing in the +y direction
o.setSourceCenterAndTipPos(0,0,0,0,1,0)

Sk_U = 1.
dwell_h = 1.
q = o.calc_to_points(posArr)
#make any negative values 1: these are likely in the source and result of G(r,theta) breakdown.
index = np.where(q < 0)
q[index] = 1

img_VS2000_vendor = np.reshape(q, (ny,nx), order='C')

##################################################################
####HEBD data for GMPlus################################
rootDir = r'H:\src\TG43'
sourceDir = rootDir + r"\sources\Ir192_GMPlus"
frthetafile = sourceDir+r"\Ir192_GMPlus_frtheta.txt"
grfile = sourceDir+r"\Ir192_GMPlus_gr.txt"
sourcedatafile = sourceDir+r"\Ir192_GMPlus_source_data.txt"
o = jkcm_TG43_calc()

o.import_aniso_table(frthetafile)
o.import_gr_table(grfile)
o.import_source_data(sourcedatafile)
#position the source with tip pointing in the +y direction
o.setSourceCenterAndTipPos(0,0,0,0,1,0)

Sk_U = 1.
dwell_h = 1.
q = o.calc_to_points(posArr)
#make any negative values 1: these are likely in the source and result of G(r,theta) breakdown at radii less then active source length. 
index = np.where(q < 0)
q[index] = 1

img_GMPlus_HEBD = np.reshape(q, (ny,nx), order='C')

#%matplotlib

#https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
cmap = plt.cm.rainbow
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap_bounds = [-20,-10,-5,-3,-1,0,1,3,5,10,20]
norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)

#compare GMPlus with vendor VS2000
from pylab import *
Z = 100*(img_GMPlus_HEBD/img_VS2000_vendor-1.)
#cmap = cm.get_cmap('bwr', 11)
im = plt.imshow(Z, cmap=cmap, norm=norm)
plt.colorbar(ticks=cmap_bounds, boundaries=cmap_bounds)
plt.title("100*(GMPlus - VS2000)/VS2000")
plt.gca().set_yticklabels([15,15,10,5,0,-5,-10,-15])
plt.gca().set_xticklabels([0,0,5,10,15])
plt.xlabel("x/cm")
plt.ylabel("y/cm")


im = plt.imshow(Z[200:401,0:101], cmap=cmap, norm=norm)
plt.title("100*(GMPlus - VS2000)/VS2000")
plt.gca().set_yticklabels([5,5,3.75,2.5,1.25,0,-1.25,-2.5,-3.75,-5])
plt.gca().set_xticklabels([0,0,2.5,5])
plt.xlabel("x/cm")
plt.ylabel("y/cm")
plt.colorbar(ticks=cmap_bounds, boundaries=cmap_bounds)



