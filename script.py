# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:19:48 2016

@author: justin

Use this script to read in mdata files and generate anisotropy table.
I have to remove the first entry in the data for mesh because it was a zero sized bin [0,0]
"""
import glob
from Import_MCNPX_output import *

directory = "/home/justin/I-125 MCNPX/mikell 2micron smesh results"
fileList = glob.glob("{0}/mdata*".format(directory))

#import the mdata smesh results
o = add_in_quadrature(fileList) #answer yes to implicit


#Reproduce the F(r,theta) in Mourtada et al AgX100 TG43 parameters Brachytherapy 2011.
#0.25 is integrated from 0 to 0.5 ( which will be tabulated as 0 degree bin)
thetalist = [0.25,1,2,3,5,7,10,12,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]
dthetalist = np.repeat(0.05,24)
rlist = [0.25,0.5,1,2,3,5,7]
drlist =[0.005,0.005,0.005,0.02,0.02,0.02,0.02]

q = list_F_r_theta(o, rlist, thetalist, drlist, dthetalist)
print(rlist)
for i in np.arange(len(q[0,:])):
    print("{0},{1}".format(thetalist[i],np.array_str(q[:,i], max_line_width=800, precision=3)))

#generate additional F(r,theta) at new requested radii
rlist = [0.25,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,7]
drlist =[0.005,0.005,0.005,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02]
q = list_F_r_theta(o, rlist, thetalist, drlist, dthetalist)
print(rlist)
for i in np.arange(len(q[0,:])):
    print("{0},{1}".format(thetalist[i],np.array_str(q[:,i], max_line_width=800, precision=3)))