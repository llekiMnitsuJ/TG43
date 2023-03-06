# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:11:52 2022

@author: mikell
"""

import numpy as np
import sys
#below is to import the class jkcm_TG43_calc
sys.path.append(r"H:\src\TG43")
                
from jkcm_TG43_calc import *

###INSTRUCTIONS
# ensure you have the correct source files:
    # GMPlus HEBD:
    #         Ir192_GMPlus_frtheta.txt
    #         Ir192_GMPlus_gr.txt
    #         Ir192_GMPlus_source_data.txt
# 1) load the jkcm_TG43_calc class first. 
# 2) adjust the rootDir and sourceDir if needed. 
# 3) then run this script. 
# 4) the following files will be created:
#       "GMPlus-single-source-HEBD-annual-qa_JKM.csv"
#       "GMPLus-dual-source-HEBD-annual-qa_JKM.csv"

##################################################################
####single source calculation for RadCalc commissioning GMPlus######
##################################################################
####HEBD consensus data for GMPlus##########################

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

yposArray = np.array(
    [
     [0,0.5,0],
     [0,1,0],
     [0,2,0],
     [0,5,0],
     [0,8,0],
     [0,1.5,0],
     [0,3,0]]
    )

xposArray = np.array(
    [
     [0.5,0,0],
     [1,0,0],
     [2,0,0],
     [5,0,0],
     [8,0,0],
     [1.5,0,0],
     [3,0,0]]
    )

#match excel printout with RadCalc
Sk_U = 40700.
dwell_h = 100./3600.
s = o.generate_handcalc_paramaters_table_for_points(xposArray, Sk_U, dwell_h)
yposArray = np.array(
    [
     [0,0.5,0],
     [0,1,0],
     [0,2,0],
     [0,5,0],
     [0,8,0]]
    )
s += "\n" + o.generate_handcalc_paramaters_table_for_points(yposArray, Sk_U, dwell_h, header=False)
s += "\n" + o.generate_handcalc_paramaters_table_for_points(-1*yposArray, Sk_U, dwell_h, header=False)
xposArray = np.array(
    [
     [0.5,0,0],
     [1,0,0],
     [2,0,0],
     [5,0,0],
     [8,0,0]]
    )
s += "\n" + o.generate_handcalc_paramaters_table_for_points(yposArray+xposArray, Sk_U, dwell_h, header=False)
s += "\n" + o.generate_handcalc_paramaters_table_for_points(-1*yposArray+xposArray, Sk_U, dwell_h, header=False)
print(s)
f =  open("GMPlus-single-source-HEBD-annual-qa_JKM.csv", 'w')
f.write(s)
f.close()


##################################################################
####dual source calculation for RadCalc commissioning GMPLus######
##################################################################
####HEBD consensus data for GMPLus##########################

sourceDir = rootDir + r"\sources\Ir192_GMPlus"
frthetafile = sourceDir+r"\Ir192_GMPlus_frtheta.txt"
grfile = sourceDir+r"\Ir192_GMPlus_gr.txt"
sourcedatafile = sourceDir+r"\Ir192_GMPlus_source_data.txt"

o = jkcm_TG43_calc()

o.import_aniso_table(frthetafile)
o.import_gr_table(grfile)
o.import_source_data(sourcedatafile)
#position the source at y=0.5 with tip pointing in the +y direction
o.setSourceCenterAndTipPos(0,0.5,0,0,1.5,0)

xposArray = np.array(
    [
     [0.5,0,0],
     [1,0,0],
     [2,0,0],
     [3,0,0],
     [4,0,0],
     [6,0,0]]
    )

#match excel printout with RadCalc
Sk_U = 40700.
dwell_h = 100./3600.
s = o.generate_handcalc_paramaters_table_for_points(xposArray, Sk_U, dwell_h)

#change position and recalc again
#position the source at y=-0.5 with tip pointing in the +y direction
o.setSourceCenterAndTipPos(0,-0.5,0,0,1.5,0)
s += "\n" + o.generate_handcalc_paramaters_table_for_points(xposArray, Sk_U, dwell_h, header=False)

print(s)
f =  open("GMPLus-dual-source-HEBD-annual-qa_JKM.csv", 'w')
f.write(s)
f.close()
