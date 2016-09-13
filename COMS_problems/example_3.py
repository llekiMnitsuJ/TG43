# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:55:30 2016

@author: jusmikel
"""

#import the libraries
#runfile('G:/data/src/TG43/jkcm_TG43_calc.py', wdir='G:/data/src/TG43')
#runfile('G:/data/src/TG43/jkcm_samemodel_multisource_TG43.py', wdir='G:/data/src/TG43')

#location of your TG43 data
rootDir="G:\data\src\TG43"
sourceDir=rootDir+"\sources\I125A"
frthetafile = sourceDir+"\I125A_frtheta.txt"
grfile = sourceDir+"\I125A_gr.txt"
sourcedatafile = sourceDir+"\I125A_source_data.txt"
 
#create object for calculation
o = jkcm_samemodel_multisource_TG43()
#initialize object with TG43 data
o.initializeTG43tables(frthetafile,grfile,sourcedatafile)
 
#display source information
print(o.jkcm_TG43_calc_obj)
 
#file containing source locations (problem specified 10x10mm tumor, so add 1 mm margin and use 12x12COMS)
source_positions_filename = rootDir+"\COMS_plaques\COMS_12mm_plaque.txt"
 
#import source positions
o.importSources(source_positions_filename)

#calc at point of interest (tumor apex)
dose_cGy_arr_per_Sk_per_effT = o.calc_at_point([0,0,0.28])
q = o.jkcm_TG43_calc_obj
effT_h = q.calc_eff_time(100)
dose_Gy_per_Sk = np.sum(dose_cGy_arr_per_Sk_per_effT)*effT_h/100
Sk = 85/dose_Gy_per_Sk

#calc other points of interest
es_dose_Gy = np.sum(o.calc_at_point([0,0,-0.1]))*Sk*effT_h/100

is_dose_Gy = np.sum(o.calc_at_point([0,0,0]))*Sk*effT_h/100

c5_dose_Gy = np.sum(o.calc_at_point([0,0,0.5]))*Sk*effT_h/100

ta_dose_Gy = np.sum(o.calc_at_point([0,0,0.28]))*Sk*effT_h/100

eo_dose_Gy = np.sum(o.calc_at_point([0,0,1.1]))*Sk*effT_h/100

or_dose_Gy = np.sum(o.calc_at_point([0,0,2.2]))*Sk*effT_h/100

#return answers
print("Sk for each seed = {0} U".format(Sk))
print("External Sclera = {0} Gy".format(es_dose_Gy))
print("Internal Sclera = {0} Gy".format(is_dose_Gy))
print("COMS 5 mm = {0} Gy".format(c5_dose_Gy))
print("Tumor apex = {0} Gy".format(ta_dose_Gy))
print("Eye Origin = {0} Gy".format(eo_dose_Gy))
print("Opposite Retina = {0} Gy".format(or_dose_Gy))