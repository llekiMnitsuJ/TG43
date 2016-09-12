# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 18:11:01 2016

@author: jusmikel
"""

#import the libraries
#runfile('G:/data/src/TG43/jkcm_TG43_calc.py', wdir='G:/data/src/TG43')
#runfile('G:/data/src/TG43/jkcm_samemodel_multisource_TG43.py', wdir='G:/data/src/TG43')

#location of your TG43 data
rootDir="G:\data\src\TG43\sources\I125A"
frthetafile = rootDir+"\I125A_frtheta.txt"
grfile = rootDir+"\I125A_gr.txt"
sourcedatafile = rootDir+"\I125A_source_data.txt"
 
#create object for calculation
o = jkcm_samemodel_multisource_TG43()
#initialize object with TG43 data
o.initializeTG43tables(frthetafile,grfile,sourcedatafile)
 
#display source information
print(o.jkcm_TG43_calc_obj)
 
#file containing source locations
source_positions_filename = "G:\data\src\TG43\COMS_plaques\COMS_20mm_plaque.txt"
 
#import source positions
o.importSources(source_positions_filename)

#calc at point of interest
dose_cGy_arr_per_Sk_per_effT = o.calc_at_point([0,0,0.32])
dose_Gy_per_Sk = np.sum(dose_cGy_arr_per_Sk_per_effT)*q.calc_eff_time(101)/100
Sk_per_seed = 85/dose_Gy_per_Sk
 
 #return answer
print("Sk per seed = {0} U".format(Sk_per_seed))