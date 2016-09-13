# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 18:11:01 2016

@author: jusmikel
"""
import datetime
import numpy as np
#import the libraries
#runfile('G:/data/src/TG43/jkcm_TG43_calc.py', wdir='G:/data/src/TG43')
#runfile('G:/data/src/TG43/jkcm_samemodel_multisource_TG43.py', wdir='G:/data/src/TG43')

#location of your TG43 data
rootDir="G:\data\src\TG43"
sourceDir=rootDir+"\sources\I125A_consensus"
frthetafile = sourceDir+"\I125A_frtheta.txt"
grfile = sourceDir+"\I125A_gr.txt"
sourcedatafile = sourceDir+"\I125A_source_data.txt"
 
 
#create object for calculation
o = jkcm_samemodel_multisource_TG43()
#initialize object with TG43 data
o.initializeTG43tables(frthetafile,grfile,sourcedatafile)
 
#display source information
print(o.jkcm_TG43_calc_obj)
 
#file containing source locations
source_positions_filename = rootDir+"\COMS_plaques\COMS_20mm_plaque.txt"
 
#import source positions
o.importSources(source_positions_filename)

#set the source strength for all the seeds
o.setStrengthsInU(2.05)

#zero the source strength for seed number 5 to simulate notched plaque
o.source_strength_dict[5] = 0

#calc at point of interest
dose_cGy_arr_per_effT = o.calc_at_point([0,0,0.32])
q = o.jkcm_TG43_calc_obj
dose_Gy_per_effT = np.sum(dose_cGy_arr_per_effT)/100

physical_time_h = q.calc_wall_time_from_eff_time(85/dose_Gy_per_effT)
 
 #return answer
print("physical time required: {0} h".format(physical_time_h))
implant_datetime = datetime.datetime(2015, 5, 2, 8, 30, 0)
print("implant_datetime: {0}".format(implant_datetime))
removal_datetime = implant_datetime + datetime.timedelta(hours=physical_time_h)
print("removal_datetime: {0}".format(removal_datetime))