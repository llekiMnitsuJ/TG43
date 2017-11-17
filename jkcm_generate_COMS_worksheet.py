# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:00:46 2017

@author: jusmikel

This file reads does the following:
    1) reads TG43 source data
    2) reads source positions in a COMS eyeplaque
    3) writes an excel file XLSX file
"""


import numpy as np
import xlsxwriter

#import the libraries
#runfile('G:/data/src/TG43/jkcm_TG43_calc.py', wdir='G:/data/src/TG43')
#runfile('G:/data/src/TG43/jkcm_samemodel_multisource_TG43.py', wdir='G:/data/src/TG43')


def create_top_rows(ws, comsname, titleformat, inputformat):
                             
    ws.write("A1", comsname, titleformat)
    ws.merge_range("A3:B3", "Physicist Name:")
    ws.merge_range("A4:B4", "Calc Date and Time:")
    ws.merge_range("E3:F3", "Patient Name:")
    ws.merge_range("E4:F4", "Patient Reg#:")
             
    ws.merge_range('C3:D3','', inputformat)
    ws.merge_range('C4:D4','', inputformat)
    ws.merge_range('G3:H3','', inputformat)
    ws.merge_range('G4:H4','', inputformat)
        
    return



outfile = r"C:\Users\jusmikel\Desktop\test.xlsx"
wb = xlsxwriter.Workbook(outfile)
#add format for title 
titleformat = wb.add_format()
titleformat.set_bold()
titleformat.set_font('Arial')
titleformat.set_font_size(22)

inputformat = wb.add_format()
inputformat.set_bg_color('#FFBB99')

                         
                         
ws = wb.add_worksheet()
create_top_rows(ws,"10 mm Eyeplaque Hand Calc (IAI-125A)", titleformat, inputformat)
    


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
 
#file containing source locations (problem specified 10x10mm tumor, so add 1 mm margin and use 12x12COMS)
source_positions_filename = rootDir+"\COMS_plaques\COMS_10mm_plaque.txt"
 
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