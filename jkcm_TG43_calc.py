# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:36:16 2016

@author: J Mikell
        

@
        
@References
    Nath et al, "Dosimetry of interstitial brachytherapy sources: Recommendations
    of the AAPM Radiation Therapy Committee Task Group No. 43," Med. Phys. 22, 209–234 (1995).

    Rivard et al, "Update of AAPM Task Group No. 43 Report: A revised AAPM protocol
    for brachytherapy dose calculations". Med. Phys. 31 3
, March 2004
    
    IROC brachytherapy source registry http://rpc.mdanderson.org/RPC/home.htm
    IROC brachytherapy AAPM TG43U! Consensus Brachytherapy Dosimetry Datasets
    
    Perez-Catalyud et al, Dose Calculation for Photon-Emitting Brachytherapy Sources 
    with Average Energy Higher than 50 keV: Full Report of the AAPM and ESTRO High Energy 
    Brachytherapy Dosimetry Working Group. Report 229. August 2012. (https://www.aapm.org/pubs/reports/)
    
"""

import numpy as np
import os
import re
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt

class jkcm_TG43_calc:
    """Use this class to 1) import 2D TG43 data and then use it to perform 2D TG43 calculations
    as a second check. 
    
    It requires the following:
    1) 2D anisotropy table (should be in its own csv text file)
    2) 1D g(r) table (should be in its own csv text file)
    #The following items will be in the same csv file
    3) effective source length in cm
    4) physical source length in cm
    5) dose rate constant (cGy/h/U)
    6) source radionuclide
    7) source name/model
    8) seed diameter in cm
    #The following will be input by the user
    9) Air kerma strength of source (U)
    10) time that source is in position (if temporary)
    
    
    *It uses geometry function, F(r,theta), and g(r) based on 
    the line source, not point source."""
    
    def __init__(self):
        
        self.seed_length_cm = None  
        self.eff_source_length_cm = None
        self.seed_diameter_cm = None
        self.dose_rate_constant_cGy_per_h_per_U = None
        self.source_name_model = None
        self.radionuclide = None
        
        self.source_data_filename = None
        self.aniso_filename = None
        self.gr_filename = None
                
        self.aniso_table_thetas_degree = np.zeros([1])
        self.aniso_table_radii_cm = np.zeros([1])
        self.aniso_table = np.zeros([1,1])

        self.g_r_radii_cm = np.zeros([1])
        self.g_r_poly = np.zeros([1]) #represents gr as  polynomial
        
        
        

        #These values will be changed when performing an n-seed calculation
        self.source_center = np.zeros([3])
        self.source_tip = np.zeros([3]) #physical tip of source
        self.source_bottom = np.zeros([3]) #physical bottom of source        

    
    def _calcCenter(self, arr):
        return(0.5*arr[0:-1]+0.5*arr[1:])
    def xc(self):        
        return self._calcCenter(self.xb)
    def yc(self):
        return self._calcCenter(self.yb)
    def zc(self):
        return self._calcCenter(self.zb)
    def nx(self):
        return len(self.xb) - 1
    def ny(self):
        return len(self.yb) - 1
    def nz(self):
        return len(self.zb) - 1
    def nxyz(self):
        return self.nx()*self.ny()*self.nz()
        
        
    def remove_comments(self, string_list, comment_char):
        """ This removes all extra white space and
        strings starting with the comment_char from input the string_list"""
        #remove empty lines

        new_list = string_list.copy()        
        print("Removing extra whitespace..")
        for i in np.arange(len(new_list)):
            new_list[i] = " ".join(new_list[i].split())        
        
        pattern = re.compile("^[ ]*$")
        remove_set = set()
        for i in new_list:
            q = pattern.match(i)
            if(q != None):
                remove_set.add(i)
            else:
                q = re.search(comment_char, i)
                if( q != None):
                    if(q.span()[0] == 0):
                        remove_set.add(i)
                    else:
                        i = i[0:q.span()[0]]
        
        for i in remove_set:
            new_list.remove(i)
    
        return(new_list)
        
    
    def import_aniso_table(self, filename):
        """ 
        This will import a text file. 
        The text file format consists of comments indicated by "#"
        In addition it looks for three fields in the following order:
        
        "radius_x:" where x is a unit of length e.g. mm. 
        
        "theta_x:" where x is a unit of angle e.g. deg.
        
        "table:" The table will be listed such that a given row represents a given angle.
        
        Data will be separated by whitespace delimiters. 
        For F(r,theta) tables that only go to 90 degrees, the table is assumed
        to be symmetric about 90 degrees; this is relevant for loose cylindrical seeds.
        For cases where cables or the source intefere with reporting, empty values in tables
        should be listed as -1. 
        
         
        Example: 
        o = jkcm_TG43_calc()
        filename = "G:\data\src\TG43\I125A_frtheta.txt"
        o.import_aniso_table(filename)
        """
        
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        mylines = self.remove_comments(mylines, "#")
     
        #find the required parameters
       
        #find the nx,ny,nz line
        radius_pat = re.compile('^[ ]*radius_(.*):.*')
        theta_pat = re.compile('^[ ]*theta_(.*):.*')
        table_pat = re.compile('^[ ]*table:.*')
        r_index=-1
        r_units=None
        th_index=-1
        th_units=None
        ta_index=-1
        for i in range(len(mylines)):
            q = radius_pat.match(mylines[i])
            if(q != None):
                r_index = i
                r_units = q.group(1)
                
                break
        for i in range(r_index, len(mylines)):
            q = theta_pat.match(mylines[i])
            if (q != None):
                th_index = i
                th_units = q.group(1)
                break
        
        for i in range(th_index, len(mylines)):
            q = table_pat.match(mylines[i])
            if (q != None):
                ta_index = i
                break
        assert r_index != -1, "did not find radius field"
        assert th_index != -1, "did not find theta field"
        assert ta_index != -1, "did not find table field"
        
        print("r_index:{0}\nth_index:{1}\nta_index{2}\n".format(r_index, th_index, ta_index))
    
        r_list = mylines[r_index:th_index]
        th_list = mylines[th_index:ta_index]
        ta_list = mylines[ta_index:]
        
        r_str = " ".join(r_list)
        r_str = r_str.split(sep=":")
        r_str = r_str[1]
        r_arr = np.array(r_str.split(), dtype=np.float)
        
        th_str = " ".join(th_list)
        th_str = th_str.split(sep=":")
        th_str = th_str[1]
        th_arr = np.array(th_str.split(), dtype=np.float)
        
        ta_str = " ".join(ta_list)
        ta_str = ta_str.split(sep=":")
        ta_str = ta_str[1]
        ta_arr = np.array(ta_str.split(), dtype=np.float)
        
        assert len(ta_arr) == len(r_arr)*len(th_arr), "inconsistent frtheta table size and radii and thetas!"
        ta_arr = ta_arr.reshape(len(th_arr), len(r_arr))
        
        return({"mylines":mylines, "r_arr":r_arr, "th_arr":th_arr, "frtheta_arr":ta_arr,
                "r_units":r_units, "th_units":th_units})
                
        

  
            
            
        
           
           
       
       
            