# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 16:44:44 2016

@author: jusmikel
"""

import csv
import numpy as np
import os
import re
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
from scipy import interpolate

class jkcm_samemodel_multisource_TG43:
    """This class is used when you have multiple seeds or dwell positions 
    of the same model implanted inside a patient. For example eye plaques, prostate implants,
    GYN HDR with multiple dwell positions.
    
    The steps are 
        1) first load the TG43 parameters for your source
        2) load the multiple source positions and dwell times
        3) input positions at which to calculate dose. 
    """
    def __init__(self):
        self.jkcm_TG43_calc_obj = jkcm_TG43_calc()  
        #These dictionaries all use the source id as a key. 
        #The source id is a unique integer.
        self.source_center_dict = {} #each looked up value returns a 3 element array
        self.source_tip_dict = {} #each looked up value returns a 3 element array
        self.source_dwell_time_dict = {} #each looked up value returns a single scalar
        self.source_strength_dict = {} #each looked up value returns a single scalar
        self.dwell_time_units = "h"
        self.source_table_filename = None
    
    def listSources(self):
        """This prints out the sources in order of source ID."""
        keys = self.source_center_dict.keys()
        sort_keys = np.sort(list(keys))
        print("sourceID, xc ,yc , zc, xt, yt, zt, time, Sk\n")
        for i in sort_keys:
            c = self.source_center_dict[i]
            t = self.source_tip_dict[i]
            dwell = self.source_dwell_time_dict[i]
            Sk = self.source_strength_dict[i]
            print("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}\n".format(i,c[0],c[1],c[2],t[0],t[1],t[2],dwell, Sk))
    
    def initializeTG43tables(self, frthetafile, grfile, sourcedatafile):
        self.jkcm_TG43_calc_obj.import_aniso_table(frthetafile)
        self.jkcm_TG43_calc_obj.import_gr_table(grfile)
        self.jkcm_TG43_calc_obj.import_source_data(sourcedatafile)
     
    def calc_at_point(self, pos):
        """This loops through all the sources changing the source position and direction and then
        calculating the dose from each source at the given point pos. It returns an array of length N 
        that corresponds to the dose at pos from the sorted source ID."""
        keys = self.source_center_dict.keys()
        sort_keys = np.sort(list(keys))
        result = np.zeros(len(sort_keys))
        for i in np.arange(len(sort_keys)):
            c = self.source_center_dict[sort_keys[i]]
            t = self.source_tip_dict[sort_keys[i]]
            time = self.source_dwell_time_dict[sort_keys[i]]
            strength = self.source_strength_dict[sort_keys[i]]
            self.jkcm_TG43_calc_obj.setSourceCenterAndTipPos(c[0],c[1],c[2], t[0],t[1],t[2])
            result[i] = self.jkcm_TG43_calc_obj._calc_to_point(pos)*time*strength
        
        return(result)
                
        
    def importSources(self, filename):
        """ This imports sources from a text file that is of the following format:
        Probably not a bad idea to inherit from this class and override this function
        for a specific format...
        
        I currently import a COMS file source format, but throw away the angle and 
        assign unity to source strength and dwell time.
        
        Comments begin with #.
        Each row in the file is a single source.
        I use whitespace delimiters. The last two columns (dwelltime) and (airkermastrength) are optional.
        sourceid sourcecenterx sourcecentery sourcecenterz sourcetipx sourcetipy sourcetipz angle 
        """
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        self.source_table_filename = filename
        mylines = data.splitlines()
        mylines = self.jkcm_TG43_calc_obj.remove_comments(mylines, "#")
        print("header in file: {0}".format(mylines[0]))
        for row in mylines[1:]:
            q = row.split()
            qid = np.int(q[0])
            qxc = q[1]
            qyc = q[2]
            qzc = q[3]
            qxt = q[4]
            qyt = q[5]
            qzt = q[6]
            #ignore the other points in records
            self.source_center_dict[qid] = np.array([qxc,qyc,qzc], dtype=np.float)/10.
            self.source_tip_dict[qid] = np.array([qxt,qyt,qzt], dtype=np.float)/10.
            self.source_dwell_time_dict[qid] = 1.
            self.source_strength_dict[qid] = 1.
    
    def setStrengthsInU(self, Sk):
        """ Sk should be in U. It is applied to all seeds in the collection"""
        for i in self.source_strength_dict.keys():
            self.source_strength_dict[i] = Sk
    
            
        
        
        

        
        