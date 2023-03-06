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
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy import interpolate

class jkcm_TG43_calc:
    """Use this class to 1) import 2D TG43 data and then use it to perform 2D TG43 calculations
    as a second check for a single source. See jkcm_samemodel_multisource_TG43 for TG43 calculations
    with multiple sources of the same model (i.e. same TG43 parameters).
    
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
        self.g_r_filename = None
                
        self.aniso_table_thetas_degree = None
        self.aniso_table_radii_cm = None
        self.aniso_table = None
        self.aniso_interp_table_obj = None


        self.g_r_radii_cm = None
        self.g_r_poly = None #represents gr as  polynomial
        self.g_r_table = None
        self.g_r_interp_table_obj = None
        
        

        #These values will be changed when performing an n-seed calculation
        self.source_center = np.zeros([3])
        self.source_tip = np.zeros([3]) #physical tip of source
        self.source_bottom = np.zeros([3]) #physical bottom of source            
        
        #These are constants, and will likely eventually be moved out of here
        self.half_life_h_dict = {}        
        self.half_life_h_dict["I-125"] = 59.4*24.
        self.half_life_h_dict["Y-90"] = 64.1
        self.half_life_h_dict["Ir-192"] = 73.8*24.
        self.half_life_h_dict["Cs-131"] = 9.7*24.
        self.half_life_h_dict["Pd-103"] = 17*24.
        
    
    def __str__(self):
        def align_field_num(name,value,valueString=False):
            string_width=45
            num_width=15
            if(valueString==False):
                tempStr="{:<"+str(string_width)+"s}{:>"+str(num_width)+"f}"
            else:
                tempStr="{:<"+str(string_width)+"s}{:>"+str(num_width)+"s}"
            return(tempStr.format(name,value)+"\n")
        
        s = align_field_num("seed_length_cm", self.seed_length_cm)
        s = s + align_field_num("eff_source_length_cm", self.eff_source_length_cm)
        s = s + align_field_num("seed_diameter_cm", self.seed_diameter_cm)
        s = s + align_field_num("dose_rate_constant_cGy_per_h_per_U", self.dose_rate_constant_cGy_per_h_per_U)
        s = s + align_field_num("source_name_model", self.source_name_model, valueString=True)
        s = s + align_field_num("radionuclide", self.radionuclide, valueString=True)
        s = s + align_field_num("source_data_filename", self.source_data_filename, valueString=True)
        s = s + align_field_num("aniso_filename", self.aniso_filename, valueString=True)
        s = s + align_field_num("g_r_filename", self.g_r_filename, valueString=True)
        s = s + "{:<45s}{:>10f}{:>10f}{:>10f}".format("source_center_cm",self.source_center[0], self.source_center[1], self.source_center[2])+"\n"
        s = s + "{:<45s}{:>10f}{:>10f}{:>10f}".format("source_tip_cm",self.source_tip[0], self.source_tip[1], self.source_tip[2])+"\n"
        s = s + "{:<45s}{:>10f}{:>10f}{:>10f}".format("source_bottom_cm",self.source_bottom[0], self.source_bottom[1], self.source_bottom[2])+"\n"
        
        return(s)
    
    def print_anisotropy_table(self):
        s = " radius(cm),, "+','.join(["%.3f" % i for i in self.aniso_table_radii_cm]) + '\n'
        ts = "theta(deg),"+"{0:.3f},".format(self.aniso_table_thetas_degree[0])+','.join(["%.3f" % i for i in self.aniso_table[0,:]]) +'\n'
        s += ts
        for i in np.arange(len(self.aniso_table_thetas_degree)-1):
            ts = ",{0:.3f},".format(self.aniso_table_thetas_degree[i+1])+','.join(["%.3f" % i for i in self.aniso_table[i+1,:]]) +'\n'
            s += ts
        print(s)
        
    def calc_eff_time(self, time_in_hours, infinite=False):
        """This returns the effective time in hours of the implant. 
        It does this by performing evaluating the analytic expression representing
        a simple exponential integral"""
        
        half_life_h_dict = self.half_life_h_dict
        
        
        assert self.radionuclide in half_life_h_dict.keys(), "{0} is not listed radionuclides!"
        
        hl_h = half_life_h_dict[self.radionuclide]
        mu = np.log(2)/hl_h
        if(infinite == True):
            eff_t = 1/mu
        else:
            eff_t = -1/mu * (np.exp(-1*mu*time_in_hours) - 1)
        
        return(eff_t)

    def calc_wall_time_from_eff_time(self, eff_time_in_hours):
        """Returns the real or physical wall time in hours required based on 
        the effective time input in hours. """
        half_life_h_dict = self.half_life_h_dict
        assert self.radionuclide in half_life_h_dict.keys(), "{0} is not listed radionuclides!"
        
        hl_h = half_life_h_dict[self.radionuclide]
        mu = np.log(2)/hl_h
        physical_time_h = -1/mu * np.log(1-mu*eff_time_in_hours)
        return(physical_time_h)
        
    
    def eval_g_r_table(self, r):
        if(r < np.min(self.g_r_radii_cm)):
            return(self.g_r_table[0])
        elif( r > np.max(self.g_r_radii_cm)):
            return(self.g_r_table[-1])
        else:
            return(self.g_r_interp_table_obj(r))
            
    def _build_g_r_interp_table(self, kind='linear'):
        assert type(self.g_r_table) != type(None), "please import the gr table first!"
        assert type(self.g_r_radii_cm) != type(None), "please import the gr table first!"
        #I catch the bounds error and correct via nearest neighbor extrapolation
        self.g_r_interp_table_obj = interpolate.interp1d(self.g_r_radii_cm, self.g_r_table, kind=kind, bounds_error=True)
        print("finished building interpolation object for g_r table evaluation!")
        
    def eval_frtheta(self,r,theta):
        return(self.aniso_interp_table_obj(r,theta))

    def _build_frtheta_interp_table(self, kind="linear", bounds_error=False):
        assert type(self.aniso_table) != type(None), "please import the anisotropy table first!"
        assert type(self.aniso_table_radii_cm) != type(None), "please import the anisotropy table first!"
        assert type(self.aniso_table_thetas_degree) != type(None), "please import the anisotropy table first!"
    
        r_arr = self.aniso_table_radii_cm
        th_arr = self.aniso_table_thetas_degree        
        
        rr,tt = np.meshgrid(r_arr,th_arr)
        self.aniso_interp_table_obj = interpolate.interp2d(r_arr, th_arr, self.aniso_table, kind=kind, bounds_error=bounds_error)
        print("finished building interpolation object for anisotropy evaluation\n")
        
    def G_r_theta(self, r, theta, theta_epsilon=0.001):
        """ From Perez-Calatayud et al Medical Physics, Vol. 39, No. 5, May 2012. 
        
        r: radius from source center to calculation point in cm
        theta: angle betweeen source center to source tip and source center to point of calculation in degrees
        effL: effective source length in cm
        """

        effL = self.eff_source_length_cm
        assert effL != None, "please import source data parameters first!"
        
        assert theta >= -1*theta_epsilon, "theta:{0} must be >= 0"
        assert theta <= 180+theta_epsilon, "theta:{0} must be <= 180"
        if (theta <= theta_epsilon):
            print("assuming theta({0}) is 0 degrees\n".format(theta))
            return(1./(r*r - effL*effL/4))
        if (np.abs(theta-180) <= theta_epsilon ):
            print("assuming theta({0}) is 180 degrees\n".format(theta))
            return(1./(r*r - effL*effL/4))

        angle = np.radians(theta)
        
        rcosth = r*np.cos(angle)
        term1 = np.arccos((rcosth - effL/2.)/np.sqrt(r*r + np.power(effL/2.,2) - effL*rcosth))
        term2 = np.arccos((rcosth + effL/2.)/np.sqrt(r*r + np.power(effL/2.,2) + effL*rcosth))
        numerator = term1 - term2
        denominator = effL*r*np.sin(angle)
        result = numerator/denominator
        return(result)
        

    def _calc_to_point(self, pos, verbose=0, paramMap=False):
        """perform a TG43 calculation
        pos: an np array of length 3. pos must be in the same units as self.source_center.
        This is typically cm. 
        
        paramMap: a boolean. When true, will return a map/dict containing all the information
        used in the calculation. 
        
        doserate = (Sk)*(drc)*G(r,theta)/G(1,90)*g(r)*F(r,theta)
        
        This function returns the doserate/(Sk) at pos.
        """
        #determine r
        r = self._length(pos, self.source_center)
        
        #get distance along and away from source
        source_vec = (self.source_tip - self.source_bottom)/self._length(self.source_tip, self.source_bottom)
        point_vec = pos - self.source_center
        along = np.dot(point_vec, source_vec)
        along_vec = along*source_vec
        away_vec = point_vec - along_vec
        away = self._length(point_vec, along_vec)
        
        #determine theta
        theta = np.degrees(np.arccos(along/r))
        if (verbose > 0):
            print("r,theta = {0},{1}".format(r,theta))
        
        grtheta = self.G_r_theta(r,theta)
        gr0theta0 = self.G_r_theta(1,90)
        frtheta = self.aniso_interp_table_obj(r,theta)
        gr = self.eval_g_r_table(r)
        drc = self.dose_rate_constant_cGy_per_h_per_U
        
        result = drc*(grtheta/gr0theta0)*gr*frtheta[0]
        
        
        paramMapOut = {}
        paramMapOut['source_name_model'] = self.source_name_model
        paramMapOut['frtheta_file'] = self.aniso_filename
        paramMapOut['gr_filename'] = self.g_r_filename
        paramMapOut['source_data_filename']=self.source_data_filename
        paramMapOut['calcpt_pos'] = pos
        paramMapOut['src_tip_pos'] = self.source_tip
        paramMapOut['src_center_pos'] = self.source_center
        paramMapOut['src_bottom_pos'] = self.source_bottom
        paramMapOut['src_unit_vec'] = source_vec
        paramMapOut['calcpt_minus_src_center_vec'] = point_vec
        paramMapOut['along_distance'] = along
        paramMapOut['along_vec'] = along_vec
        paramMapOut['away_distance'] = away
        paramMapOut['away_vec'] = away_vec
        paramMapOut['theta'] = theta
        paramMapOut['radius'] = r
        paramMapOut['grtheta'] = grtheta
        paramMapOut['gr0theta0'] = gr0theta0
        paramMapOut['frtheta'] = frtheta[0]
        paramMapOut['gr'] = gr
        paramMapOut['doserateconstant'] = drc
        paramMapOut['D(r,theta)/(U-h)'] = result
        
        
        
        if(verbose > 0):
            print("###########calc point: {0}".format(pos))
            print("###########src center: {0}".format(self.source_center))
            print("###########src tip: {0}".format(self.source_tip))
            print("###########calculated away length: {0}".format(away))
            print("###########calculated along length: {0}".format(along))
            print("###########calculated r: {0}".format(r))
            print("###########calculated theta(deg): {0}".format(theta))
            print("###########calculated G(r,theta): {0}".format(grtheta))
            print("###########calculated G(r0,theta0): {0}".format(gr0theta0))
            print("###########calculated F(r,theta): {0}".format(frtheta))            
            print("###########calculated g(r): {0}".format(gr))
            print("###########calculated drc: {0}".format(drc))
        
        if(paramMap == False):
            return(result)
        else:
            return(paramMapOut)

    def calc_to_points(self, arr, verbose=0):
        """
        This calculates the doserate from the current source position to the 
        points defined in the arr. 
        
        arr: an Nx3 array representing spatial coordinates for each point of interest.
        
        Returns a an array of length n corresponding to the dose rate at each point 
        from the current source position and orientation.
        """
        assert np.size(np.shape(arr)) == 2, "expected Nx3 array"
        assert np.shape(arr)[1] == 3, "expected Nx3 array"
        numPoints = np.shape(arr)[0]
        resultArr = np.zeros(numPoints)
        
        for i in np.arange(numPoints):
            resultArr[i] = self._calc_to_point(arr[i,:],verbose=verbose)
        
        return(resultArr)
        
  
    def generate_handcalc_paramaters_table_for_points(self, arr, Sk_U, time_h, verbose=0, header=True):
        assert np.size(np.shape(arr)) == 2, "expected Nx3 array"
        assert np.shape(arr)[1] == 3, "expected Nx3 array"
        numPoints = np.shape(arr)[0]
        
        s = self._generate_handcalc_parameters_table_for_point(arr[0,:], Sk_U, time_h, verbose=verbose)
        for i in np.arange(numPoints-1):        
            s += "\n" + self._generate_handcalc_parameters_table_for_point(arr[i+1,:], Sk_U, time_h, verbose=verbose, header=False)
        return(s)
        
  
    def _generate_handcalc_parameters_table_for_point(self, pos, Sk_U, time_h, verbose=0, header=True):
            
            q = self._calc_to_point(pos, verbose=verbose, paramMap=True)
            
            h = "source model"
            r = q['source_name_model']
            
            h += ",calcpt X/cm, calcpt Y/cm, calcpt Z/cm"
            x = q['calcpt_pos']
            r += ",{0},{1},{2}".format(x[0],x[1],x[2])
    
            h += ", src_bottom X/cm, src_bottom Y/cm, src_bottom Z/cm"
            x = q['src_bottom_pos']
            r += ",{0},{1},{2}".format(x[0],x[1],x[2])
            
            h += ", src_center X/cm, src_center Y/cm, src_center Z/cm"
            x = q['src_center_pos']
            r += ",{0},{1},{2}".format(x[0],x[1],x[2])
            
            h += ", src_tip X/cm, src_tip Y/cm, src_tip Z/cm"
            x = q['src_tip_pos']
            r += ",{0},{1},{2}".format(x[0],x[1],x[2])
            
            h += ", src_unit_vec X/cm, src_unit_vec Y/cm, src_unit_vec Z/cm"
            x = q['src_unit_vec']
            r += ",{0},{1},{2}".format(x[0],x[1],x[2])
            
            
            h+= ",along_vec_x/cm,along_vec_y/cm,along_vec_z/cm, along_length/cm"
            x = q['along_vec']
            r+= ",{0},{1},{2},{3}".format(x[0],x[1],x[2],q['along_distance'])
            
            h+= ",away_vec_x/cm,away_vec_y/cm,away_vec_z/cm,away_length/cm"
            x = q['away_vec']
            r+= ",{0},{1},{2},{3}".format(x[0],x[1],x[2],q['away_distance'])
            
            h += ", src_center_to_calcpoint radius/cm"
            r += ",{0}".format(q['radius'])
            
            h += ", angle between (src_tip - src_center) to calcpt / degrees"
            r += ",{0}".format(q['theta'])
            
            h += ", G(r;theta), G(r_0;theta_0), g(r), F(r;theta), lambda (cGy/U), Sk(cGy*cm**2/h), dwell time(h), D(r;theta)/(U-h), D(r;theta)"
            r += ",{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(q['grtheta'],q['gr0theta0'],q['gr'],q['frtheta'],q['doserateconstant'],Sk_U,time_h, q['D(r,theta)/(U-h)'], q['D(r,theta)/(U-h)']*Sk_U*time_h)
            
            if(header == True):
                r = h +"\n" + r
            
            return(r)
                        

            
            
            
        
    def _calcCenter(self, arr):
        return(0.5*arr[0:-1]+0.5*arr[1:])

    def _length(self, v1, v2):
        x = v1 - v2
        x = np.power(x,2)
        x = np.sqrt(np.sum(x))
        return(x)
        
    def setSourceCenterPos(self, x0,y0,z0):
        self.source_center[0] = x0
        self.source_center[1] = y0
        self.source_center[2] = z0

    def setSourceCenterAndTipPos(self, x0,y0,z0, xt,yt,zt):
        self.setSourceCenterPos(x0,y0,z0)
        self.source_tip[0] = xt
        self.source_tip[1] = yt
        self.source_tip[2] = zt
        delta = self.source_tip - self.source_center
        self.source_bottom = self.source_center - delta
        self._checkAndSetSourceLengthBasedOnCenterOfEnds()

    def _unitVectorOfSource(self):
        length = self._length(self.source_tip, self.source_bottom)
        unit_vector = (self.source_tip - self.source_bottom)/length
        return(unit_vector)
        
    def _checkAndSetSourceLengthBasedOnCenterOfEnds(self):
        epsilon = 0.001
        length = self._length(self.source_tip, self.source_bottom)
        if ( (length > (self.seed_length_cm + epsilon)) or (length < (self.seed_length_cm + epsilon))):
            print("source length was outside tolerance, updating based on center!!!\n")
            unit_vector = (self.source_tip - self.source_bottom)/length
            self.source_center =  0.5*self.source_tip + 0.5*self.source_bottom
            self.source_tip = self.source_center + 0.5*self.seed_length_cm*unit_vector
            self.source_bottom = self.source_center - 0.5*self.seed_length_cm*unit_vector      
            
        
    
    def setSourceTipAndBottomPos(self, xt,yt,zt, xb,yb,zb):
        self.source_tip[0] = xt
        self.source_tip[1] = yt
        self.source_tip[2] = zt
        self.source_bottom[0] = xb
        self.source_bottom[1] = yb
        self.source_bottom[2] = zb        
        self.source_center = 0.5*self.source_tip + 0.5*self.source_bottom
        self._checkAndSetSourceLengthBasedOnCenterOfEnds()


        
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
    
    def import_gr_table(self, filename):
        """
        This will import a text file referred to by filename.
        The test file format consists of comments indicated by "#".
        
        The data is just two columns of text with headers:

        radius_cm g_r
        0.5 1.08
        ...

        Example:
        o = jkcm_TG43_calc()
        filename = "G:\data\src\TG43\I125A_gr.txt"
        o.import_gr_table(filename)
        """
        print("trying to import gr table from: {0}".format(filename))
    
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        mylines = self.remove_comments(mylines, "#")
     
        #find the required parameter strings in file
        pat = re.compile('^[ ]*radius_(.*) g_r[ ]*$')
        r_index=-1
        r_units=None
        for i in range(len(mylines)):
            q = pat.match(mylines[i])
            if(q != None):
                r_index = i
                r_units = q.group(1)
                break
        

        assert r_index != -1, "did not find radius field"
        
        
        print("r_index:{0}\n".format(r_index))
    
        ta_list = mylines[(r_index+1):]
        
        parse_str = " ".join(ta_list)
        parsed = parse_str.split()
        r_arr = parsed[0::2]
        gr_arr = parsed[1::2]
        r_arr = np.array(r_arr, dtype=np.float)
        gr_arr = np.array(gr_arr, dtype=np.float)
        
        
        
        assert len(r_arr) == len(gr_arr), "inconsistent grtheta table size and radii!"
        assert str.upper(r_units) == "CM", "requiring radii in units of cm for now, current units:{0}".format(r_units)
       
        self.g_r_radii_cm = r_arr
        self.g_r_filename = filename
        self.g_r_table = gr_arr
        
        self._build_g_r_interp_table()
        
        #return({"mylines":mylines, "r_arr":r_arr, "gr_arr":gr_arr,"r_units":r_units})
        
                
    def import_aniso_table(self, filename, extrapolation="nn"):
        """         
        This will import a text file referred to by filename.
        The text file format consists of comments indicated by "#"
        In addition it looks for three fields in the following order:
        
        
        "radius_c:" where x is a unit of length e.g. cm. Currently cm is required.
        
        "theta_x:" where x is a unit of angle e.g. deg. Currently deg is required.
        
        "table:" The table will be listed such that a given row represents a given angle.
        
        Data will be separated by whitespace delimiters. 
        For F(r,theta) tables that only go to 90 degrees, the table is assumed
        to be symmetric about 90 degrees; this is relevant for loose cylindrical seeds.
        For cases where cables or the source intefere with reporting, empty values in tables
        should be listed as -1. 
        

        filename: in, string, the text file containing the anisotropy table
        
        extrapolation: in, string, default = "nn" for nearest neighbor extrapolation. 
        This follows TG43 recommendations for extrapolating. Appendix C in TG43-U1.
         
        Example: 
        o = jkcm_TG43_calc()
        filename = "G:\data\src\TG43\I125A_frtheta.txt"
        o.import_aniso_table(filename)
        """
        
        print("trying to import anisotropy table from: {0}".format(filename))
    
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        mylines = self.remove_comments(mylines, "#")
     
        #find the required parameter strings in file
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

        assert str.upper(r_units) == "CM", "requiring radii in units of cm for now"

        assert str.upper(th_units) == "DEG", "requiring theta in units of degrees"
       
       #now check if theta goes to 90 or 180 (update frtheta and theta accordingly)
        if(max(th_arr) == 90):
            print("input table only goes to 90 degree, reflecting across 90 to generate complete table!")
            temp_arr = np.zeros((2*len(th_arr)-1, len(r_arr)), dtype=np.float)
            temp_arr[0:(len(th_arr)),:] = ta_arr
            temp_arr[(len(th_arr)-1):, :] = ta_arr[::-1,:]
            ta_arr = temp_arr
            
            temp_arr = np.zeros(2*(len(th_arr))-1, dtype=np.float)
            temp_arr[0:len(th_arr)] = th_arr
            temp_arr[(len(th_arr)-1):] = (180 - th_arr)[::-1]
            th_arr = temp_arr
       
        #now replace -1 values with nearest neighbor extrapolation through radius
        for i in np.arange(len(th_arr)):
           myrow = ta_arr[i,:]
           index_na = np.where(myrow == -1)
           index_exist = np.where(myrow != -1)
           exist_val = myrow[index_exist[0][0]]           
           if (index_na[0].size != 0):               
               print("updating anisotropy table with nearest neighbor extrapolation\n")
               print("replacing -1 with {0}\n".format(exist_val))
               for j in index_na:
                   myrow[j] = exist_val
           
            
        self.aniso_table_thetas_degree = th_arr
        self.aniso_table_radii_cm = r_arr
        self.aniso_table = ta_arr
        self.aniso_filename = filename
        
        self._build_frtheta_interp_table()
        
        #return({"mylines":mylines, "r_arr":r_arr, "th_arr":th_arr, "frtheta_arr":ta_arr,
        #        "r_units":r_units, "th_units":th_units})
                
        
    def import_source_data(self, filename):
        """         
        This will import a text file referred to by filename.
        The text file format consists of comments indicated by "#"
        In addition it looks for the following fields in the following order:
        
        
        "seed_length_cm:" the length of the physical seed in cm.
        
        "effective_source_length_cm:" the effective length of the radioactive source component in cm.
        
        "seed_diameter_cm:" The physical diameter of the seed in cm.
        
        "dose_rate_constant_cGy_per_U_per_h:" the dose rate constant of the source in cGy/U/h.
        
        "seed_model_name:" a string representing the seed model e.g ("IsoAid-I125A")
        
        "seed_radionuclide:" a string representing the radionuclide e.g. ("Ir-192")
        
        Data will be contained on the same line following the field identifier. 
        

        filename: in, string, the text file containing the anisotropy table
        
         
        Example: 
        o = jkcm_TG43_calc()
        filename = "G:\data\src\TG43\I125A_source_data.txt"
        o.import_source_data(filename)
        """
        
        print("trying to import source data from: {0}".format(filename))
    
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        mylines = self.remove_comments(mylines, "#")
     
        #find the required parameter strings in file
        seed_length_cm_pat = re.compile('^[ ]*seed_length_cm:(.*)')
        eff_source_length_cm_pat = re.compile('^[ ]*effective_source_length_cm:(.*)')
        seed_diameter_cm_pat = re.compile('^[ ]*seed_diameter_cm:(.*)')
        dose_rate_constant_cGy_per_U_per_h_pat = re.compile('^[ ]*dose_rate_constant_cGy_per_U_per_h:(.*)')
        seed_model_name_pat = re.compile('^[ ]*seed_model_name:(.*)')
        seed_radionuclide_pat = re.compile('^[ ]*seed_radionuclide:(.*)')
        
        sl_index=-1
        esl_index = -1
        sd_index = -1
        drc_index = -1
        smn_index = -1
        sr_index = -1
        
        self.seed_length_cm = None  
        self.eff_source_length_cm = None
        self.seed_diameter_cm = None
        self.dose_rate_constant_cGy_per_h_per_U = None
        self.source_name_model = None
        self.radionuclide = None        
        
        for i in range(len(mylines)):
            q = seed_length_cm_pat.match(mylines[i])
            if(q != None):
                sl_index = i                
                self.seed_length_cm = np.float(q.group(1))
                break
        assert sl_index != -1, "did not find seed_length_cm field"
            
        for i in range(sl_index, len(mylines)):
            q = eff_source_length_cm_pat.match(mylines[i])
            if (q != None):
                esl_index = i
                self.eff_source_length_cm = np.float(q.group(1))
                break
        assert esl_index != -1, "did not find effective source length cm field"
        
        for i in range(esl_index, len(mylines)):
            q = seed_diameter_cm_pat.match(mylines[i])
            if (q != None):
                sd_index = i
                self.seed_diameter_cm = np.float(q.group(1))
                break
        assert sd_index != -1, "did not find seed diameter cm field"
        
        for i in range(sd_index, len(mylines)):
            q = dose_rate_constant_cGy_per_U_per_h_pat.match(mylines[i])
            if (q != None):
                drc_index = i
                self.dose_rate_constant_cGy_per_h_per_U = np.float(q.group(1))
                break
        assert drc_index != -1, "did not find dose_rate_constant field"
        
        for i in range(drc_index, len(mylines)):
            q = seed_model_name_pat.match(mylines[i])
            if (q != None):
                smn_index = i
                self.source_name_model = q.group(1).strip()
                break
        assert smn_index != -1, "did not find seed_model_name field"        
        
        for i in range(smn_index, len(mylines)):
            q = seed_radionuclide_pat.match(mylines[i])
            if (q != None):
                sr_index = i
                self.radionuclide = str.upper(q.group(1)).strip()
                break
        assert sr_index != -1, "did not find seed_radionuclide field"        
        

        self.source_data_filename = filename
        

        
                
        

  
            
            
  
            
            
        
           
           
       
       
            