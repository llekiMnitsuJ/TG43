# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 21:32:56 2016

@author: justin
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

Use this class to read in and process rmesh data output by mcnpx.
"""

import numpy as np
import os
import re
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt


def Import_MCNPX_output(filename,
                        DOSE_FACTOR=1.0, 
                        MESH_NUM=1,
                        VERBOSE=0
                        )
""" 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;AUTHOR: Justin Mikell, jmikell@mdanderson.org, justin[Dot]mikell@gmail.com
;NAME: Import_MCNPX_output
;PURPOSE: to import a mesh tally (1D,2D, or 3D) from an MCNPX output file
;         
;INPUTS:        
;                
;CATEGORY: Monte Carlo, MCNPX
;CALLING SEQUENCE: resultStruct = Import_MCNPX_output(filename [,/VERBOSE]) 
;KEYWORDS:    
;           DOSE_FACTOR: if set will multiply the tally data by this factor (default value of 1.0)
;           /POS_VOLUME_OBJ: if set will return a pos volume object
;           /VERBOSE: if set will print out information to console
;           
;RESTRICTIONS:
;EXAMPLE:
;-
;MODIFICATION HISTORY:
;Created on 8-06-2010 by Justin Mikell
;;Added removal of extra white space 8-19-2010 JM
;;Added number of particles simulated to return struct 8-28-2010 JM
;TODO:  
;       -add support for multiple meshes
;       -add support for 1D
;       -add support 2D
;       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; docformat = 'rst'
;+
; :Author: Justin Mikell
;           justin[Dot]mikell(AT]gmail.com
;
; :Description:
;    Describe the procedure.
;
; :Uses:
;			List the procedures/functions that this method calls here.
;
; :Params:
;    filename: in, optional, type=string
;    filename represents an MCNPX mesh tally outputfile with the following
;    format::
;    
;       line 1: comment
;       line 2: number_of_meshes number_of_particles_simulated
;         line 3: ignore particle_type?
;         line 4: 3rd integer is number of x boundaries (nx+1) 
;                 4th integer is number of y boundaries (ny+1)
;                 5th integer is number of z boundaries (nz+1)
;         line 5: ignore energy_range?
;         line 6: ignore
;       (repeat 3-6 for number of meshes)
;       line 8: all the xb's are on this line separated by white space (in cm)
;       line 9: all the yb's are on this line separated by white space (in cm)
;       line 10: all the zb's are on this line separated by white space (in cm) (could also be angles depending on mesh)
;       line 11 to line (ny*nz+10): nx*ny*nz values of output in MeV/(g*photon)
;               each column represents x, each line represents y, ny lines implies change in z and reset y.
;       line (ny*nz+11) to line (2*ny*nz+10): nx*ny*nz values of fractional uncertainty
;
; :Returns:  {Tally_xyz:values,Unc_xyz:values, xb:xb,yb:yb,zb:zb, nps:nps}
;      Tally_xyz is nx*ny*nz elements and represents the data stored at 
;      mesh(i,j,k).
;      Unc_xyz is nx*ny*nz elements and represents the fractional uncertainty 
;      at mesh(i,j,k). (e.g. 0.01 means 1%)
;      xb,yb,zb represent the voxel boundaries. (in cm) 
;      (The voxels may not be evenly spaced or the same size.)
;      nps is the number of particles simulated.
;
; :Categories:
;
; :Keywords:
;    DOSE_FACTOR
;    POS_VOLUME_OBJ
;    MESH_NUM
;    VERBOSE
;
; :Examples:
; 		Please put an example::
;
; 		;example IDL code
; 		; output success
;
; 	regular comment again
;
; :History:
; 		Created by Justin Mikell on Sep 21, 2010
            Ported to python by Justin Mikell on July 21, 2016
            No POS_VOLUME_OBJ support yet.
;"""
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        n_lines = len(mylines)
        print("Importing data from :{0}".format(filename))
        if(VERBOSE > 0):
            print("first lines of file preceding data are...")
            for i in np.arange(9):
                print("{0}".format(mylines[i]))
            
        
                
