# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 21:32:56 2016

@author: justin

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
                        VERBOSE=0):
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
    ;            List the procedures/functions that this method calls here.
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
    ;         Please put an example::
    ;
    ;         ;example IDL code
    ;         ; output success
    ;
    ;     regular comment again
    ;
    ; :History:
    ;         Created by Justin Mikell on Sep 21, 2010
                Ported to python by Justin Mikell on July 21, 2016
                No POS_VOLUME_OBJ support yet.
    """ 
    f = open(filename, 'r')      
    data = f.read()
    f.close()
    mylines = data.splitlines()
    n_lines = len(mylines)
    print("Importing data from :{0}".format(filename))
    if(VERBOSE > 0):
       print("first lines of file preceding data are...")
       for i in np.arange(9):
          print("{0}:{1}".format(i,mylines[i]))
            
    print("Removing extra whitespace..")
    for i in np.arange(len(mylines)):
        mylines[i] = " ".join(mylines[i].split())
    if(VERBOSE > 0):
       print("first lines after removing extra whitespace...")
       for i in np.arange(9):
          print("{0}:'{1}'".format(i,mylines[i]))
  
    #get the number of particles simulated
    npsline = mylines[1].split()
    nps = np.double(npsline[1])
    total_meshes = np.uint(npsline[0])
    if (MESH_NUM > total_meshes):
        print("Warning mesh_num: {0} not found. Setting mesh_num=1".format(MESH_NUM))
        MESH_NUM=1 #1=-based counting for mesh
    
    if(VERBOSE > 0):
        print("nps: {0}".format(nps))
        print("total_meshes: {0}".format(total_meshes))
    
    #read in all the mesh dimensions
    mesh_dimensions_dict = {}
    for i in np.arange(total_meshes):           
        boundsline = mylines[np.int(3+4*i)].split()
        if(VERBOSE > 0):
            print("reading bounds for mesh {0}".format(i))
            print("boundsline:{0}".format(boundsline))
        nx = np.long(boundsline[2])-1
        ny = np.long(boundsline[3])-1
        nz = np.long(boundsline[4])-1
        
        #this will help with 1D and 2D arrays
        if (nx == 0):
            nx = nx+1
        if (ny == 0):
            ny = ny+1
        if (nz == 0):
            nz = nz+1
        
        mesh_dimensions_dict[i] = [nx,ny,nz]
    
    #Each mesh consists of 3+2*ny*nz lines
    # 3 comes from xb,yb,zb
    # 2*nx*ny comes from the dose values and their uncertainties
    first_xb_line = 7+4*(total_meshes-1) #this is past all the "header" info
    xb_line = np.long(first_xb_line)
    print("xb_line:{0}".format(xb_line))
    for i in np.arange(total_meshes-1):
        xb_line = np.long(xb_line + 3 + 2*mesh_dimensions_dict[i][1]*mesh_dimensions_dict[i][2])
    
    print("xb_line:{0}".format(xb_line))
    #read in xbounds, ybounds, zbounds
    xb = np.array(mylines[xb_line].split(), dtype=np.float)
    yb = np.array(mylines[xb_line+1].split(), dtype=np.float)    
    zb = np.array(mylines[xb_line+2].split(), dtype=np.float)
    
    #now set nx,ny,nz based on the mesh we have selected to read
    nx = mesh_dimensions_dict[MESH_NUM-1][0]
    ny = mesh_dimensions_dict[MESH_NUM-1][1]
    nz = mesh_dimensions_dict[MESH_NUM-1][2]
    
    #see if ny,nz agree with rest of file
    if (np.longlong(2)*ny*nz) != (n_lines-10):
        print("This may be a spherical/cylindrical mesh with an implicit 0.")
        if(yb[ny] == 180):
            print("This appears to be a spherical/cylindrical mesh file")
        ans = 'NA'
        while ( ans.upper() != 'Y' and ans.upper() != 'N'):
            print("Treat as spherical/cylindrical with implicit 0 in polar angle? [y/n]")
            ans = input('[y/n]')
        if ans.upper() == 'Y':
            yb = np.append(0, yb)
            ny = ny + 1
    
    #create the funcrz
    print("reading in dose values now...")
    tally_xyz = np.ndarray([nx,ny,nz], dtype=np.double)
    for k in np.arange(nz):
        for j in np.arange(ny):
            temp = np.array(mylines[np.long(ny*k+j+(xb_line+3))].split(), dtype=np.double) #data starts 3 lines after xb
            if(VERBOSE > 0):
                print("Line number: {0}".format(np.long(ny*k+j+(xb_line+4))))
                print("total elements on line: {0}".format(len(temp)))
            
            if(len(temp) != nx):
                print("number of imported elements doesnt match expected!!!")
                exit(1)
            tally_xyz[:,j,k] = temp
    
    #now read the uncertainty    
    print("reading in uncertainty values...")
    unc_xyz = np.ndarray([nx,ny,nz], dtype=np.double)
    for k in np.arange(nz):
        for j in np.arange(ny):
            temp = np.array(mylines[np.long(ny*k+j+ny*nz+(xb_line+3))].split(), dtype=np.double)
            if(VERBOSE > 0):
                print("Line number: {0}".format(np.long(ny*k+j+ny*nz+(xb_line+4))))
                print("total elements on line: {0}".format(len(temp)))
            
            if(len(temp) != nx):
                print("number of imported elements doesnt match expected!!!")
                exit(1)
            unc_xyz[:,j,k] = temp
    
    tally_xyz = tally_xyz*DOSE_FACTOR
    return({'tally_xyz':tally_xyz, 'unc_xyz':unc_xyz, 'xb':xb, 'yb':yb, 'zb':zb, 'nps':nps})
