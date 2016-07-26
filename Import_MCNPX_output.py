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
import math
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
    
def calc_centers_from_bounds(arr):
    return (0.5*arr[0:-1] + 0.5*arr[1:])
    
def add_in_quadrature(fileList):
    """ This adds the results and uncertainty in quadrature assuming equal weights.
    fileList: a list of filenames containing mdata to add/average together"""
    w = 1/len(fileList)    
    o = Import_MCNPX_output(fileList[0])
    doseArr = w*o['tally_xyz']
    uncArr2 = np.power(w*o['tally_xyz']*o['unc_xyz'],2)
    for i in np.arange(1,len(fileList)):
        o=Import_MCNPX_output(fileList[i])
        doseArr = doseArr+w*o['tally_xyz']
        uncArr2 = uncArr2+np.power(w*o['tally_xyz']*o['unc_xyz'],2)
    
    return({'tally_xyz':doseArr,
            'unc_xyz':np.sqrt(uncArr2)/doseArr,
            'xb':o['xb'],
            'yb':o['yb'],
            'zb':o['zb'],
            'nps':o['nps']})
    

def geometry_r_theta(r,theta,L=0.35):
    """This calculates the geometry function according to TG43. 
    r: radius in cm
    theta: angle in degrees
    L: active source length in cm
    
    @Reference:
        Perez-Calatayud et al "Dose Calculation for Photon-Emitting Brachytherapy Sources with Average
        Energy Higher than 50 keV: Full Report of the AAPM and ESTRO"
        Report of the High Energy Brachytherapy Source Dosimetry (HEBD) Working Group
        2012
        
    """
    if (theta == 0):
        return(1/(r*r-L/2*L/2))
    else:    
        thetaR = math.radians(theta)
        beta1 = np.arccos( (r*math.cos(thetaR) - L/2)/math.sqrt(r*r+L/2*L/2-L*r*math.cos(thetaR)))
        beta2 = np.arccos( (r*math.cos(thetaR) + L/2)/math.sqrt(r*r+L/2*L/2+L*r*math.cos(thetaR)))
        return((beta1-beta2)/(L*r*math.sin(thetaR)))

def F_r_theta(smesh_tally, r, theta, L=0.35, dr=0.001, dtheta=0.001):
    """ calculate F(r,theta) given a smesh tally from MC.
    smesh_tally: the tally_xyz output from Import_MCNPX_output
    r: the radius in cm, where you want to evaluate the F(r,theta)
    theta: the angle in degrees where you want to evaluate F(r,theta)
    L: length (cm) of the active source
    dr: tolerance for checking for radius equivalence in cm 
        This is used to identify index in smesh_tally
    dtheta: tolerance for checking for theta equivalence in cm
        This is used to identify index in smesh_tally
    """
    rc = calc_centers_from_bounds(smesh_tally['xb'])
    thetac = calc_centers_from_bounds(smesh_tally['yb'])
    doseArr = smesh_tally['tally_xyz']
    
    #find r index
    rindex = (rc < (r+dr)) & (rc > (r-dr))    
    thetaindex = (thetac < (theta+dtheta)) & (thetac > (theta-dtheta))    
    theta90index = (thetac < (90+dtheta)) & (thetac > (90-dtheta))
    assert sum(rindex) != 0, "did not find rindex, check your r value and smesh_tally..."    
    assert sum(rindex) == 1, "found more than 1 rindex, try decreasing dr..."
    assert sum(thetaindex) != 0, "did not find thetaindex, check your theta value and smesh_tally..."    
    assert sum(thetaindex) == 1, "found more than 1 thetaindex, try decreasing dtheta..."
    assert sum(theta90index) != 0, "did not find 90index, check your theta value and smesh_tally..."    
    assert sum(theta90index) == 1, "found more than 1 90index, try decreasing dtheta..."
    
    print("Asked for F({0},{1})".format(r,theta))    
    print("Found F({0},{1})".format(rc[rindex],thetac[thetaindex]))
    D_rtheta = doseArr[rindex,thetaindex,0]
    D_r90 = doseArr[rindex, theta90index, 0]
    G_r90 = geometry_r_theta(rc[rindex], thetac[theta90index], L=L)
    G_rtheta = geometry_r_theta(rc[rindex], thetac[thetaindex], L=L)
    
    return(D_rtheta/D_r90*G_r90/G_rtheta)
    

def list_F_r_theta(smesh_tally, rlist, thetalist, drlist, dthetalist):
    """ This just prints out a table of F(r,theta) defined by the list for you"""
    res = np.ndarray([len(rlist), len(thetalist)], dtype=np.double)    
    for i in np.arange(len(rlist)):
        for j in np.arange(len(thetalist)):
            res[i,j] = F_r_theta(smesh_tally, rlist[i], thetalist[j], dr=drlist[i], dtheta=dthetalist[j])
            
    return(res)
    
    