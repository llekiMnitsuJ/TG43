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

class jkcm_mcnpx_rmesh:
    """This class represents an rmesh object output by mcnpx."""
    def __init__(self):
        self.tally_values = np.zeros([1])
        self.unc_values = np.zeros([1])
        self.xb = np.zeros([2])
        self.yb = np.zeros([2])
        self.zb = np.zeros([2])
        self.nps = np.longlong(0)
    
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
        
    
                
    def griddata_yplane(self, xVec, zVec, yindex=0):
        """return the output from griddata at the y index supplied.
        xVec: positions to plot the data (e.g. np.linspace(-10,10,200))
        yVec: positions to plot the data"""
        xc = self.xc()
        zc = self.zc()
        nx = self.nx()
        ny = self.ny()
        nz = self.nz()        
        e = np.reshape(np.repeat(xc, nz), [nx, nz])
        e1 = np.reshape(np.repeat(zc, nx), [nx, nz], order='F')
        imgPlane = griddata(e.flatten(), e1.flatten(), o.tally_values[:,yindex,:].flatten(), xVec, zVec, interp='linear')
        return(imgPlane)
        
    def import_from_mdata_ascii(self, filename):
        """ Example 
        o = jkcm_mcnpx_rmesh()
        filename = "/home/justin/001m"
        o.import_from_mdata_ascii(filename)
        """
        
        f = open(filename, 'r')        
        data = f.read()
        f.close()
        mylines = data.splitlines()
        #get nps associated with mdata
        self.nps = np.longlong(mylines[0].split()[5])
        
       
        #find the nx,ny,nz line
        p = re.compile('^f .*')
        j=0
        for i in np.arange(len(mylines)):
           if(p.match(mylines[i])):
               j=i
               break
       
        q = mylines[j].split()
        nxyz = np.longlong(q[1])
        nx = np.long(q[3])
        ny = np.long(q[4])
        nz = np.long(q[5])
       
        #allocate your arrays
        self.tally_values= np.zeros([nxyz])
        self.unc_values=np.zeros([nxyz])
        self.xb = np.zeros([nx+1])
        self.yb = np.zeros([ny+1])
        self.zb = np.zeros([nz+1])
       
        #read in xbounds
        t=0
        for i in np.arange((j+1), len(mylines)):
           t= t + len(mylines[i].split())
           if(t == (nx+1)):
               break
        temp=''
        for k in np.arange((j+1),i+1):
            temp += mylines[k]
        self.xb = np.array(temp.split(), dtype=np.double)
        print("min/max xb: {0},{1}".format(min(self.xb),max(self.xb)))
        #read in ybounds
        j=i
        t=0
        for i in np.arange((j+1), len(mylines)):
           t= t + len(mylines[i].split())
           if(t == (ny+1)):
               break
        temp=''
        for k in np.arange((j+1),i+1):
            temp += mylines[k]
        self.yb = np.array(temp.split(), dtype=np.double)
        print("min/max yb: {0},{1}".format(min(self.yb),max(self.yb)))
        #read in zbounds
        j=i
        t=0
        for i in np.arange((j+1), len(mylines)):
           t= t + len(mylines[i].split())
           if(t == (nz+1)):
               break
        temp=''
        for k in np.arange((j+1),i+1):
            temp += mylines[k]
        self.zb = np.array(temp.split(), dtype=np.double)
        print("min/max zb: {0},{1}".format(min(self.zb),max(self.zb)))
        
        #advance to tally values
        p = re.compile('^vals.*')
        j=0
        for i in np.arange(len(mylines)):
           if(p.match(mylines[i])):
               j=i
               break
        #read everything into an array
        tempArr = np.zeros([2*nxyz])
        si=0
        fi=0
        for i in np.arange((j+1),len(mylines)):
            temp = np.asarray(mylines[i].split(), dtype=np.double)
            fi = si + len(temp)
            tempArr[si:fi] = temp
            si = fi
        
        #separate absobred dose and uncertainty
        self.unc_values = tempArr[1::2].reshape([nx,ny,nz],order='F')
        self.tally_values = tempArr[0::2].reshape([nx,ny,nz], order='F')
        #xc = 0.5*self.xb[0:-1] + 0.5*self.xb[1:]
        #e= np.reshape(np.repeat(xc,nz),[nx,nz])
        #e1=np.reshape(np.repeat(zc,nx), [nx,nz], order='F')
            
            
        
           
           
       
       
            