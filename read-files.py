# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 20:23:10 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



spectrum = np.genfromtxt("file\hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat")
R_lam = 1 - spectrum[:,1]/spectrum[:,3]

class Line:
    def __init__(self,ID,lam_0, f, gamma):
        self.ID = ID
        self.lam_0 = float(lam_0)
        self.f = float(f)
        self.gamma = float(gamma)
    def set_index(self,i):
        self.index = i

lines = []

with open('file/Ly_FeII_NiII.dat','rt') as atomfile:
    for line in atomfile.readlines():
        print(line)
        if line[0] != '#':
            lines.append(Line(*line.rsplit()[0:4]))
#        pass


#with open("file\hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat",'rt') as spectrumfile:
#    for line in spectrumfile.readlines():
#        print(line)

#fig, ax = plt.subplots()
#ax.plot(spectrum[:,0],spectrum[:,1]/spectrum[:,3])
#ax.set_xlim(1250,1270)
#ax.set_ylim(0,1.1)





        
        
        
R = integrate.trapz(1-spectrum[:,1]/spectrum[:,3],spectrum[:,0])   

    
        
def find_index(line,spectrum):
    i = 0
    while spectrum[i,0] < line.lam_0:
        i += 1
    line.set_index(i)
    

for line in lines:
    find_index(line,spectrum)


def eq_width(line):
    W = 0
    dW = 0
    i = line.index
    while dW>=0:
        print(i)
        dW = (R_lam[i] + R_lam[i+1])/2 * (spectrum[i+1:0]-spectrum[i,0])
        W += dW
        i += 1
    return W

W = eq_width(lines[2])












#from glob import glob

#from astropy.wcs import WCS
#from astropy.io import fits

# glob searches through directories similar to the Unix shell
#filelist = glob('file/*.fits')
# sort alphabetically - given the way the filenames are
# this also sorts in time
#filelist.sort()

#sp = fits.open(filelist[0])
#sp.info()

#header = sp[0].header

#wcs = WCS(header)


#from PyAstronomy import pyasl
#wvl, flx = pyasl.read1dFitsSpec(filelist[0])


#from specutils.io import read_fits
#myspec = read_fits.read_fits_spectrum1d(filelist[0])











