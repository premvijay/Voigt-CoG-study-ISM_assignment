# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 20:23:10 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



spectrum = np.genfromtxt("file\hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat")


class Line:
    def __init__(self,ID,lambda_0, f, gamma):
        self.ID = ID
        self.lambda_0 = lambda_0
        self.f = f
        self.gamma = gamma

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











