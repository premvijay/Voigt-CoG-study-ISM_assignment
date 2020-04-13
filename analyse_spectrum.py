# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 20:23:10 2020

@author: premv
"""

import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#from scipy import integrate
import read_atoms


spectrum = np.genfromtxt("file\hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat")
R_lam = 1 - spectrum[:,1]/spectrum[:,3]

        
class Line_spec_data(read_atoms.Line):
    def set_integral_range(self,start,stop):
        self.start = start
        self.stop = stop
    def set_overlap(self,overlap):
        self.overlap = overlap
        
    def find_center_index(self,spectrum):
        i = 0
        while spectrum[i,0] < self.lam_0:
            i += 1
        self.center_index = i
        return i
        
    def find_eq_width(self,spectrum):
        R_lam = 1 - spectrum[:,1]/spectrum[:,3]
        W = 0
        if self.overlap is None:
            i = self.find_center_index(spectrum)
            while R_lam[i]>0:
        #        print(i)
                dW = (R_lam[i] + R_lam[i+1])/2 * (spectrum[i+1,0]-spectrum[i,0])
        #        print(dW)
                W += dW
                i += 1
            stop = spectrum[i,0] 
                
            i = self.find_center_index(spectrum)
            while R_lam[i]>0:
        #        print(i)
                dW = (R_lam[i] + R_lam[i-1])/2 * (spectrum[i,0]-spectrum[i-1,0])
        #        print(dW)
                W += dW
                i -= 1
            start = spectrum[i,0]
            
        elif self.overlap == 'Left':
            i = self.find_center_index(spectrum)
            while R_lam[i]>0:
        #        print(i)
                dW = (R_lam[i] + R_lam[i+1])/2 * (spectrum[i+1,0]-spectrum[i,0])
        #        print(dW)
                W += dW
                i += 1
            stop = spectrum[i,0]
            start = - stop + 2*self.lam_0
            W *= 2
            
        elif self.overlap == 'Right':
            i = self.find_center_index(spectrum)
            while R_lam[i]>0:
        #        print(i)
                dW = (R_lam[i] + R_lam[i-1])/2 * (spectrum[i,0]-spectrum[i-1,0])
        #        print(dW)
                W += dW
                i -= 1
            start = spectrum[i,0]
            stop = - start + 2*self.lam_0
            W *= 2
            
        self.set_integral_range(start,stop)
        self.W = W
        self.W_by_lam = W / self.lam_0
        return W
    
    def find_eq_width_cont_fix(self,spectrum,flux_0):
        W = 0
        R_lam_fix = 1 - spectrum[:,3]/flux_0
        dW = 10
        i = self.find_center_index(spectrum)
        while dW>=0:
    #        print(i)
            dW = (R_lam_fix[i] + R_lam_fix[i+1])/2 * (spectrum[i+1,0]-spectrum[i,0])
    #        print(dW)
            W += dW
            i += 1
        stop = spectrum[i,0]    
        
        dW = 0
        i = self.find_center_index(spectrum)
        while dW>=0:
    #        print(i)
            dW = (R_lam_fix[i] + R_lam_fix[i-1])/2 * (spectrum[i,0]-spectrum[i-1,0])
    #        print(dW)
            W += dW
            i -= 1
        start = spectrum[i,0]
            
        self.set_integral_range(start,stop)
        self.W = W
        self.W_by_lam = W / self.lam_0
        return W
    
    def set_N(self,N):
        self.N = N



lines = read_atoms.read_lines('file/atom_identified.dat')

for ID in lines:
    for line in lines[ID]:
        line.__class__ = Line_spec_data
        line.set_overlap(None)
        
lines['NI'][1].set_overlap('Right')
lines['NI'][2].set_overlap('Left')
lines['SII'][2].set_overlap('Right')
lines['SiII'][2].set_overlap('Left')
lines['CI'][1].set_overlap('Right')


      
def find_index(line,spectrum):
    i = 0
    while spectrum[i,0] < line.lam_0:
        i += 1
    line.set_index(i)
    
for ID in lines:
    for line in lines[ID]:
        line.find_center_index(spectrum)


for ID in lines:
    if ID != 'HI':
        for line in lines[ID]:
            line.find_eq_width(spectrum)
#            print(vars(line))
    else:
        lines['HI'][0].find_eq_width_cont_fix(spectrum,1.6e-14)

#def find_eq_width_manual(spectrum,integ_range):
#    spectrum_cut = spectrum[integ_range[0]:integ_range[1]]
#    return integrate.trapz(1-spectrum_cut[:,1]/spectrum_cut[:,3],spectrum_cut[:,0])

#W1 = find_eq_width_manual(spectrum,)

#fig, axes = plt.subplots(2)
#
#axes[0].plot(spectrum[:,0],spectrum[:,1]/spectrum[:,3])
#axes[0].set_xlim(1140,1150)
#axes[0].set_ylim(0,1.1)
#            
#
#for k in [1,2,3]:
#    axes[0].vlines(lines[k].start,0,1)
#    axes[0].vlines(lines[k].stop,0,1)
#
#axes[1].plot(spectrum[:,0],spectrum[:,1]/spectrum[:,3])
#axes[1].set_xlim(1600,1615)
#axes[1].set_ylim(0,1.1)
#            
#
#for k in [4,5]:
#    axes[1].vlines(lines[k].start,0,1)
#    axes[1].vlines(lines[k].stop,0,1)
   
    
#if __name__=='__main__':   
#    fig, ax = plt.subplots(figsize=(12,6))
#    
#    ax.plot(spectrum[:,0],spectrum[:,1]/spectrum[:,3])
#    ax.plot(spectrum[:,0],spectrum[:,3],'-')
#    ax.set_ylim(0,1.2)
#    
#    for ID in lines:
#        for line in lines[ID]:
#            ax.vlines(line.start,0,1.2)
#            ax.vlines(line.stop,0,1.2)
#            ax.text(line.lam_0,0.2,line.ID)
#    
#    with PdfPages('multipage_pdf.pdf') as pdf:
#        lam_start = 1125
#        while lam_start<1780:
#            lam_end= lam_start + 28
#            ax.set_xlim(lam_start,lam_end)
#            pdf.savefig()
#            lam_start = lam_start+25
#    
#
#
#    #ax.plot(spectrum[:,0],spectrum[:,3])
#    #ax.set_xlim(1250,1270)
#    #ax.set_ylim(0,1.1)
#        
#    #for line in lines[:]:
#    #    print(vars(line))
#    
#    for ID in lines:
#        for line in lines[ID]:
#            print("{},{:.3f},{:.4e},{:.3e},{:.4f}".format(line.ID,line.lam_0,line.f,line.gamma,line.W))
        






















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











