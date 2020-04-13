# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:04:16 2020

@author: premv
"""

import numpy as np
from scipy import integrate
from scipy.special import wofz
import matplotlib.pyplot as plt
import astropy.units as unit
import astropy.constants as const

c = float(const.c/(unit.m/unit.s))

#from time import time

def H(a,x):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM a
    and Gaussian component 2 * sigma**2 = one.
    """
    return np.real(wofz(x + 1j*a))


sigma_0_dim = np.pi * const.e.gauss**2 / const.m_e /const.c
sigma_0 = np.float(sigma_0_dim.decompose()/unit.m**2*unit.s)



linear_factor = float(np.pi* const.e.gauss**2/const.m_e / const.c**2 / unit.cm**2 *unit.Angstrom)

def W_by_lam_linear(N_lam_f):
    return linear_factor * (N_lam_f)

def W_by_lam_saturated(N_lam_f,b):
    b = b * 1e3#unit.km / unit.s
    N_lam_f_SI = N_lam_f * 1e-6# unit.cm**-2 * unit.Angstrom
    tau_0 = sigma_0  * N_lam_f_SI / (b)
#    print(tau_0) / np.sqrt(np.pi) 2*np.pi *  * np.sqrt(np.pi) 
    return 2 * b / c * np.sqrt(np.log(tau_0))

def W_by_lam_damped(N_lam_f,line):
    return 2*np.sqrt(N_lam_f * linear_factor * line.gamma * line.lam_0 * 1e-10 /c/np.pi**3)
    




def R(lam1,N,b,line):
    lam_0 = line.lam_0 *1e-10#* unit.Angstrom
    om_0 = 2 * np.pi * c / lam_0
    f = line.f
    gamma = line.gamma #* unit.s**-1
    
    N = N * 1e4 #unit.cm**-2
    b = b * 1e3 #unit.km/unit.s
    
    lam = lam1 * 1e-10 + lam_0
    
    a =  np.float((gamma / ( 4 * om_0) * c / b)) 
    om = 2 * np.pi * c / lam
    u = np.array(((om - om_0)/ om_0 * c / b))
    
    sigma = 2 * np.sqrt(np.pi) * lam * f * sigma_0 / b / (2*np.pi) * H(a,u)
    tau = sigma * N
    return 1-np.exp(-tau)

def eq_width(N,b,line):
    return integrate.quad(R, -np.inf,np.inf,args=(N,b,line))[0]


def get_voigt(N,b,line,start,stop):
    lam = np.linspace(start,stop,4000)
    lam1 = lam - line.lam_0
    
    return lam,R(lam1,N,b,line)

def eq_width_trapz(N,b,line,start,stop):    
    voigt = get_voigt(N,b,line,start,stop)
#    plt.plot(lam,1-R_list)
    return float(integrate.trapz(1-voigt[1], voigt[0]))
  
def plot_CoG(N_list,b_list,line_list,gen_plot=True,save_file=None,check_approx=False,legend_id=False,fig_size=(8.5,7)):
    plt.figure(dpi=200,figsize=fig_size)
    for line in line_list:
        for b in b_list:
            W = np.zeros(len(N_list))
            for i in range(len(N_list)):
                W[i] = (eq_width(N_list[i],b,line))
        #    W = np.array(W) 
            if gen_plot:
                plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="{}-{:.0f}-line with b={} km/s".format(line.ID,line.lam_0,b) if legend_id else "b={} km/s".format(b))
            else:
                plt.plot(N_list, W, label="{}-{:.0f}-line with b={} km/s".format(line.ID,line.lam_0,b) if legend_id else "b={} km/s".format(b))
    #            pass
            if check_approx:
                N_lam_f = N_list * line.lam_0 * line.f
                plt.plot(N_lam_f[5:16], W_by_lam_saturated(N_lam_f[5:16],b),'--',label="saturated b = {} km/s".format(b))
        if check_approx:
            plt.plot(N_lam_f[6:], W_by_lam_damped(N_lam_f[6:],line),'--',label="damped approx")

    if check_approx:
#        N_lam_f = N_list * line.lam_0 * line.f
        plt.plot(N_lam_f[:10], W_by_lam_linear(N_lam_f[:10]),'--',label="linear approx")  

    plt.legend(framealpha=0.2)    
    plt.xscale('log')
    plt.yscale('log')
    if gen_plot:
        plt.xlabel(r'$ N \lambda$ f $(cm^{-2} \AA)$')
        plt.ylabel(r'$ W_{\lambda} / \lambda $')
        plt.title("Curve of Growth")
        if save_file is not None:
            plt.savefig(save_file+"_gen.pdf",bbox_inches='tight')
    else:
        plt.xlabel(r'$ N$ $(cm^{-2})$')
        plt.ylabel(r'$ W_{\lambda} (\AA)$')
        plt.title("Curve of Growth")
        if save_file is not None:
            plt.savefig(save_file+".pdf",bbox_inches='tight')
    #plt.figure(figsize=(13,8))
#    plt.close()
    return plt.gcf(), plt.gca()

#N_list = np.logspace(11,22,23) #* unit.cm**-2
##NH = NH_list[9]
##NH = 1 * 10**12.5 * unit.cm**-2
##b_list = np.logspace(0,2,10) #* unit.km/unit.s
#b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
##b = b_list[1]
#
#from analyse_spectrum import lines
#line = lines['HI'][0]

#for b in b_list:
#    W = []
#    for N in N_list:
#        W.append(eq_width(N,b,line))
#    W = np.array(W) 
##    plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="b={} km/s".format(b))
#    plt.plot(N_list, W, label="b={} km/s".format(b))









