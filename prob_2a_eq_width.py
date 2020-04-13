# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:09:31 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from analyse_spectrum import spectrum,lines

fig, ax = plt.subplots(figsize=(12,6))
    
ax.plot(spectrum[:,0],spectrum[:,1]/spectrum[:,3])
ax.plot(spectrum[:,0],spectrum[:,3],'-')
ax.set_ylim(0,1.2)

for ID in lines:
    for line in lines[ID]:
        ax.vlines(line.start,0,1.2)
        ax.vlines(line.stop,0,1.2)
        ax.text(line.lam_0,0.2,line.ID)

#with PdfPages('multipage_pdf.pdf') as pdf:
#    lam_start = 1125
#    while lam_start<1780:
#        lam_end= lam_start + 28
#        ax.set_xlim(lam_start,lam_end)
#        pdf.savefig()
#        lam_start = lam_start+25



#ax.plot(spectrum[:,0],spectrum[:,3])
#ax.set_xlim(1250,1270)
#ax.set_ylim(0,1.1)
    
#for line in lines[:]:
#    print(vars(line))

for ID in lines:
    for line in lines[ID]:
        print("{},{:.3f},{:.4e},{:.3e},{:.4f},{:.4e}".format(line.ID,line.lam_0,line.f,line.gamma,line.W,line.W_by_lam))
        























