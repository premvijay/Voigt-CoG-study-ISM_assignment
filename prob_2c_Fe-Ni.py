# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:27:48 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as unit
import astropy.constants as const


from analyse_spectrum import lines, spectrum
from voigt_and_W import eq_width

lines_Fe = lines['FeII']

