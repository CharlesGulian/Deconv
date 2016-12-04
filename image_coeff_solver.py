# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 00:13:11 2016

@author: charlesgulian
"""


import os
os.chdir('/home/cgulian2/Deconv')
curr_dir = os.getcwd()
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

# ===============================================================================
# Getting image file paths

# Image directory:
image_dir = '/home/DATA/Charlie/STRIPE82_330-360_AlignCropped/test7'
image_files = glob.glob(os.path.join(image_dir,'*alignCropped.fits'))
image_files.remove('/home/DATA/Charlie/STRIPE82_330-360_AlignCropped/test7/fpC-4927-x4127-y118_stitched_alignCropped.fits')
#print type(image_files)
#print image_files[0:5]

# ===============================================================================
# Solving for coefficients a and b (alpha and beta)

for i in range(len(image_files)):

    header = fits.getheader(image_file)
    m0 = 0.0
    m1 = 20.0
    c0 = header['flux0']
    c1 = header['flux20']
    b = header['softbias']
    s = header['sky']
    
    # -2.5*np.log10(a*(c0 - b - s)) + b - m0
    # -2.5*np.log10(a*(c1 - b - s)) + b - m1
    
    def equations(p):
        a, b = p
        return (-2.5*np.log10(a*(c0 - b - s)) + b - m0,-2.5*np.log10(a*(c1 - b - s)) + b - m1)
    
    a, b =  fsolve(equations, (1.0, 0.0))
    
    print equations((a, b)) # Should be 0,0