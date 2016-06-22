#!/usr/bin/env python

"""
Created on Fri Jun  3 16:55:20 2016

@author: charlesgulian
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# ** MUST FIX WITH ARGPARSE or OPTPARSE


img_tag = 'fpC-007202-r2-0199'
dir_name = os.getcwd() + '/AstroImages/'
filename = dir_name+img_tag+'.fits'
hdulist = fits.open(filename)  # open a FITS file
hdr = hdulist[0].header
hdulist.close()

exptime = float(hdr['exptime'])
flux20 = float(hdr['flux20'])
const = 10**8
# f0 = flux20*exptime*const
f0 = flux20*const # Base flux value (f0)

MAG_ZEROPOINT = 2.5*np.log10(f0)
print MAG_ZEROPOINT
