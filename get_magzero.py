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

plot = False
if plot: 
    f = open(dir_name+'fpC-007202-r2-0199.cat','r')
    cat_hdr_size = 12
    
    PRINT = True
    for i in range(cat_hdr_size):
        if PRINT:
            print f.readline()
        else:
            f.readline()
    
    data = []
    magnitudes = []
    for line in f:
        line = line.strip()
        column = line.split()
        magnitude = float(column[4])
        source = {}
        source['NUMBER'] = int(column[0])
        source['MAG_BEST'] = float(column[4])
        magnitudes.append(magnitude + MAG_ZEROPOINT)
        data.append(source)
        
    data, magnitudes = np.array(data), np.array(magnitudes)
    n, bins, other = plt.hist(magnitudes, 35, facecolor='green', alpha=0.75)
    plt.title('MAG_BEST + MAG_ZEROPOINT: MAG_ZEROPOINT = ' + str(MAG_ZEROPOINT)[0:7])
    plt.xlabel('MAG_BEST')
    plt.ylabel('Frequency (N)')
    #plt.axis([])
    #plt.grid(True)
    plt.show()