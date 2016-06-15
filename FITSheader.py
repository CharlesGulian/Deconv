# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 16:55:20 2016

@author: charlesgulian
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

dir_name = '/Users/annepstein/Documents/Astro Images/'
filename = dir_name+'fpC-007202-r2-0199.fits'
hdulist = fits.open(filename)  # open a FITS file
hdr = hdulist[0].header
hdulist.close()

exptime = float(hdr['exptime'])
flux20 = float(hdr['flux20'])
const = 10**8
# f0 = flux20*exptime*const
f0 = flux20*const # Base flux value (f0)

MAG_ZEROPOINT = 2.5*np.log10(f0)


f = open(dir_name+'fpC-007202-r2-0199.cat','r')
cat_hdr_size = 11

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
