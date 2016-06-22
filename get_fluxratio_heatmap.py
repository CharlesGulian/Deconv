# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:17:33 2016

@author: charlesgulian
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

dir_name = os.getcwd()

#test_filename1 = dir_name+'/AstroImages/fpC-6484-x4078-y134_stitched_alignCropped.fits'
#test_filename2 = dir_name+'/AstroImages/fpC-6584-x4447-y144_stitched_alignCropped.fits'

test_filename1 = dir_name+'/AstroImages/fpC-4868-x4211-y138_stitched_alignCropped.fits'
test_filename2 = dir_name+'/AstroImages/fpC-4874-x4052-y124_stitched_alignCropped.fits'


hdulist1 = fits.open(test_filename1)
hdr1 = hdulist1[0].header
data1 = hdulist1[0].data

hdulist2 = fits.open(test_filename2)
hdr2 = hdulist2[0].header
data2 = hdulist2[0].data

if np.shape(data1) != np.shape(data2):
    print 'Error: data arrays are of different sized'

data1,data2 = data1.astype(float),data2.astype(float)
fluxRatio = np.divide(data1,data2)

fluxRatio_mean = np.mean(fluxRatio)
fluxRatio_std = np.std(fluxRatio)
fluxRatio_meanSubtracted = fluxRatio - fluxRatio_mean

fluxRatio_meanSubtracted_sigmaClipped = fluxRatio_meanSubtracted
for i in range(np.shape(fluxRatio)[0]):
    for j in range(np.shape(fluxRatio)[1]):
        
        if np.abs(fluxRatio_meanSubtracted_sigmaClipped[i,j]) >= 2*fluxRatio_std:
            fluxRatio_meanSubtracted_sigmaClipped[i,j] = 2*fluxRatio_std*np.sign(fluxRatio_meanSubtracted_sigmaClipped[i,j])
        else:
            next

heatmap = plt.pcolormesh(fluxRatio_meanSubtracted_sigmaClipped)
plt.axis([0, np.shape(fluxRatio)[1], 0, np.shape(fluxRatio)[0]])
plt.xlabel('X_IMAGE')
plt.ylabel('Y_IMAGE')
plt.title('Flux Ratio Heat Map; mean subtracted and 2-sigma clipped (fpC-4868-x4211-y138 & fpC-4874-x4052-y124)')
plt.savefig('/Users/annepstein/Work/Deconv/Flux Ratio Heat Map_mean subtracted and 2-sigma clipped (fpC-4868-x4211-y138 & fpC-4874-x4052-y124).png')
plt.show()

# fpC-6484-x4078-y134 & fpC-6584-x4447-y144     Top
# fpC-4868-x4211-y138 & fpC-4874-x4052-y124     Bottom