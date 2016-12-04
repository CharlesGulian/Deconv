# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 10:47:09 2016

@author: charlesgulian
"""

# A program to scale custom co-added images with SDSS co-added image
import numpy as np
from astropy.io import fits

class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
        self.masked = None

# Co-added images
coaddedImage1 = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','_SDSS')
coaddedImage2 = Image('AstroImages/Coadd/custom_coadd_median.fits','Coadded','_Custom_Median')
coaddedImage3 = Image('AstroImages/Coadd/custom_coadd_mean.fits','Coadded','_Custom_Mean')

imgs = [coaddedImage1,coaddedImage2,coaddedImage3]

for img in imgs:
    imageData1 = fits.getdata(coaddedImage1.filename)
    if img.ID == '_SDSS':
        S = 1.0
        imageData2 = imageData1
    else:
        imageData2 = fits.getdata(img.filename)
        S = np.sum(np.multiply(imageData1,imageData2))/np.sum(np.multiply(imageData2,imageData2))
    
    header = fits.getheader(img.filename)
    header['pixel_scale_factor'] = S
    
    fits.writeto(img.filename,imageData2,header=header,clobber=True)