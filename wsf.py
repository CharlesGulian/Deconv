# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 12:35:27 2017

@author: charlesgulian
"""
import numpy as np
from astropy.io import fits

class Image:

    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
        self.masked = None

# SDSS Co-add
coaddedImage1 = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','_SDSS')

def write_scaling_factor(img):
    
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