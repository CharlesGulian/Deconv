# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 10:47:09 2016

@author: charlesgulian
"""

# A program to scale custom co-added images with SDSS co-added image
import numpy as np
from astropy.io import fits
import fits_tools

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

# Deconvolved image
deconvolvedImage = Image('AstroImages/Deconvolved/deconv.fits','Deconvolved','Deconvolved')
# Trimming deconvolved image to 1600x1600
deconvolvedImageData = fits_tools.getPixels(deconvolvedImage.filename)
TRIM = True
if (np.shape(deconvolvedImageData) != (1600,1600)) and TRIM:
    print 'Trimming deconvolved image'
    padSize = np.array(list(np.shape(deconvolvedImageData))) - np.array([1600,1600])
    deconvolvedImageData_trimmed = deconvolvedImageData[padSize[0]/2:np.shape(deconvolvedImageData)[0]-padSize[0]/2,
                                                        padSize[1]/2:np.shape(deconvolvedImageData)[1]-padSize[1]/2]
    temp = (deconvolvedImage.filename).replace('.fits','_trimmed.fits')
    fits.writeto(temp,deconvolvedImageData_trimmed,fits.getheader(deconvolvedImage.filename),clobber=True)
    deconvolvedImage = Image(temp,'Deconvolved','')

imgs = [coaddedImage1,coaddedImage2,coaddedImage3,deconvolvedImage]
S_array = []
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
    S_array.append(S)
