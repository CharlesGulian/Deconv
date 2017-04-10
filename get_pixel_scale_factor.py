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
deconvolvedImage1 = Image('AstroImages/Deconvolved/deconv.fits','Deconvolved','Deconvolved')
deconvolvedImage2 = Image('AstroImages/Deconvolved/normal_deconv.fits','Deconvolved','Normal')
deconvolvedImage3 = Image('AstroImages/Deconvolved/transposed_deconv.fits','Deconvolved','Transposed')
deconvolvedImages = [deconvolvedImage1,deconvolvedImage2,deconvolvedImage3]

imgs = [coaddedImage2,coaddedImage3,deconvolvedImage2,deconvolvedImage3]
S_array = []
for img in imgs:
    
    if img.category == 'Deconvolved':
        deconvolvedImageData = fits_tools.getPixels(img.filename)
        TRIM = True
        if (np.shape(deconvolvedImageData) != (1600,1600)) and TRIM:
            print 'Trimming deconvolved image'
            padSize = np.array(list(np.shape(deconvolvedImageData))) - np.array([1600,1600])
            deconvolvedImageData_trimmed = deconvolvedImageData[padSize[0]/2:np.shape(deconvolvedImageData)[0]-padSize[0]/2,
                                                                padSize[1]/2:np.shape(deconvolvedImageData)[1]-padSize[1]/2]
            temp = img.filename
            fits.writeto(temp,deconvolvedImageData_trimmed,fits.getheader(img.filename),clobber=True)
            tempimg = Image(temp,'Deconvolved',img.ID)
            img = tempimg
        
    imageData1 = fits.getdata(coaddedImage1.filename)
    imageData2 = fits.getdata(img.filename)
    S = np.sum(np.multiply(imageData1,imageData2))/np.sum(np.multiply(imageData2,imageData2))
    
    header = fits.getheader(img.filename)
    header['PIXSCALE'] = S
    print S
    
    fits.writeto(img.filename,imageData2,header=header,clobber=True)
    S_array.append(S)
