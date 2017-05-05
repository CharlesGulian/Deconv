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
deconvolvedImage4 = Image('AstroImages/Deconvolved/transposed_dir_deconv.fits','Deconvolved','Transposed_Dir')
deconvolvedImage5 = Image('AstroImages/Deconvolved/twice_transposed_deconv.fits','Deconvolved','Twice_Transposed')
deconvolvedImage6 = Image('AstroImages/Deconvolved/shuffled_deconv.fits','Deconvolved','Shuffled')
deconvolvedImage7 = Image('AstroImages/Deconvolved/shuffled_seed1_deconv.fits','Deconvolved','Shuffled_Seed1')
deconvolvedImage8 = Image('AstroImages/Deconvolved/initialized_complete_deconv.fits','Deconvolved','Initialized_Complete')
deconvolvedImage9 = Image('AstroImages/Deconvolved/initialized_rm0-9_deconv.fits','Deconvolved','Initialized_rm0-9')
deconvolvedImage10 = Image('AstroImages/Deconvolved/initialized_rm10-19_deconv.fits','Deconvolved','Initialized_rm10-19')
deconvolvedImage11 = Image('AstroImages/Deconvolved/initialized_rm20-29_deconv.fits','Deconvolved','Initialized_rm20-29')
deconvolvedImage12 = Image('AstroImages/Deconvolved/og_initialized_complete_deconv.fits','Deconvolved','OG_Initialized_Complete')
deconvolvedImage13 = Image('AstroImages/Deconvolved/og_initialized_smoothed_complete_deconv.fits','Deconvolved','OG_Initialized_Smoothed_Complete')


imgs = [deconvolvedImage13]
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
