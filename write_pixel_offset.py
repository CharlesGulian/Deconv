# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:37:01 2016

@author: charlesgulian
"""

# Add (x,y) pixel offset to .FITS header of an image

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
        self.masked = None
        
#image_file = 'AstroImages/Bad/fpC-5781-x25627-y293_stitched_alignCropped.fits'
#new_image_file = ('AstroImages/TEST_IMAGE.fits')

def write_pixel_offset(x_offset,y_offset,image_file,new_image_file=None):
    # Add (x,y) pixel offset to .FITS header of an image
    header = fits.getheader(image_file)
    header['x_offset'] = x_offset
    header['y_offset'] = y_offset
    if new_image_file == None:
        new_image_file = image_file
    fits.writeto(new_image_file,fits.getdata(image_file),header,clobber=True)

#write_pixel_offset(1,2,image_file,new_image_file=new_image_file)

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


#write_pixel_offset(0,0,deconvolvedImage2.filename)
#write_pixel_offset(0,0,deconvolvedImage3.filename)
write_pixel_offset(0,0,deconvolvedImage13.filename)