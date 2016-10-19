# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:37:01 2016

@author: charlesgulian
"""

# Add (x,y) pixel offset to .FITS header of an image

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

image_file = 'AstroImages/Bad/fpC-5781-x25627-y293_stitched_alignCropped.fits'
new_image_file = ('AstroImages/TEST_IMAGE.fits')

def write_pixel_offset(x_offset,y_offset,image_file,new_image_file=None):
    # Add (x,y) pixel offset to .FITS header of an image
    header = fits.getheader(image_file)
    header['x_offset'] = x_offset
    header['y_offset'] = y_offset
    if new_image_file == None:
        new_image_file = image_file
    fits.writeto(new_image_file,fits.getdata(image_file),header,clobber=True)

write_pixel_offset(1,2,image_file,new_image_file=new_image_file)


