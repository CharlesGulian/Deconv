# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 23:30:18 2016

@author: charlesgulian
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def create_mask(image_file,mask_file):
    
    # Get image data, create blank (zeros) array
    image = fits.getdata(image_file)
    mask = np.zeros(np.shape(image))
    
    # Define indices where image pixel values are greater than 5-sigma threshold
    inds = np.where(image >= np.mean(image) + 5.0*np.std(image))
    # Set mask = 1.0 at these indices, save mask to new file
    mask[inds] = 1.0
    fits.writeto(mask_file,mask,clobber=True)
'''
# Test image:
test_image_file = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
test_image = fits.getdata(test_image_file)

mask_file = 'AstroImages/Masks/test_obj_region_mask.fits'
create_mask(test_image_file, mask_file)

mask = fits.getdata(mask_file)
masked_image = np.multiply(test_image,mask)
plt.subplot(2,1,1)
plt.pcolormesh(test_image)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolormesh(test_image-masked_image)
plt.colorbar()
plt.show()
'''

    
    