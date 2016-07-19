# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:34:39 2016

@author: charlesgulian
"""
import os
curr_dir = os.getcwd()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pysex
import sex_stats
import fits_tools
#import sex_config
#import do_config

# Image deconvolution project:
# Main script for data analysis, image comparison, photometric statistics, and more

# Good image comparison
goodImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
goodImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'
goodImage3 = 'AstroImages/Good/fpC-4868-x4211-y138_stitched_alignCropped.fits'
goodImage4 = 'AstroImages/Good/fpC-6383-x5176-y121_stitched_alignCropped.fits'

goodImgs = [goodImage1,goodImage2,goodImage3,goodImage4]

# Bad image comparison:
badImage1 = 'AstroImages/Bad/fpC-5759-x24775-y300_stitched_alignCropped.fits' # Modest gradient from top to bottom
badImage2 = 'AstroImages/Bad/fpC-6548-x24940-y302_stitched_alignCropped.fits' # Modest gradient from top to bottom
badImage3 = 'AstroImages/Bad/fpC-5781-x25627-y293_stitched_alignCropped.fits' # Very weak gradient from bottom left to top right
badImage4 = 'AstroImages/Bad/fpC-7140-x24755-y270_stitched_alignCropped.fits' # Weak gradient from bottom left to top right

badImgs = [badImage1,badImage2,badImage3,badImage4]

for image_file in goodImgs:
    
    mask_file = image_file.replace('.fits','_mask.fits').replace('Good','Masks').replace('Bad','Masks')
    if os.path.exists(mask_file):
        fits_tools.maskImage(image_file,mask_file,image_bias=1000.0)
        
        
        