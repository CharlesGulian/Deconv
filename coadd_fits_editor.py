# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 21:35:06 2016

@author: charlesgulian
"""

import os
os.chdir('/home/cgulian2/Deconv')
curr_dir = os.getcwd()
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# ===============================================================================
# Getting image file paths for co-add frames

# Image directory:
image_dir = '/home/DATA/STRIPE82_330-360_AlignCropped/test7'
image_files = glob.glob(os.path.join(image_dir,'*alignCropped.fits'))
image_files.remove('/home/DATA/STRIPE82_330-360_AlignCropped/test7/fpC-4927-x4127-y118_stitched_alignCropped.fits')
#print type(image_files)
#print image_files[0:5]

# ===============================================================================
# Getting image file paths for unaligned frames to extract header information

# Unaligned image directory
og_dir = 'home/DATA/STRIPE82_330-360'

# Search for co-add frame counterparts in unaligned image directory
og_files = []
for image_file in image_files:
    tag = image_file[0:8]
    tag = tag.replace('-','-00')
    og_file = glob.glob(os.path.join(og_dir,tag+'*.fit'))[0]
    og_files.append(og_file)

# ===============================================================================
# Iteratively opening headers, grabbing keywords, writing to corresponding co-add frame headers
for i in range(len(og_files)):
    og_file = og_files[i]
    image_file = image_file[i]
    
    og_header = fits.getheader(og_file)
    flux0 = og_header['flux0']
    flux20 = og_header['flux20']
    bias = og_header['softbias']
    skylevel = og_header['sky']
    
    header = fits.getheader(image_file)
    header['flux0'] = flux0
    header['flux20'] = flux20
    header['softbias'] = bias
    header['sky'] = skylevel
    
    image_file_copy = image_file.replace('.fits','_copy.fits')
    fits.writeto(image_file_copy,fits.getdata(image_file),header,clobber=True)
    
    del(header)
    del(og_header)



    