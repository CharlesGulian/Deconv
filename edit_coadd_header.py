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

#print image_files

#print type(image_files)
#print image_files[0:5]

# ===============================================================================
# Getting image file paths for unaligned frames, extracting data from header, rewriting to headers of aligned frames

# Unaligned image directory
og_dir = '/home/DATA/STRIPE82_330-360'

# Search for co-add frame counterparts in unaligned image directory
og_files = []
bad_image_files = []
for image_file in image_files:
    image_tag = os.path.split(image_file)[1]
    tag = image_tag[0:8]
    tag = tag.replace('-','-00')
    #print os.path.join(og_dir,tag+'*.fit')
og_file_list = glob.glob(os.path.join(og_dir,tag+'*.fit'))
    if len(og_file_list) == 0:
        print 'Error: non-aligned image not found'
        continue
    # Select an image from unaligned image directory with appropriate header information
    print 'Selecting corresponding image from unaligned directory'
    for i in range(len(og_file_list)):
        og_file = og_file_list[i]
        og_header = fits.getheader(og_file)

        print '\nChecking image {0} header for appropriate keywords'.format(str(i))
        print og_file

        try:
            flux0 = og_header['flux0']
            flux20 = og_header['flux20']
            bias = og_header['softbias']
            skylevel = og_header['sky']
        except KeyError:
            del(og_header)
            continue
        break
    if i == len(og_file_list)-1:
        print 'WARNING: No suitable image found for {0}'.format(image_file)
        bad_image_files.append(image_file)
        continue

    header = fits.getheader(image_file)
    header['flux0'] = flux0
    header['flux20'] = flux20
    header['softbias'] = bias
    header['sky'] = skylevel

    print 'Header data:'
    print 'Sky level: {0} | Flux 20: {1}'.format(skylevel, flux20)

    new_image_file = image_file.replace('DATA','DATA/charlie')
    fits.writeto(new_image_file,fits.getdata(image_file),header,clobber=True)

    og_files.append(og_file)
    del(og_header)
    del(header)

print len(og_files)
print len(bad_image_files)
print len(image_files)



