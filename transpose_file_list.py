# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 01:23:22 2017

@author: charlesgulian
"""

from astropy.io import fits
import os
import glob

LIST = True
if LIST:
    image_list_file = 'image_subset_list.txt'
    with open(image_list_file,'r') as f:
        lines = f.readlines()
        image_files = []
        for i in range(len(lines)):
            image_file = lines[i].strip()
            image_files.append(image_file)
else:
    image_dir = '/home/DATA/STRIPE82_330-360_AlignCropped/test7/'
    transposed_dir = image_dir+'transposed/'
    image_files = glob.glob(os.path.join(image_dir,'*alignCropped.fits'))

for image_file in image_files:
    image_data = fits.getdata(image_file)
    image_data_transposed = image_data.T
    image_file_transposed = image_file.replace(image_dir,transposed_dir)
    fits.writeto(image_file_transposed,image_data_transposed,header=fits.getheader(image_file),clobber=False)
    
