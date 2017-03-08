# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 21:33:10 2017

@author: charlesgulian
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits

# Cutout new co-adds from stitched co-adds 

# /home/DATA/STRIPE82_330-360_Stitched/

curr_dir = os.getcwd()
img_dir =  os.path.join(curr_dir,'AstroImages','Stitched_Coadds')

stitched_image_file =  os.path.join(img_dir,'fpC-004874-2.0-221.0-287.0_stitched.fits')
stitched_image = fits.getdata(stitched_image_file)

region_center_dict = {}
# Sparse regions
region_center_dict[1] = [1024,65466]
region_center_dict[2] = [1072,77946]
region_center_dict[3] = [1054,53218]
# Crowded regions
region_center_dict[4] = [1046,47134]

def cut(x,y,image):
    # Cuts out a 1600X1600 region centered at (x,y)    
    cutout_image = image[x-800:x+800,y-800:y+800]
    return cutout_image
    
img_list = []
for i in range(4):
    # Initially swap X,Y from DS9 image coordinates
    y,x = region_center_dict[i+1]
    img = cut(x,y,stitched_image)
    img_list.append(img)
    # Swap back X,Y to match DS9 image coordinates
    x,y = y,x

    # Create new .fits image for cutout region
    cutout_image_file = os.path.join(curr_dir,'AstroImages','Coadd','fpC-4874-x{0}-y{1}_cutout.fits'.format(x,y)) 
    
    header = fits.Header()
    header['OG_IMG'] = os.path.split(stitched_image_file)[1]
    header['CUT_X'] = x
    header['CUT_Y'] = y    
    
    #its.writeto(cutout_image_file,img,header=header,clobber=True)
    
