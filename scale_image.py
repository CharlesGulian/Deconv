# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 16:43:43 2016

@author: charlesgulian
"""

import os
from astropy.io import fits

def scale_image(image_file):
    
    header = fits.getheader(image_file)
    data = fits.getdata(image_file)
    try:
        flux20 = header['flux20']
        sky = header['sky']
        softbias = header['softbias']
    except KeyError:
        print 'Error: Missing keyword in image header'
    
    alpha = 1e-8/(flux20)
    flux = alpha*(data - softbias - sky)
    return flux

curr_dir = os.getcwd()
image_dict = {}

complete_header_images = []
with open(os.path.join(curr_dir,'complete_header_imagelist.txt'),'r') as g:
    g.readline()
    lines = g.readlines()
    for line in lines:
        print line[0:len(line)-1]
        complete_header_images.append(line[0:len(line)-1])
        
        line_split = line.split('/')
        image_file = line_split[len(line_split)-1]
        image_file = image_file[0:len(image_file)-1]
        complete_header_images.append(image_file)

for image in complete_header_images:
    image_dict[image] = True


    
