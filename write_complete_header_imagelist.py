# -*- coding: utf-8 -*-

import glob
import os

image_files = glob.glob(os.path.join('/home/DATA/charlie/STRIPE82_330-360_AlignCropped/test7/','*.fits'))
#print len(image_files)

with open('complete_header_imagelist.txt','w') as g:
    g.write('#List of images with complete header\n')
    for i in range(len(image_files)):
    	g.write('{}\n'.format(image_files[i]))
