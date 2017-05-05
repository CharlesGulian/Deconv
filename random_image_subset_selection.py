# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 01:17:57 2017

@author: charlesgulian
"""

# Random image subset selection:

import os
import glob
import numpy as np

curr_dir = os.getcwd()

file_dir = '/home/DATA/STRIPE82_330-360_AlignCropped/test7'
directory = False
if directory:
    files = glob.glob(os.path.join(file_dir,'*alignCropped.fits'))
else:
    image_list_file = os.path.join(curr_dir,'imagelist.txt')
    files = []
    with open(image_list_file,'r') as g:
        lines = g.readlines()
        for line in lines:
            files.append(line.strip())
            
# Make dictionary
file_dict = {}
for i in range(len(files)):
    file_dict[i+1] = files[i]

# Choose random subset of size N, write to new dict
N = len(files) - 1
subset_dict = {}
indices = np.random.randint(1,len(file_dict),N)
for i in range(N):
    subset_dict[i+1] = file_dict[indices[i]]

with open('image_subset_list.txt')
    
        
