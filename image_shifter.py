# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 01:17:32 2016

@author: charlesgulian
"""

import os
import matplotlib.pyplot as plt
import numpy as np

d1,d2 = 8,8
Image = np.random.rand(d1,d2)
Blank = np.zeros([d1,d2])
Gauss_star = Blank
Gauss_star[4:7,4:7] = 1.0;
Gauss_star[5,5] = 2.0;
Image += Gauss_star
New = Blank
print Image
print Blank

dx,dy = -1,-1
'''
offset = [dx,dy]
for i in range(d1):
    for j in range(d2):
        if (i-dx < d1) and (i-dx >= 0) and (j+dy < d2) and (j+dy >= 0):
            New[i,j] = Image[i-dx,j+dy] # Why does this work?
            
print New
#'''

def shift_image(image,x_offset,y_offset):
    
    dims = np.shape(image) # Image dimensions
    dim1,dim2 = dims[0],dims[1]
    
    blank = np.zeros(dims) + 1e-8 # Define blank array to receive new image data
    shifted_image = blank
    
    dx,dy = x_offset,y_offset
    for i in range(dim1):
        for j in range(dim2):
            if (i+dx < dim1) and (i+dx >= 0) and (j+dy < dim2) and (j+dy >= 0):
                shifted_image[i,j] = image[i+dx,j+dy] # Why does this work?
    
    return shifted_image

New = shift_image(Image,1,1)
#print New


import matplotlib
cmap = matplotlib.cm.Greys_r
plt.subplot(2,1,1)
plt.pcolormesh(Image,cmap=cmap)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolormesh(New,cmap=cmap)
plt.colorbar()
plt.show()