# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:36:10 2016

@author: charlesgulian
"""

import os
import matplotlib.pyplot as plt
import numpy as np
#from astropy.io import fits

dir_name = os.getcwd()

n = open(dir_name + '/img_name_compare.txt','r')
img_tag1 = (n.readline()).strip()
img_tag2 = (n.readline()).strip()

    
#filename = dir_name+'/AstroImages/'img_tag+'.fits'
#hdulist = fits.open(filename)  # open a FITS file
#fits_hdr = hdulist[0].header
#hdulist.close()

f1 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag1+'_comparison.cat','r')
lines = f1.readlines()
cat_hdr_size = len(lines[len(lines)-1].split()) # Catalog header size; number of output columns in default.param
f1.seek(0) # Reset pointer to start of file
lines = None

params = []
PRINT = False
if PRINT:
    for i in range(cat_hdr_size):
        param = (f1.readline()).split()[2] # Get parameter name
        print param
        params.append(param)
else:
    for i in range(cat_hdr_size):
        param = (f1.readline()).split()[2] # Get parameter name
        params.append(param)

data1 = []
for line in f1:
    line = line.strip()
    column = line.split()
    source = {}
    for i in range(len(params)):
        source[params[i]] = float(column[i])
    source[params[0]] = int(column[0]) # First column = object #
    data1.append(source)
    
f2 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag2+'_comparison.cat','r')
data2 = []
for i in range(cat_hdr_size):
    f2.readline()
for line in f2:
    line = line.strip()
    column = line.split()
    source = {}
    for i in range(len(params)):
        source[params[i]] = float(column[i])
    source[params[0]] = int(column[0]) # First column = object #
    data2.append(source)

if len(data1) != len(data2):
    print "Warning: data arrays are of different lengths"

flux1,flux2 = [],[]
x,y = [],[]
for i in range(len(data1)):
    Obj1,Obj2 = data1[i],data2[i]
    if abs(Obj1['FLUX_BEST']/Obj2['FLUX_BEST']) > 20.0:
        next # Ignore outliers
    flux1.append(Obj1['FLUX_BEST'])
    flux2.append(Obj2['FLUX_BEST'])
    x.append(Obj1['X_IMAGE'])
    y.append(Obj1['Y_IMAGE'])

flux1,flux2 = np.array(flux1),np.array(flux2)
flux_ratio = np.divide(flux1,flux2)
x,y = np.array(x),np.array(y)

# Plotting flux ratios in 3D:

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, flux_ratio, c='b', marker='o')

ax.set_xlabel('X_IMAGE')
ax.set_ylabel('Y_IMAGE')
ax.set_zlabel('Flux Ratio')

plt.show()