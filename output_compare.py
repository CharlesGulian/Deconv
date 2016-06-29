# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:36:10 2016

@author: charlesgulian
"""

import os
import matplotlib.pyplot as plt
import numpy as np
#from astropy.io import fits

# Later add command line functionality

dir_name = os.getcwd()

n = open(dir_name + '/img_name_compare.txt','r')
img_tag1 = (n.readline()).strip()
img_tag2 = (n.readline()).strip()

    
#filename = dir_name+'/AstroImages/'img_tag+'.fits'
#hdulist = fits.open(filename)  # open a FITS file
#fits_hdr = hdulist[0].header
#hdulist.close()

f1 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag1+'_compare.cat','r')
lines = f1.readlines()
cat_hdr_size = len(lines[len(lines)-1].split()) # Catalog header size; number of output columns in default.param
f1.seek(0) # Reset pointer to start of file
lines = None

params = []
PRINT = True
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
    data1.append(source)
    
f2 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag2+'_compare.cat','r')
data2 = []
for i in range(cat_hdr_size):
    f2.readline()
for line in f2:
    line = line.strip()
    column = line.split()
    source = {}
    for i in range(len(params)):
        source[params[i]] = float(column[i])
    data2.append(source)

if len(data1) != len(data2):
    print "Warning: data arrays are of different lengths"

for i,j in data1.items():
    