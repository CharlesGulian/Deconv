#!/usr/bin/env python
"""
Created on Tue Jun 14 14:06:37 2016

@author: charlesgulian
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

dir_name = '/Users/annepstein/Documents/Astro Images/'

#filename = dir_name+'fpC-007202-r2-0199'+'.fits'
#hdulist = fits.open(filename)  # open a FITS file
#fits_hdr = hdulist[0].header
#hdulist.close()

n = open(dir_name + 'img_name.txt','r')
n.readline() # Throw away first line of .txt file
for line in n:
    img_tag = line.strip()
    

    f = open(dir_name+img_tag+'.cat','r')
    cat_hdr_size = 11
    
    PRINT = True
    if PRINT:
        for i in range(cat_hdr_size):
            print f.readline()
    else:
        for i in range(cat_hdr_size):
            f.readline()
    
    data = []
    for line in f:
        line = line.strip()
        column = line.split()
        source = {}
        source['NUMBER'] = int(column[0])
        source['X_IMAGE'] = float(column[6])
        source['Y_IMAGE'] = float(column[7])
        source['A_IMAGE'] = float(column[8])
        source['B_IMAGE'] = float(column[9])
        source['THETA_IMAGE'] = float(column[10])
        data.append(source)
        
    g = open(dir_name+img_tag+'.reg','w')
    # Initializing .reg file
    g.write('# Region file format: DS9 version 4.1\n')
    g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    g.write('image\n')
    for i in range(len(data)):
        Obj = data[i]
        Line = 'ellipse('+str(Obj['X_IMAGE'])+','+str(Obj['Y_IMAGE'])+','+str(Obj['A_IMAGE'])+','+str(Obj['B_IMAGE'])+','+str(Obj['THETA_IMAGE'])+')\n'
        g.write(Line)
    g.close()
    
print os.getcwd()