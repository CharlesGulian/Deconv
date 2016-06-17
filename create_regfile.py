#!/usr/bin/env python
"""
Created on Tue Jun 14 14:06:37 2016

@author: charlesgulian
"""
import os
#import matplotlib.pyplot as plt
#import numpy as np
#from astropy.io import fits

dir_name = os.getcwd()

#filename = dir_name+img_tag+'.fits'
#hdulist = fits.open(filename)  # open a FITS file
#fits_hdr = hdulist[0].header
#hdulist.close()

n = open(dir_name + '/img_name.txt','r')
n.readline() # Ignore first line of .txt file
for line in n:
    img_tag = line.strip()
    
    f = open(dir_name+'/Results/'+img_tag+'.cat','r')
    lines = f.readlines()
    cat_hdr_size = len(lines[len(lines)-1].split()) # Number of output columns in default.param
    f.seek(0) # Reset pointer to start of file
    lines = None
    
    params = []
    PRINT = False
    if PRINT:
        for i in range(cat_hdr_size):
            param = (f.readline()).split()[2] # Get parameter name
            print param
            params.append(param)
    else:
        for i in range(cat_hdr_size):
            param = (f.readline()).split()[2] # Get parameter name
            params.append(param)
    
    data = []
    for line in f:
        line = line.strip()
        column = line.split()
        source = {}
        for i in range(len(params)):
            source[params[i]] = float(column[i])
        source[params[0]] = int(column[0]) # First column = object #
        data.append(source)
        
    g = open(dir_name+'/Results/'+img_tag+'.reg','w')
    # Initializing .reg file
    g.write('# Region file format: DS9 version 4.1\n')
    g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    g.write('image\n')
    for i in range(len(data)):
        Obj = data[i]
        Line = 'ellipse('+str(Obj['X_IMAGE'])+','+str(Obj['Y_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['A_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['B_IMAGE'])+','+str(Obj['THETA_IMAGE'])+')\n'
        g.write(Line)
    g.close()
    
#print os.getcwd()