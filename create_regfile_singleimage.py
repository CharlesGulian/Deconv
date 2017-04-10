# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 17:20:02 2017

@author: charlesgulian
"""

#!/usr/bin/env python
"""
Created on Tue Jun 14 14:06:37 2016

@author: charlesgulian
"""
import os
import numpy as np
import sex_stats
#import matplotlib.pyplot as plt
#import numpy as np
#from astropy.io import fits

curr_dir = os.getcwd()
def create_regfile(catalog_file):
    f = open(catalog_file,'r')
    lines = f.readlines()
    cat_hdr_size = len(lines[len(lines)-1].split()) # Catalog header size; number of output columns in default.param
    f.seek(0) # Reset pointer to start of file
    lines = None
    
    params = []
    PRINT = True
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
        source['NUMBER'] = int(source['NUMBER']) # First column = object #
        data.append(source)
        
    g = open(os.path.join(curr_dir,'AstroImages','Matthias_Test','regfile.reg'),'w')
    # Initializing .reg file
    g.write('# Region file format: DS9 version 4.1\n')
    g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    g.write('image\n')
    for i in range(len(data)):
        Obj = data[i]
        Line = 'ellipse('+str(Obj['X_IMAGE'])+','+str(Obj['Y_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['A_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['B_IMAGE'])+','+str(Obj['THETA_IMAGE'])+')\n'
        g.write(Line)
    g.close()
    return params

create_regfile(os.path.join(curr_dir,'catalog.cat'))