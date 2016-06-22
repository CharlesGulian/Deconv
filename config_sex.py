#!/usr/bin/env python

"""
Created on Fri Jun 17 16:28:41 2016

@author: charlesgulian
"""

# This script is passed arguments from the command line of a .fits image, .sex config file, .param output parameter file

import os
dir_name = os.getcwd()
import argparse
from astropy.io import fits
import numpy as np
'''
parser = argparse.ArgumentParser()
parser.add_argument('image_file',help='Input .FITS image file name',type=str)
parser.add_argument('config_file',help='Input .sex configuration file name',type=str)
parser.add_argument('param_file',help='Input .param output catalog parameter file name',type=str)

args = parser.parse_args() # Get arguments

c = open(dir_name+'/'+args.config_file,'r')
p = open(dir_name+'/'+args.param_file,'r')
'''
Q = open()

# Number of lines/params in config_file and param_file constant
# We can write a dictionary from each file:

config_dict = {}
param_dict = {}

# Cleaning .sex config_file and writing dictionary
for line in c:
    line = line.strip()
    #print line.split()
    if (line.split() != [] and ((line.split())[0])[0] != '#'):
        config_dict[(line.split())[0]] = (line.split())[1]
        #print line.split()[0],line.split()[1]
    else:
            next
c.close()

for line in p:
    line = (line.strip())
    param_dict[((line.split())[0])[1::]] = False # True/False to determine output status

p.close()

# Is this useful? Not sure, but it does provide a fast and convenient way of searching for a certain
# parameter in config_file or param_file and changing its value or status

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Get image filename "tag" (without .fits extension)
img_tag = args.image_file[0:len(args.image_file)-5]

# MAG_ZEROPOINT computed from .FITS header keywords:
hdulist = fits.open(dir_name+'/AstroImages/'+args.image_file)
fits_hdr = hdulist[0].header # Get header from .FITS image
hdulist.close()

flux20 = float(fits_hdr['flux20'])
const = 10**8
f0 = flux20*const # Base flux value (f0)
config_dict['MAG_ZEROPOINT'] = 2.5*np.log10(f0)

# Specifying other parameters:
config_dict['CATALOG_NAME'] = dir_name+'/Results/'+img_tag+'.cat'
config_dict['CHECKIMAGE_NAME'] = dir_name+'/Results/'+img_tag+'_checkimg.fits'

#print len(config_dict), len(param_dict)
c_new = open(dir_name+'/'+args.config_file,'r+')
for key,val in config_dict.items():
    c_new.write(key+' '+val)
c.seek(0)
print c.readlines()


