#!/usr/bin/env python

"""
Created on Fri Jun 17 16:28:41 2016

@author: charlesgulian
"""

import os
dir_name = os.getcwd()
curr_dir = dir_name
from astropy.io import fits
import numpy as np

class configure:
    
    default_params = ['NUMBER','FLUX_APER','FLUXERR_APER','MAG_APER','MAGERR_APER','KRON_RADIUS',\
'X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE']
    #default_params = []
    
    def __init__(self,image_file,config_file,param_file,dual=False,default_params=default_params):
        
        self.config_file,self.param_file = config_file,param_file
        self.dual = dual
        
        if self.dual:
            print 'Configuring SExtractor for dual image mode'
            self.image_file1,self.image_file2 = (image_file).split(',')
        else:
            self.image_file = image_file
                
        self.config_dict = {}
        self.param_dict = {}
        
        # Opening and cleaning .sex config_file, writing dictionary:
        with open(os.path.join(curr_dir,config_file),'r') as c:
            for line in c:
                line = line.strip()
                #print line.split()
                if (line.split() != [] and ((line.split())[0])[0] != '#'):
                    self.config_dict[(line.split())[0]] = (line.split())[1]
                    #print line.split()[0],line.split()[1]
                else:
                        next
    
        # Opening and cleaning .param param_file, writing dictionary:
        with open(os.path.join(curr_dir,param_file),'r') as p:
            for line in p:
                line = (line.strip())
                if ((line.split())[0])[0] == '#': # If first character of line is '#'
                    self.param_dict[((line.split())[0])[1::]] = False # True/False to determine output status
                else:
                    self.param_dict[((line.split())[0])[0::]] = False # True/False to determine output status
        self.output_params = default_params
    
    def reconfigure(self,config_param,value):
        
        if type(config_param) == str:
            self.config_dict[config_param] = value
        
        elif type(config_param) == list:
            for i in range(len(config_param)):
                self.config_dict[config_param[i]] = value[i]
    
    def write_param_file(self,new_param_file=None):
        
        if new_param_file == None:
            new_param_file = self.param_file
        
        for param in self.output_params:
            self.param_dict[param] = True
          
        with open(os.path.join(curr_dir,new_param_file),'w') as p:
            for i,j in self.param_dict.iteritems():
                if j:
                    p.write(str(i)+'\n')
                else:
                    p.write('#'+str(i)+'\n') # Comment out the parameter but still write it to the new file
        
    def write_config_file(self,new_config_file=None,new_param_file=None):
        
        if new_param_file == None:
            new_param_file = self.param_file
        self.config_dict['PARAMETERS_NAME'] = os.path.join(curr_dir,new_param_file)
        
        temp = new_param_file # To avoid confusion in the next line
        self.write_param_file(new_param_file=temp)
        
        if new_config_file == None:
            new_config_file = self.config_file
        
        with open(os.path.join(curr_dir,new_config_file),'w') as c:
            for i,j in self.config_dict.iteritems():
                c.write(str(i)+' '+str(j)+'\n')
            c.close()
                
    def default_config(self):
        
        if self.dual:
            img_tag1 = os.path.split(self.image_file1)[1]
            img_tag2 = os.path.split(self.image_file2)[1]
            img_tag1 = img_tag1[0:(len(img_tag1)-len('.fits'))]
            img_tag2 = img_tag2[0:(len(img_tag2)-len('.fits'))]
        else:
            img_tag = os.path.split(self.image_file)[1]
            img_tag = img_tag[0:(len(img_tag)-len('.fits'))]
        
        # Extracting .FITS header information for configuration file
        if self.dual:
            hdulist = fits.open(os.path.join(curr_dir,self.image_file2)) # Use .fits header from measurement (second) image
        else:
            hdulist = fits.open(os.path.join(curr_dir,self.image_file))
            
        fits_hdr = hdulist[0].header # Get header from .FITS image
        hdulist.close()
        try:
            flux20 = float(fits_hdr['flux20'])
            const = 10**8
            f0 = flux20*const # Base flux value (f0)
            self.reconfigure('MAG_ZEROPOINT',2.5*np.log10(f0))
        except KeyError:
            print 'sex_config.py: KeyError: .FITS header keyword "flux20" not found'            
            self.reconfigure('MAG_ZEROPOINT',28.05)
        
        self.reconfigure('PARAMETERS_NAME',os.path.join(curr_dir,self.param_file))
        self.reconfigure('CHECKIMAGE_TYPE','APERTURES')
        
        if self.dual:
            self.reconfigure('CATALOG_NAME',os.path.join(curr_dir,'Results',(img_tag1+'_'+img_tag2+'_compare.cat')))
            self.reconfigure('CHECKIMAGE_NAME',os.path.join(curr_dir,'Results','CheckImages',(img_tag1+'_'+img_tag2+'_compare_'+self.config_dict['CHECKIMAGE_TYPE']+'.fits')))
        else:
            self.reconfigure('CATALOG_NAME',os.path.join(curr_dir,'Results',(img_tag+'.cat')))
            self.reconfigure('CHECKIMAGE_NAME',os.path.join(curr_dir,'Results',(img_tag+'_'+self.config_dict['CHECKIMAGE_TYPE']+'.fits')))   
        
        # Finding and applying binary mask file to image
        MASK = True
        if self.dual:
            mask_file = self.image_file1.replace('.fits','_mask.fits').replace('Good','Masks').replace('Bad','Masks')
        else:
            mask_file = self.image_file.replace('.fits','_mask.fits').replace('Good','Masks').replace('Bad','Masks')
        if MASK and os.path.exists(mask_file):      
            self.reconfigure('WEIGHT_IMAGE',mask_file)
            self.reconfigure('WEIGHT_TYPE','MAP_WEIGHT')
            
        # Adjust detection threshold
        self.reconfigure('DETECT_THRESH',2.5)
            
# NEXT STEPS:
# 1) Dual image mode options (Done)
# 2) Default configuration options (Done)
# 3) Rewrite pysex.py to use sex_config 
      

""" # MASTER COMMENT LINE#
# For argparse use:
'''
parser = argparse.ArgumentParser()
parser.add_argument('image_file',help='Input .FITS image file name',type=str)
parser.add_argument('config_file',help='Input .sex configuration file name',type=str)
parser.add_argument('param_file',help='Input .param output catalog parameter file name',type=str)
parser.add_argument('-d','--dual',help='Configure for dual-image (comparison) mode; enter comma-separated image ' +\
                     'filenames (e.g. image1.fits,image2.fits)',action='store_true')

args = parser.parse_args() # Get command line arguments

if args.dual:
    print 'Configuring SExtractor for dual image mode'
    image_file1,image_file2 = (args.image_file).split(',')
else:
    image_file = args.image_file

c = open(dir_name+'/'+args.config_file,'r')
p = open(dir_name+'/'+args.param_file,'r')
'''

# For testing use:
'''
filenames = []
Q = open(dir_name+'/do_config_testfiles.txt','r')
for line in Q:
    filenames.append(line.strip())
Q.close()
c = open(dir_name+'/'+filenames[1],'r')
p = open(dir_name+'/'+filenames[2],'r')
#'''





config_dict = {}
param_dict = {}

# Cleaning .sex config_file and writing dictionary:
for line in c:
    line = line.strip()
    #print line.split()
    if (line.split() != [] and ((line.split())[0])[0] != '#'):
        config_dict[(line.split())[0]] = (line.split())[1]
        #print line.split()[0],line.split()[1]
    else:
            next
c.close()

# Cleaning .param param_file and writing dictionary:
for line in p:
    line = (line.strip())
    if ((line.split())[0])[0] == '#': # If first character of line is '#'
        param_dict[((line.split())[0])[1::]] = False # True/False to determine output status
    else:
        param_dict[((line.split())[0])[0::]] = False # True/False to determine output status
p.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Is this useful? Not sure, but it does provide a fast and convenient way of searching for a certain
# parameter in config_file or param_file and changing its value or status
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specifying output parameters:
output_params = ['NUMBER','FLUX_BEST','FLUXERR_BEST','MAG_BEST','MAGERR_BEST','KRON_RADIUS',\
'X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE']
for param in output_params:
    param_dict[param] = True

# Writing new output catalog parameters file from param_dict:
p = open(dir_name+'/'+args.param_file,'w')
for i,j in param_dict.iteritems():
    if j:
        p.write(i+'\n')
    else:
        p.write('#'+i+'\n') # Comment out the parameter but still write it to the new file
p.close()

# Speciying configuration parameters:

# Getting image filename "tag" (without .fits extension)
# For argparse use:
if args.dual:
    img_tag1 = os.path.split(os.path.abspath(image_file1))[1]
    img_tag2 = os.path.split(os.path.abspath(image_file2))[1]
    print img_tag1
    print img_tag2
    img_tag1 = img_tag1[0:(len(img_tag1)-len('.fits'))]
    img_tag2 = img_tag2[0:(len(img_tag2)-len('.fits'))]
else:
    img_tag = os.path.split(os.path.abspath(image_file))[1]
    img_tag = img_tag[0:(len(img_tag)-len('.fits'))]

#img_tag = (filenames[0])[0:len(filenames[0])-5] # For testing use

# Compute MAG_ZEROPOINT from .FITS header keywords (can possible compute other config paramaters from .FITS header):
# For argparse use:
if args.dual:
    hdulist = fits.open(os.path.join(dir_name,image_file2)) # Use .fits header from measurement (second) image
else:
    hdulist = fits.open(os.path.join(dir_name,image_file))
#hdulist = fits.open(dir_name+'/AstroImages/'+filenames[0]) # For testing use
fits_hdr = hdulist[0].header # Get header from .FITS image
hdulist.close()
try:
    flux20 = float(fits_hdr['flux20'])
    const = 10**8
    f0 = flux20*const # Base flux value (f0)
    config_dict['MAG_ZEROPOINT'] = str(2.5*np.log10(f0))
except KeyError:
    print 'do_config.py: KeyError: .FITS header keyword "flux20" not found'
    config_dict['MAG_ZEROPOINT'] = str(28.05)

# Specifying other configuration parameters:
config_dict['PARAMETERS_NAME'] = dir_name+'/'+args.param_file
config_dict['CHECKIMAGE_TYPE'] = 'APERTURES'
if args.dual: # Dual image mode configuration parameters
    config_dict['CATALOG_NAME'] = dir_name+'/Results/'+img_tag1+'_'+img_tag2+'_compare.cat'
    config_dict['CHECKIMAGE_NAME'] = dir_name+'/Results/'+img_tag2+'_compare_'+config_dict['CHECKIMAGE_TYPE']+'.fits'
    print config_dict['CATALOG_NAME']
else: # Single image mode configuration parameters
    config_dict['CATALOG_NAME'] = dir_name+'/Results/'+img_tag+'.cat'
    config_dict['CHECKIMAGE_NAME'] = dir_name+'/Results/'+img_tag+'_'+config_dict['CHECKIMAGE_TYPE']+'.fits'

# Writing new configuration file from config_dict:
c = open(dir_name+'/'+args.config_file,'w')
for i,j in config_dict.iteritems():
    c.write(i+' '+j+'\n')
c.close()


#"""
