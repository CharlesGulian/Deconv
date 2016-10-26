# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 23:51:07 2016

@author: charlesgulian
"""

import os
import glob
import subprocess

sex_path = '/usr/local/bin/sex'
curr_dir = os.getcwd()

# All image filenames must be absolute paths or relative to current working directory

def _call(args):
    try:
        subprocess.call(args,shell=False)
    except subprocess.CalledProcessError:
        print 'Error: sex_call CalledProcessError'
        pass # Handle errors in the called executable
    except OSError:
        print 'Error: OSError'
        pass # Executable not found

def call_sex(imgName,config_file=None,args_ext=[]):
    # Call SExtractor from Python script
    
    args = [sex_path,imgName]
    
    if not os.path.exists(imgName):
        if os.path.exists(os.path.abspath(imgName)):
            imgName = os.path.abspath(imgName)
        else:
            imgName = os.path.join(curr_dir,imgName)
        args[1] = imgName
        
    if config_file != None:
        if not os.path.exists(config_file):
            config_file = os.path.join(curr_dir,config_file)
        args.extend(['-c',config_file])
        
    args.extend(args_ext)
    _call(args)

def loop(img_dir,config_files):
    # Python version of sextractor_loop.sh
    if not os.path.exists(img_dir):
        img_dir = os.path.join(curr_dir,img_dir)
    print 'Image directory: ',img_dir
        
    # List names of .FITS files in img_dir
    imgNames = glob.glob(os.path.join(img_dir,'*.fits'))
    print len(imgNames),' .FITS files found in this directory'
    
    for i in range(len(imgNames)):
        imgName,config_file = imgNames[i],config_files[i]
        print ('Running SExtractor on current image: {}').format(imgName)
        call_sex(imgName,config_file=config_file)    
  
def compare(image_file1,image_file2,config_file1,config_file2,masked_image_file1=None):
    
    if masked_image_file1 != None:
        print('Detecting from {}, measuring from {}').format(masked_image_file1,image_file1)
        call_sex(masked_image_file1+','+image_file1,config_file=config_file1)
    
        print('Detecting from {}, measuring from {}').format(masked_image_file1,image_file2)
        call_sex(masked_image_file1+','+image_file2,config_file=config_file2)
    else:
        print('Detecting from {}, measuring from {}').format(image_file1,image_file1)
        call_sex(image_file1+','+image_file1,config_file=config_file1)
    
        print('Detecting from {}, measuring from {}').format(image_file1,image_file2)
        call_sex(image_file1+','+image_file2,config_file=config_file2)

#testImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
#testImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'

#compare(testImage1,testImage2)
