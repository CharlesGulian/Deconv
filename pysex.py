# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 23:51:07 2016

@author: charlesgulian
"""

import os
import glob
import subprocess
import sex_config

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

def loop(img_dir):
    # Python version of sextractor_loop.sh
    if not os.path.exists(img_dir):
        img_dir = os.path.join(curr_dir,img_dir)
    print 'Image directory: ',img_dir
        
    # List names of .FITS files in img_dir
    imgNames = glob.glob(os.path.join(img_dir,'*.fits'))
    print len(imgNames),' .FITS files found in this directory'

    for imgName in imgNames:
        
        fig = sex_config.configure(imgName,'default.sex','default.param')
        fig.default_config()
        fig.write_config_file(new_config_file='copy.sex',new_param_file='copy.param')
        
        print ('Running SExtractor on current image: {}').format(imgName)
        call_sex(imgName,config_file='copy.sex')    
  
def compare(imgName1,imgName2):
    
    fig = sex_config.configure(imgName1+','+imgName1,'default.sex','default.param',dual=True)
    fig.default_config()
    fig.write_config_file(new_config_file='copy_compare.sex',new_param_file='copy_compare.param')
    
    MASK = True
    mask_file = os.path.join(curr_dir,'AstroImages','Masks',imgName1.replace('.fits','_mask.fits').replace('AstroImages/Good/',''))    
    if MASK and os.path.exists(mask_file):        
        fig.reconfigure('WEIGHT_IMAGE',mask_file)
        fig.reconfigure('WEIGHT_TYPE','MAP_WEIGHT')
        #pass
    fig.write_config_file(new_config_file='copy_compare.sex',new_param_file='copy_compare.param')
    
    
    print('Detecting from {}, measuring from {}').format(imgName1,imgName1)
    call_sex(imgName1+','+imgName1,config_file='copy_compare.sex')
    
    fig = sex_config.configure(imgName1+','+imgName2,'default.sex','default.param',dual=True)
    fig.default_config()
    
    mask_file = os.path.join(curr_dir,'AstroImages','Masks',imgName1.replace('.fits','_mask.fits').replace('AstroImages/Good/',''))
    if MASK and os.path.exists(mask_file):      
        fig.reconfigure('WEIGHT_IMAGE',mask_file)
        fig.reconfigure('WEIGHT_TYPE','MAP_WEIGHT')
        #pass
    fig.write_config_file(new_config_file='copy_compare.sex',new_param_file='copy_compare.param')

    print('Detecting from {}, measuring from {}').format(imgName1,imgName2)
    call_sex(imgName1+','+imgName2,config_file='copy_compare.sex')


#testImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
#testImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'

#compare(testImage1,testImage2)
