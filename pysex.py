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
        pass # handle errors in the called executable
    except OSError:
        print 'Error: OSError'
        pass # executable not found

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
        
        #print 'Copying configuration and output catalog parameter files'
        _call(['cp','default.sex','copy.sex'])
        _call(['cp','default.param','copy.param'])
        
        # Do configuration with ./do_config.py imgName copy.sex copy.param
        # Rewrite this section after modifying do_config.py
        
        #print 'Configuring with do_config.py'
        do_config_args = [os.path.join(curr_dir,'do_config.py'),imgName,'copy.sex','copy.param']
        _call(do_config_args)
        
        print ('Running SExtractor on current image: {}').format(imgName)
        call_sex(imgName,config_file='copy.sex')
        
        
def compare(imgName1,imgName2):
    
    #print 'Copying configuration and output catalog parameter files'
    _call(['cp','default.sex','copy_compare.sex'])
    _call(['cp','default.param','copy_compare.param'])
    
    # Do configuration with ./do_config.py imgName copy.sex copy.param
    # Rewrite this section after modifying do_config.py
    #print 'Configuring with do_config.py'
    do_config_args = [os.path.join(curr_dir,'do_config.py'),imgName1+','+imgName1,\
                    'copy_compare.sex','copy_compare.param','-d']
    _call(do_config_args)
    
    print('Detecting from {}, measuring from {}').format(imgName1,imgName1)
    call_sex(imgName1+','+imgName1,config_file='copy_compare.sex')
    
    _call(['cp','default.sex','copy_compare.sex'])
    _call(['cp','default.param','copy_compare.param'])
    
    #print 'Configuring with do_config.py'
    do_config_args = [os.path.join(curr_dir,'do_config.py'),imgName1+','+imgName2,\
                    'copy_compare.sex','copy_compare.param','-d']
    _call(do_config_args)

    print('Detecting from {}, measuring from {}').format(imgName1,imgName2)
    call_sex(imgName1+','+imgName2,config_file='copy_compare.sex')


#testImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
#testImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'

#compare(testImage1,testImage2)