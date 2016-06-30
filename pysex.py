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

imgName = os.path.join(curr_dir,'AstroImages/','fpC-007202-r2-0152.fits')

def _check_call(args):
    try:
        subprocess.check_call(args)
    except subprocess.CalledProcessError:
        print 'Error: sex_call CalledProcessError'
        pass # handle errors in the called executable
    except OSError:
        print 'Error: OSError'
        pass # executable not found

def call_sex(imgName,config_file=None,args_ext=[]):
    # Call SExtractor
    args = [sex_path,imgName] 
    
    if not os.path.exists(imgName):
        imgName = os.path.join(curr_dir,'AstroImages/',imgName)
        args[1] = imgName
        
    if config_file != None:
        if not os.path.exists(config_file):
            config_file = os.path.join(curr_dir,config_file)
        args.extend(['-c',config_file])
        
    args.extend(args_ext)
    _check_call(args)

def loop_sex(img_dir):
    # Python version of sextractor_loop.sh
    if not os.path.exists(img_dir):
        img_dir = os.path.join(curr_dir,img_dir)
    # List names of .FITS files in img_dir
    imgNames = glob.glob(img_dir+'*.fits')

    for imgName in imgNames:
        
        _check_call(['cp','default.sex','copy.sex'])
        _check_call(['cp','default.param','copy.param'])
        
        # Do configuration with ./do_config.py imgName copy.sex copy.param
        # Rewrite this section after modifying do_config.py
        do_config_args = [os.path.join(curr_dir,'do_config.py'),imgName,'copy.sex','copy.param']
        _check_call([do_config_args])
        
        call_sex(imgName,config_file='copy.sex')
        
    

'''
#echo > img_name.txt
# Clear current contents of img_name
for f in AstroImages/*.fits
do
        # Copy default versions of configuration file and output catalog parameter file
        cp default.sex copy.sex
        cp default.param copy.param

        # Rewrite configuration and output catalog parameters in Python script do_config.py
        ./do_config.py $f copy.sex copy.param

        ######
        #echo $FNAME >> img_name.txt
        # Append filename to img_name.txt
        ######

        sex $f -c copy.sex
        # Run SExtractor
        #exit
done

#./create_regfile.py
'''