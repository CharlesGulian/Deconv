# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 23:51:07 2016

@author: charlesgulian
"""

import os
import subprocess

curr_dir = os.getcwd()
#script = os.path.join(curr_dir,'sextractor_loop.sh')
imgName = os.path.join(curr_dir,'AstroImages/','fpC-007202-r2-0152.fits')
print imgName
try:
    subprocess.check_call(['sex',imgName])
except subprocess.CalledProcessError:
    pass # handle errors in the called executable
except OSError:
    pass # executable not found