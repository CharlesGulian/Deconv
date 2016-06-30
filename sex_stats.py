# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 20:31:36 2016

@author: charlesgulian
"""
# Library of statistics functions for photometric statistics and data analysis

import os
os.chdir('/Users/annepstein/Work/Deconv')
import numpy as np
import matplotlib.pyplot

# 3D linear least-squares regression analysis:
def linReg3D(x1,x2,y):
    A = np.column_stack((np.ones(x1.size),x1,x2))
    output = np.linalg.lstsq(A,y)
    return output
    
