# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:46:50 2016

@author: charlesgulian
"""

# Do multiple linear regression

import os
import numpy as np
import matplotlib.pyplot as plt

def LSS(X,Y):
    A = np.ones(np.size())    
    
    A = [ones(size(X,1),1) X];
x = inv(A'*A)*A'*Y;