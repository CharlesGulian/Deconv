# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 12:24:37 2016

@author: charlesgulian
"""

import numpy as np
import matplotlib as pyplot
from astropy.io import fits

# Image co-add program
# Takes a list of N .FITS image files as input, outputs a new image whose individual pixels
# are equal to the median of the corresponding N pixels from the input images;
# i.e. let co-add be matrix C, kth image be I_k: C_ij = median({I_k_ij}) : k = 1,2,...,N

