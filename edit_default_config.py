# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 10:12:21 2016

@author: charlesgulian
"""

import numpy as np
import matplotlib.pyplot as plt

class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
        self.masked = None

# Default mask file
maskFile = 'AstroImages/Masks/fpC-6484-x4078-y134_stitched_alignCropped_mask.fits'

# Good images    
goodImage1 = Image('AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits','Good','1')
goodImage2 = Image('AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits','Good','2')
goodImage3 = Image('AstroImages/Good/fpC-4868-x4211-y138_stitched_alignCropped.fits','Good','3')
goodImage4 = Image('AstroImages/Good/fpC-6383-x5176-y121_stitched_alignCropped.fits','Good','4')

# Bad images
badImage1 = Image('AstroImages/Bad/fpC-5759-x24775-y300_stitched_alignCropped.fits','Bad','1')
badImage2 = Image('AstroImages/Bad/fpC-6548-x24940-y302_stitched_alignCropped.fits','Bad','2')
badImage3 = Image('AstroImages/Bad/fpC-5781-x25627-y293_stitched_alignCropped.fits','Bad','3')
badImage4 = Image('AstroImages/Bad/fpC-7140-x24755-y270_stitched_alignCropped.fits','Bad','4')

# Co-added images
coaddedImage1 = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','_SDSS')
coaddedImage2 = Image('AstroImages/Coadd/custom_coadd_median.fits','Coadded','_Custom_Median')
coaddedImage3 = Image('AstroImages/Coadd/custom_coadd_mean.fits','Coadded','_Custom_Mean')

import sex_config
fig = sex_config.configure(coaddedImage2,'default.sex','default.param')
fig.write_config_file(new_config_file='default_copy.sex',new_param_file='param_copy.param')
