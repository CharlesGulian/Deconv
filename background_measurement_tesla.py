# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 18:45:14 2017

@author: charlesgulian
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fits_tools_tesla as fits_tools
from astropy.io import fits
from scipy.stats import linregress as linreg
import glob
import os

os.chdir('/home/cgulian2/Deconv')
curr_dir = os.getcwd()
new_dir = os.path.join(curr_dir,'Figures','Jan20')

# ===============================================================================
# Get images from directory

# Image directory:
image_dir = '/home/DATA/charlie/STRIPE82_330-360_AlignCropped/test7'
image_files = glob.glob(os.path.join(image_dir,'*alignCropped.fits'))
#image_files.remove('/home/DATA/charlie/STRIPE82_330-360_AlignCropped/test7/fpC-4927-x4127-y118_stitched_alignCropped.fits')
    
# ===============================================================================
# Mask sources; measure background

def create_mask(imageData):
    
    # Get image data, create blank (zeros) array
    mask = np.ones(np.shape(imageData))
    
    # Define indices where image pixel values are greater than 5-sigma threshold
    inds = np.where(imageData >= np.mean(imageData) + 0.05*np.std(imageData))
    # Set mask = 1.0 at these indices, save mask to new file
    mask[inds] = 0.0
    return mask


for m in range(len(image_files)):

    image_file = image_files[m]
    image_tag = os.path.split(image_file)[1]
    tag = image_tag[0:8]
    
    imageData = fits.getdata(image_file)
    header = fits.getheader(image_file)
    try:
        flux20 = header['flux20']
        sky = header['sky']
        softbias = header['softbias']
    except KeyError:
        print 'Error: Missing keyword in image header'
	continue

    alpha = 1e-8/(flux20)
    flux = alpha*(imageData - softbias - sky)
    imageData = flux

    mask = create_mask(imageData)
    temp = imageData*mask
    imageData = temp
    
    imageSum = np.sum(imageData)
    M = 1
    N_list = [16,64,800]
    linewidths = [12.0,5.0,0.5]
    alphas = [0.3,0.4,0.5]
    colors = ['y','g','b']
    for k in range(len(N_list)):
        N = N_list[k]
        imgBinDict = fits_tools.binImage(imageData,M,N)
        points = np.zeros(N)
        y_values = np.zeros(N)
        for i in range(N):
            imgBin = imgBinDict[M,i+1]
            #point = (np.sum(imgBin)*N)/imageSum
            #points[i] = point
            point = np.sum(imgBin)/len(np.where(imgBin != 0.0)[0])
            points[i] = point

            y_value = (1600./N)*(i+0.5)
            y_values[i] = y_value
            
        plt.plot(y_values,points,color=colors[k],alpha=alphas[k],linewidth=linewidths[k])
    plt.axis([0.0,1600.0,0.5,1.5])
    plt.title('Normalized Background Value vs. Y-coordinate')
    plt.xlabel('Image Y-coordinate (pixels)')
    plt.ylabel('Normalized Background Value')
    slope, intercept, r_value, p_value, std_err = linreg(y_values,points)
    plt.plot(y_values,slope*y_values + intercept,'r--',linewidth=2.0)
    plt.figtext(.45,.85,'p-value = {}'.format(p_value))
    plt.savefig(os.path.join(new_dir,'Background','Data','normalized_background_value_vs_y_linreg_{}.png'.format(image_tag)))
    plt.close()
    
    print p_value
    
    # Save background image (image data * mask)
    fits.writeto(os.path.join(new_dir,'Background','Images','background_image_{0}.fits'.format(image_tag)),imageData,clobber=True)
