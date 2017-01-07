# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 21:18:54 2016

@author: charlesgulian
"""

import os
os.chdir('/Users/annepstein/Work/Deconv')
curr_dir = os.getcwd()
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib

from photutils import aperture_photometry
from photutils import CircularAperture

def binImage(pixelArray,M=3,N=3):
    '''
    - Bins pixels along image axes into MxN bins (default MxN = 3x3)  
    '''
    pixels = pixelArray
    imgDim1,imgDim2 = np.shape(pixels)
    xBinSize,yBinSize = float(imgDim1)/float(M),float(imgDim2)/float(N)
    
    imgBinDict = {} # Dictionary for storing 
    #print xBinSize,yBinSize
    
    for i in range(M):
        for j in range(N):
            imgBinDict[i,j] = pixels[int(np.ceil(i*xBinSize)):int(np.floor((i+1)*xBinSize)),\
                                    int(np.ceil(j*yBinSize)):int(np.floor((j+1)*yBinSize))]
            #print ''
            #print 'Bin: ',i,j
            #print 'Shape: ',np.shape(imgBinDict[i,j])
    return imgBinDict
    
def computeObjectFlux(x0,y0,radius,image):
    '''
    - Compute total flux within circular aperture of the given radius 
        from source at image coordinates (x0,y0)
    '''
    position = (x0,y0)
    aperture = CircularAperture(position,r=radius)
    return aperture_photometry(image,aperture)[0][0]
    
    
# getPixels() can be replaced by fits.getdata() (I did not know this)
def getPixels(image_file,delete=False):
    hdulist = fits.open(image_file)
    data = hdulist[0].data
    hdulist.close()
    if delete:
        del hdulist[0].data
    return data
    
    
def applyMask(image,mask):
    '''
    - Apply a binary mask to an array
    '''
    masked_image = np.multiply(image,mask)
    return masked_image
    

def maskImage(image_file,mask_file,masked_image_file=None,Return=False):
    '''
    - Takes a .fits image file and .fits binary mask file as input
    - Applies binary mask to .fits image data
    - Rewrites masked image to new .fits file (masked_image_file)
    '''
    
    image = fits.getdata(image_file)
    mask = fits.getdata(mask_file)
    masked_image = np.multiply(image,mask)
    
    inds = np.where(masked_image == 0.0)
    masked_image[inds] += 1e-12 # Prevent NaNs
    
    if masked_image_file == None:
        masked_image_file = image_file.replace('.fits','_masked.fits').replace('Good','MaskedImages').replace('Bad','MaskedImages')
    
    fits.writeto(masked_image_file,masked_image,fits.getheader(image_file),clobber=True)
    
    if Return:
        return masked_image

    
def shift_image(image,x_offset,y_offset):
    
    # Shifts image pixels from (x,y) to (x-x_offset),(y-y_offset)
    
    dims = np.shape(image) # Image dimensions
    dim1,dim2 = dims[0],dims[1]
    
    blank = np.zeros(dims) + 1e-8 # Define blank array to receive new image data
    shifted_image = blank
    
    dy,dx = x_offset,y_offset # These are intentionally reversed
    for i in range(dim1):
        for j in range(dim2):
            if (i+dx < dim1) and (i+dx >= 0) and (j+dy < dim2) and (j+dy >= 0):
                shifted_image[i,j] = image[i+dx,j+dy] # Why does this work?
    
    return shifted_image


def subtractBias(image_file,new_image_file=None,bias=0.0,Return=False):
    '''
    - Takes a .fits image file as input
    - Subtracts median from image data, writes new data to new image file (new_image_file)
    '''
    if new_image_file == None:
        new_image_file = image_file
    
    image = fits.getdata(image_file)
    image -= bias
    
    fits.writeto(new_image_file,image,fits.getheader(image_file),clobber=True)
    
    if Return:
        return image
        
        
def subtractMedian(image_file,new_image_file=None,Return=False):
    '''
    - Takes a .fits image file as input
    - Subtracts median from image data, writes new data to new image file (new_image_file)
    '''
    if new_image_file == None:
        new_image_file = image_file
    
    image = fits.getdata(image_file)
    image -= np.median(image)
    
    fits.writeto(new_image_file,image,fits.getheader(image_file),clobber=True)
    
    if Return:
        return image
        
def write_pixel_offset(x_offset,y_offset,image_file,new_image_file=None):
    # Add (x,y) pixel offset to .FITS header of an image

    header = fits.getheader(image_file) # Get original .FITS header
    header['x_offset'] = x_offset # Add new keywords and values to header
    header['y_offset'] = y_offset
    # If no new image file specified, default writes new header to original image header
    if new_image_file == None:
        new_image_file = image_file
    # Write header to new image
    fits.writeto(new_image_file,fits.getdata(image_file),header,clobber=True)

  
'''
# Testing:

test_image_file = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
test_image = fits.getdata(test_image_file)

catalog = 'Results/fpC-6484-x4078-y134_stitched_alignCropped_fpC-6484-x4078-y134_stitched_alignCropped_compare.cat'
import sex_stats
fig = sex_stats.data(catalog)
x = fig.get_data('X_IMAGE')
y = fig.get_data('Y_IMAGE')

xlow = np.where(x > 651.0)
xhigh = np.where(x < 658.9)
xin = np.intersect1d(xlow,xhigh)
ylow = np.where(y > 820.0)
yhigh = np.where(y < 826.0)
yin = np.intersect1d(ylow,yhigh)

obj = np.intersect1d(xin,yin)

DATA = fig.Data

x,y = 848.39102,727.23274
radius = 10.

flux = computeObjectFlux(x,y,radius,test_image)
print flux


#testMask = 'AstroImages/Masks/fpC-6484-x4078-y134_stitched_alignCropped_mask.fits'
#maskImage(testImage,testMask)
'''