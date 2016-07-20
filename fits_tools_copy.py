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
    - Compute flux of object with image coordinates (x0,y0) within circular aperture
        of given radius r
    - Find all pixels within r pixels of (x0,y0) and sum
    '''
    # Defining distance function, dist:
    dist = lambda x,y: np.sqrt((x-x0)**2 + (y-y0)**2)     
    
    # Defining search region with xmin, xmax, ymin, ymax
    xmin = int(np.floor(x0-radius))
    xmax = int(np.ceil(x0+radius))
    ymin = int(np.floor(y0-radius))
    ymax = int(np.floor(y0+radius))
    
    print xmin,ymin
    print xmax,ymax
    
    search_region = image[xmin:xmax,ymin:ymax]
    aperture = np.zeros(np.shape(search_region))
    for x in range(xmin,xmax):
        for y in range(ymin,ymax):
            if dist(x,y) <= radius:
                aperture[x-xmin,y-ymin] = search_region[x-xmin,y-ymin]
            '''    
            elif (dist(x,y) > radius) and (dist(x,y) <= radius + 1.0):
                frac = 1.0 - (dist(x,y) - radius)
                print frac
                aperture[x-xmin,y-ymin] = frac*search_region[x-xmin,y-ymin]
            '''
    plt.pcolormesh(search_region,cmap=matplotlib.cm.gray)
    plt.colorbar()
    plt.show()
    plt.close()
    plt.pcolormesh(search_region - aperture,cmap=matplotlib.cm.jet)
    plt.colorbar()
    plt.show()
    return np.sum(aperture)
    
    
def computeObjectFlux2(x0,y0,radius,image):
    
    position = (x0,y0)
    aperture = CircularAperture(position,r=radius)
    return aperture_photometry(image,aperture)[0][0]
    
    
# getPixels() can be replaced by fits.getdata() (I did not know this)
def getPixels(image_file):
    hdulist = fits.open(image_file)
    data = hdulist[0].data
    hdulist.close()
    return data
    
    
def applyMask(image,mask,imageBias=0.0):
    '''
    - Apply a binary mask to an array
    '''
    image -= imageBias
    masked_image = np.multiply(image,mask)
    return masked_image
    

def maskImage(image_file,mask_file,masked_image_file=None,imageBias=0.0,Return=False):
    '''
    - Takes a .fits image file and .fits binary mask file as input
    - Applies binary mask to .fits image data
    - Rewrites masked image to new .fits file (masked_image_file)
    '''
    
    image = fits.getdata(image_file)
    mask = fits.getdata(mask_file)
    # Correct for image bias:
    image -= imageBias    
    masked_image = applyMask(image,mask,imageBias)
    
    if masked_image_file == None:
        masked_image_file = image_file.replace('.fits','_masked.fits').replace('Good','MaskedImages').replace('Bad','MaskedImages')
    
    fits.writeto(masked_image_file,masked_image,fits.getheader(image_file),clobber=True)
    
    if Return:
        return masked_image
        
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


    
    
test_image_file = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
test_image = fits.getdata(test_image_file)

x,y = 848.40234,727.06002
radius = 10.0

#flux1 = computeObjectFlux(x,y,radius,test_image)
flux2 = computeObjectFlux2(x,y,radius,test_image)
print flux2

#testMask = 'AstroImages/Masks/fpC-6484-x4078-y134_stitched_alignCropped_mask.fits'
#maskImage(testImage,testMask)