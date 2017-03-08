# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 18:45:14 2017

@author: charlesgulian
"""

import numpy as np
import matplotlib.pyplot as plt
import fits_tools
from astropy.io import fits
from scipy.stats import linregress as linreg
import os

# ===============================================================================

class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
        self.masked = None

# Bad images
badImage1 = Image('AstroImages/Bad/fpC-5759-x24775-y300_stitched_alignCropped.fits','Bad','1')

# Co-added images
coaddedImage1 = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','_SDSS')
coaddedImage2 = Image('AstroImages/Coadd/custom_coadd_median.fits','Coadded','_Custom_Median')
coaddedImage3 = Image('AstroImages/Coadd/custom_coadd_mean.fits','Coadded','_Custom_Mean')
coaddedImage4 = Image('AstroImages/Coadd/fpC-4874-x1024-y65466_cutout.fits','Coadded','_Sparse1')
coaddedImage5 = Image('AstroImages/Coadd/fpC-4874-x1054-y53218_cutout.fits','Coadded','_Sparse2')
coaddedImage6 = Image('AstroImages/Coadd/fpC-4874-x1072-y77946_cutout.fits','Coadded','_Sparse3')
coaddedImage7 = Image('AstroImages/Coadd/fpC-4874-x1046-y47134_cutout.fits','Coadded','_Crowded1')

# Deconvolved image
deconvolvedImage = Image('AstroImages/Deconvolved/deconv.fits','Deconvolved','Deconvolved')
# Trimming deconvolved image to 1600x1600
deconvolvedImageData = fits_tools.getPixels(deconvolvedImage.filename)
TRIM = True
if (np.shape(deconvolvedImageData) != (1600,1600)) and TRIM:
    print 'Trimming deconvolved image'
    padSize = np.array(list(np.shape(deconvolvedImageData))) - np.array([1600,1600])
    deconvolvedImageData_trimmed = deconvolvedImageData[padSize[0]/2:np.shape(deconvolvedImageData)[0]-padSize[0]/2,
                                                        padSize[1]/2:np.shape(deconvolvedImageData)[1]-padSize[1]/2]
    temp = (deconvolvedImage.filename).replace('.fits','_trimmed.fits')
    fits.writeto(temp,deconvolvedImageData_trimmed,fits.getheader(deconvolvedImage.filename),clobber=True)
    deconvolvedImage = Image(temp,'Deconvolved','Deconvolved')
    
# ===============================================================================

def create_mask(imageData):
    
    # Get image data, create blank (zeros) array
    mask = np.ones(np.shape(imageData))
    
    # Define indices where image pixel values are greater than 5-sigma threshold
    inds = np.where(imageData >= np.mean(imageData) + 0.1*np.std(imageData))
    # Set mask = 1.0 at these indices, save mask to new file
    mask[inds] = 0.0
    return mask

image_list = [deconvolvedImage]
for i in range(len(image_list)):
    image = image_list[i]
    # Original code:
    imageData = fits.getdata(image.filename)
    #imageData -= 1000.0
    
    mask = create_mask(imageData)
    #temp = imageData*mask
    #imageData = temp
    
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
            #print np.where(imgBin < 0.0)
            #point = (np.sum(imgBin)*N)/imageSum
            point = np.sum(imgBin)/len(np.where(imgBin != 0.0)[0])
            print ''
            #print len(np.where(imgBin != 0.0)[0])
            #print len(np.where(imgBin != 0.0))
            #print np.shape(np.where(imgBin != 0.0))
            #print np.shape(imgBin)
            
            nonzero_inds = np.where(imgBin != 0.0)
            nonzero_imgBin = imgBin[nonzero_inds]
            
            print 'Mean = ',np.mean(nonzero_imgBin)
            print 'Standard Dev = ',np.std(nonzero_imgBin)
            print 'Max = ',np.max(nonzero_imgBin)
            print 'Min = ',np.min(nonzero_imgBin)
        
            points[i] = point
        
            y_value = (1600./N)*(i+0.5)
            y_values[i] = y_value
            
        plt.plot(y_values,points,color=colors[k],alpha=alphas[k],linewidth=linewidths[k])
        #plt.plot(y_values,1000.0*np.ones(len(y_values)),'k',linewidth=2.0)
        plt.ylim(0.0,100.0)
    plt.title('Local Background Value vs. Y-coordinate')
    plt.xlabel('Image Y-coordinate (pixels)')
    plt.ylabel('Local Mean Background Value')
    
    
    # ===============================================================================
    # Save data
    curr_dir = os.getcwd()
    new_dir = os.path.join(curr_dir,'Figures','Jan25')
    #plt.savefig(os.path.join(new_dir,'normalized_background_value_vs_y.png'))
    
    # Create and save new plot with best-fit line
    slope, intercept, r_value, p_value, std_err = linreg(y_values,points)
    plt.plot(y_values,slope*y_values + intercept,'r--',linewidth=2.0)
    plt.savefig(os.path.join(new_dir,'Background','Data','mean_background_value_vs_y_linreg{}.png'.format(image.ID)))
    plt.close()
    
    print p_value
    
    # Save background image (image data * mask)
    fits.writeto(os.path.join(new_dir,'Background','Images','background_image{0}.fits'.format(image.ID)),imageData,clobber=True)
