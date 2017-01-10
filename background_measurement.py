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

# Co-added images
coaddedImage1 = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','_SDSS')
coaddedImage2 = Image('AstroImages/Coadd/custom_coadd_median.fits','Coadded','_Custom_Median')
coaddedImage3 = Image('AstroImages/Coadd/custom_coadd_mean.fits','Coadded','_Custom_Mean')

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
    inds = np.where(imageData >= np.mean(imageData) + 0.05*np.std(imageData))
    # Set mask = 1.0 at these indices, save mask to new file
    mask[inds] = 0.0
    return mask

image = coaddedImage1
imageData = fits.getdata(image.filename)
mask = create_mask(imageData)
temp = imageData*mask
imageData = temp

imageSum = np.sum(imageData)
M = 1
N_list = [40,160,1600]
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
        point = (np.sum(imgBin)*N)/imageSum
        points[i] = point
    
        y_value = (1600./N)*(i+0.5)
        y_values[i] = y_value
        
    plt.plot(y_values,points,color=colors[k],alpha=alphas[k],linewidth=linewidths[k])
plt.title('Normalized Background Value vs. Y-coordinate')
plt.xlabel('Image Y-coordinate (pixels)')
plt.ylabel('Normalized Background Value')


# ===============================================================================
# Save data
curr_dir = os.getcwd()
new_dir = os.path.join(curr_dir,'Figures','Jan9')
#plt.savefig(os.path.join(new_dir,'normalized_background_value_vs_y.png'))

# Create and save new plot with best-fit line
slope, intercept, r_value, p_value, std_err = linreg(y_values,points)
plt.plot(y_values,slope*y_values + intercept,'r--',linewidth=2.0)
plt.savefig(os.path.join(new_dir,'Background','Data','normalized_background_value_vs_y_linreg{}.png'.format(image.ID)))
plt.close()

print p_value

# Save background image (image data * mask)
fits.writeto(os.path.join(new_dir,'Background','Images','background_image{0}.fits'.format(image.ID)),imageData,clobber=True)
