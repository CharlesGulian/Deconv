# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:34:39 2016

@author: charlesgulian
"""
import os
curr_dir = os.getcwd()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pysex
import sex_stats
import scipy.stats as spst
#import sex_config
#import do_config

# Image deconvolution project:
# Main script for data analysis, image comparison, photometric statistics, and more

# Good image comparison
goodImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
goodImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'
goodImage3 = 'AstroImages/Good/fpC-4868-x4211-y138_stitched_alignCropped.fits'
goodImage4 = 'AstroImages/Good/fpC-6383-x5176-y121_stitched_alignCropped.fits'

goodImgs = [goodImage1,goodImage2,goodImage3,goodImage4]

# Bad image comparison:
badImage1 = 'AstroImages/Bad/fpC-5759-x24775-y300_stitched_alignCropped.fits' # Modest gradient from top to bottom
badImage2 = 'AstroImages/Bad/fpC-6548-x24940-y302_stitched_alignCropped.fits' # Modest gradient from top to bottom
badImage3 = 'AstroImages/Bad/fpC-5781-x25627-y293_stitched_alignCropped.fits' # Very weak gradient from bottom left to top right
badImage4 = 'AstroImages/Bad/fpC-7140-x24755-y270_stitched_alignCropped.fits' # Weak gradient from bottom left to top right

badImgs = [badImage1,badImage2,badImage3,badImage4]

for testImage1 in goodImgs:
    for testImage2 in goodImgs:
        
        if testImage1 == testImage2:
            continue
        
        # Images to compare 
        #testImage1 = badImage1
        #testImage2 = badImage2
        
        
        output = pysex.compare(testImage1,testImage2) # (Implement/uncomment to create new comparison file)
        
        img_tag1 = (os.path.split(testImage1)[1])
        img_tag1 = img_tag1[0:len(img_tag1)-len('.fits')]
        img_tag2 = (os.path.split(testImage2)[1])
        img_tag2 = img_tag2[0:len(img_tag2)-len('.fits')]
        
        outputCat1 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag1+'_compare.cat')
        if not os.path.exists(outputCat1):
            print 'Error: first output catalog path does not exist'
        
        outputCat2 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag2+'_compare.cat')
        if not os.path.exists(outputCat2):
            print 'Error: second output catalog path does not exist'
            
        # Getting image data (pixel values for testImage1 and testImage2)
        fluxPixels1,fluxPixels2 = sex_stats.getPixelValues(testImage1),sex_stats.getPixelValues(testImage2)
        # Correcting for image bias = 1000        
        pixelBias = 1000
        fluxPixels1 -= pixelBias
        fluxPixels2 -= pixelBias
        
        fluxRatio_pixelWise = np.divide(fluxPixels1,fluxPixels2)
        
        m,n = 128,128
        fluxRatioBins = sex_stats.binImage(fluxRatio_pixelWise,M=m,N=n)
        
        fluxRatioBin_Avgs = np.zeros([m,n])
        for i in range(m):
            for j in range(n):
                fluxRatioBin_Avgs[i,j] = np.mean(fluxRatioBins[i,j])
        
        # Masking NaNs in fluxRatioBin_Avgs with red squares
        fluxRatioBin_Avgs_Masked = np.ma.array(fluxRatioBin_Avgs,mask=np.isnan(fluxRatioBin_Avgs))
        cmap = matplotlib.cm.gray
        cmap.set_bad('r',1.)
        
        plt.pcolormesh(fluxRatioBin_Avgs,cmap=cmap,vmin=np.nanmean(fluxRatioBin_Avgs)-7.5*np.nanstd(fluxRatioBin_Avgs),vmax=np.nanmean(fluxRatioBin_Avgs)+7.5*np.nanstd(fluxRatioBin_Avgs))
        plt.colorbar()
        plt.axis([0,m,0,n])
        plt.xlabel('X Bin')
        plt.ylabel('Y Bin')
        plt.title('Flux Ratio Bin Averages: {} x {}'.format(m,n))
        if not os.path.exists(os.path.join(curr_dir,'Figures','Jul14Imgs','ImgBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10]))):
                os.mkdir(os.path.join(curr_dir,'Figures','Jul14Imgs','ImgBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10])))
        plt.savefig(os.path.join(curr_dir,'Figures','Jul14Imgs','ImgBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10]),'fluxRatioBin_Avgs_{}x{}.png'.format(m,n)))
        plt.close()




""" THIS SECTION OF CODE WAS COMMENTED OUT ON July 12th, 2016; uncomment to do statistical analysis

chiSqNorm_linear = []
chiSqNorm_flat = []
rSqAdj = []
numPoints = []

for j in range(len(maxSig)): 
    
    # Clipping data
    fluxRatio_meanSubtracted_sigmaClipped = []
    fluxRatio_excess = []
    fluxAvg_clip = []
    x_clip,y_clip = [],[]
    x_exc,y_exc = [],[]
    for i in range(len(fluxRatio_meanSubtracted)):
        if np.abs(fluxRatio_meanSubtracted[i]) < maxSig[j]*fluxRatio_std:
            fluxRatio_meanSubtracted_sigmaClipped.append(fluxRatio_meanSubtracted[i])
            x_clip.append(x[i])
            y_clip.append(y[i])
            fluxAvg_clip.append(fluxAvg[i])
            
        else:
            fluxRatio_excess.append(fluxRatio_meanSubtracted[i])
            x_exc.append(x[i])
            y_exc.append(y[i])
    
    fluxRatio_meanSubtracted_sigmaClipped,fluxRatio_excess = np.array(fluxRatio_meanSubtracted_sigmaClipped),np.array(fluxRatio_excess)
    x_clip,y_clip,x_exc,y_exc = np.array(x_clip),np.array(y_clip),np.array(x_exc),np.array(y_exc)    
    numPoints.append(float(len(x_clip)))
    
    # Analyzing goodness-of-fit of 3D linear model fitted to data:

    coeffs = sex_stats.linReg3D(x_clip,y_clip,fluxRatio_meanSubtracted_sigmaClipped)[0]
    
    linearModelPoints = coeffs[0] + coeffs[1]*x_clip + coeffs[2]*y_clip
    flatModelPoints = np.ones(np.shape(fluxRatio_meanSubtracted_sigmaClipped))*fluxRatio_mean
    
    # SciPy: scipy.stats.chisquare

    #CSN_lin = spst.chisquare()    
    
    CSN_lin = sex_stats.chiSquareNormalized(fluxRatio_meanSubtracted_sigmaClipped,linearModelPoints,3)
    CSN_flat = sex_stats.chiSquareNormalized(fluxRatio_meanSubtracted_sigmaClipped,flatModelPoints,1)
    RSA = sex_stats.rSquaredAdjusted(fluxRatio_meanSubtracted_sigmaClipped,linearModelPoints,3)
    
    chiSqNorm_linear.append(CSN_lin)
    chiSqNorm_flat.append(CSN_flat)
    rSqAdj.append(RSA)
    
    plot = True # Warning: do not change to true unless length of maxSig small
    if plot:
        # Plotting source-wise flux ratio w/ colors
        plt.scatter(x_clip, y_clip, s=25*np.log10(0.1*np.array(fluxAvg_clip)), c=fluxRatio_meanSubtracted_sigmaClipped, vmin=-1.5*maxSig[j]*fluxRatio_std, vmax=1.5*maxSig[j]*fluxRatio_std, alpha=0.75)        
        plt.axis([0,1600,0,1600])
        plt.colorbar()
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')
        plt.title('Flux Ratio Color Map: sigma cutoff = '+str(maxSig[j])[0:4])
        plt.savefig((curr_dir+'/Figures/{}_{}_maxSig{}_fluxRatio_LINETEST.png'.format(img_tag1, img_tag2, str(maxSig[j])[0:4])))
        plt.close()
    
    hist = False # Warning: do not change to true unless length of maxSig small
    if hist:
        # Plotting histogram of flux ratio
        plt.hist(fluxRatio_meanSubtracted_sigmaClipped,bins=20,color='green')
        plt.title('Histogram of Flux Ratio')
        plt.ylabel(('Mean subtracted + clipped @ {} sigma').format(str(maxSig[j])[0:4]))
        plt.xlabel('Flux ratio')
        plt.savefig((curr_dir+'/Figures/Hist_{}_{}_maxSig{}_fluxRatio.png'.format(img_tag1, img_tag2, str(maxSig[j])[0:4])))
        plt.close()

# Changing lists to NumPy arrays:
# chiSqNorm_linear,chiSqNorm_flat,rSqAdj = np.array(chiSqNorm_linear),np.array(chiSqNorm_flat),np.array(rSqAdj)

# Number of data points analyzed:
numPoints = np.array(numPoints)
numPoints = numPoints*(1.0/float(len(fluxRatio)))

# Plotting reduced chi-square statistic
plt.close()
plt.plot(maxSig,chiSqNorm_linear,'r-',label='Linear model')
plt.plot(maxSig,chiSqNorm_flat,'b-',label='Flat model')
plt.plot(maxSig,numPoints,'0.35',label='Frac. of data points')
plt.legend()
plt.axis([-0.1,1.0,0.0,3.0])
plt.title('Normalized Chi-square vs. Sigma Cutoff: 3D Linear Eq. + Gaussian Noise Test')
plt.ylabel('Normalized Chi-square: Linear')
plt.xlabel('Sigma Cutoff  (# standard deviations from mean)')
plt.ylabel('Normalized Chi-square')
#plt.savefig(os.path.join(os.getcwd(),'Figures','StatAnalysis','Linear_eq_test_6'))
plt.show()

'''
# Plotting adjusted r-squared statistic

plt.plot(maxSig,rSqAdj,'k-')
plt.axis([0.0,1.0,-1.1,1.1])
plt.title('Adjusted r-Squared vs. Sigma Cutoff')
plt.xlabel('Sigma Cutoff  (# standard deviations from mean)')
plt.ylabel('Adjusted r-Squared')
plt.show()
'''
"""



''' # Testing with artificial data: random noise and linear equations with noise
# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:

#temp = np.std(flux1)*np.random.randn(flux1.size) + np.mean(flux1)
#flux1 = temp
#flux1 = np.mean(flux1)*np.ones(flux1.size) + np.random.randn(flux1.size)
#temp = (np.mean(flux2)/np.mean(flux1))*flux1 + np.std(flux2)*np.random.randn(flux2.size)
#flux2 = temp
#flux2 += 10.0*(x + y - np.mean(x) - np.mean(y))


# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:# TESTING:
#'''