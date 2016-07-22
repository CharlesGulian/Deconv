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
from astropy.io import fits

import pysex
import sex_stats
import sex_config
import fits_tools

# Image deconvolution project:
# Main script for data analysis, image comparison, photometric statistics, and more


class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID

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

# Co-added image
coaddedImage = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','')

# Deconvolved image
deconvolvedImage = Image('AstroImages/Deconvolved/deconv.fits','Deconvolved','')
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
    deconvolvedImage = Image(temp,'Deconvolved','')
# ===============================================================================
# List of images to be compared: 
comparisonImages = [deconvolvedImage,coaddedImage,goodImage2,goodImage3]

for img1 in comparisonImages:
    for img2 in comparisonImages:
        if img1 == img2:
            continue        
        
        # Make directory to save results and figures to different comparison categories (Good vs. Good, Good vs. Deconvolved, etc.)
        image1,image2 = img1.filename,img2.filename # Get image filenames
        new_dir = os.path.join(curr_dir,'Figures','Jul21',img1.category+img1.ID+'_vs_'+img2.category+img2.ID)
        if not os.path.exists(new_dir):
            os.mkdir(os.path.join(new_dir))
        
        # Set medSub = True to subtract median from each comparison image
        medSub = False
        if medSub:
            image1_medSub = image1.replace('.fits','_medSub.fits')
            image2_medSub = image2.replace('.fits','_medSub.fits')
            
            fits_tools.subtractMedian(image1,new_image_file=image1_medSub)
            fits_tools.subtractMedian(image2,new_image_file=image2_medSub)
            
            image1 = image1_medSub
            image2 = image2_medSub
        
        # Subtract bias of 1000.0 from "Good" images
        if img1.category == 'Good':
            image1_biasSub = image1.replace('.fits','_biasSub.fits')
            fits_tools.subtractBias(image1,new_image_file=image1_biasSub,bias=1000.0)
            image1 = image1_biasSub
        if img2.category == 'Good':
            image2_biasSub = image2.replace('.fits','_biasSub.fits')
            fits_tools.subtractBias(image2,new_image_file=image2_biasSub,bias=1000.0)            
            image2 = image2_biasSub
            
        
        # =======================================================================
        # Write configuration files for SExtractor comparison
        
        # Writing configuration file for first image
        fig1 = sex_config.configure(image1+','+image1,'default.sex','default.param',dual=True)
        # Do default configuration        
        fig1.default_config()
        
        # Writing onfiguration file for second image
        fig2 = sex_config.configure(image1+','+image2,'default.sex','default.param',dual=True)
        # Do default configuration
        fig2.default_config()
        
        imgs = [img1,img2]
        figs = [fig1,fig2]
        for i,img in enumerate(imgs):
            if img.category == 'Deconvolved':
                print ''
                print 'Reconfiguring deconvolved image: DETECT_THRESH = 7.0'
                print 'Turning off background estimation/subtraction'
                print ''
                figs[i].reconfigure('DETECT_THRESH',7.0)
                figs[i].reconfigure('BACK_TYPE','MANUAL')
                figs[i].reconfigure('BACK_VALUE',0.0)
                                
                
                pass
            
            if img.category == 'Coadded':
                print ''
                print 'Subtracting median from co-added image: median = {0}'.format(np.median(fits_tools.getPixels(img.filename)))
                print ''
                temp = (img.filename).replace('.fits','_medSub.fits')
                fits_tools.subtractMedian(img.filename,new_image_file=temp)
                img = Image(temp,'Coadded','')
        
        fig1.write_config_file(new_config_file='copy_compare1.sex',new_param_file='copy_compare1.param')
        fig2.write_config_file(new_config_file='copy_compare2.sex',new_param_file='copy_compare2.param')
        
        # =======================================================================
        # Compare images in SExtractor
        
        pysex.compare(image1,image2,'copy_compare1.sex','copy_compare2.sex')
        
        # =======================================================================
        # Retrieve data from output catalogs        
        
        # Get image tags
        img_tag1 = (os.path.split(image1)[1])
        img_tag1 = img_tag1[0:len(img_tag1)-len('.fits')]
        img_tag2 = (os.path.split(image2)[1])
        img_tag2 = img_tag2[0:len(img_tag2)-len('.fits')]
        
        outputCat1 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag1+'_compare.cat')
        if not os.path.exists(outputCat1):
            print 'Error: first output catalog path does not exist'
        
        outputCat2 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag2+'_compare.cat')
        if not os.path.exists(outputCat2):
            print 'Error: second output catalog path does not exist'
        
        # Create sex_stats.data objects:
        img1data = sex_stats.data(outputCat1)
        img2data = sex_stats.data(outputCat2)
        
        # Create .reg files from output catalogs
        CREATE_regFiles = False
        if CREATE_regFiles:
            img1data.create_regFile()
            img2data.create_regFile()
    
    
        # =======================================================================
        # Flux ratio analysis:
        
        # Getting object positions, 0.90-light radii, and image data (pixels)
        x,y = img1data.get_data('X_IMAGE'),img1data.get_data('Y_IMAGE')
        rads1 = img1data.get_data('FLUX_RADIUS')
        rads2 = img2data.get_data('FLUX_RADIUS')
        # Handling negative radii (which result from noisy detection or measurement images, accroding to Emmanuel Bertin, creator of SExtractor)
        neg_inds1 = np.where(rads1 < 0.0)[0] # Where rads1 is negative
        neg_inds2 = np.where(rads2 < 0.0)[0] # Where rads2 is negative
        neg_inds = np.intersect1d(neg_inds1,neg_inds2) # Where both arrays are negative
        # Radii array = rads1 (default)   
        rads = rads1
        for i in neg_inds1:
            rads[i] = rads2[i] # rads = rads2 where rads1 is negative
        for i in neg_inds:
            rads[i] = np.mean(rads1[i]) # rads = mean(rads1) where both rads1 and rads2 are negative
        
        imageData1,imageData2 = fits_tools.getPixels(image1),fits_tools.getPixels(image2)
        # Set all pixel values >= 0.0:
        if len(np.where(imageData1 < 0.0)[0]) > 0:
            print ''
            print 'Warning: negative pixel values in first image'
            print 'Subtracting minimum value of image from pixel array'
            print 'Number of negative pixels: ', len(np.where(imageData1 < 0.0)[0])
            print 'Minimum pixel value: {} -> {}'.format(np.min(imageData1),0.0)
            imageData1 -= np.min(imageData1)
        if len(np.where(imageData2 < 0.0)[0]) > 0:
            print ''
            print 'Warning: negative pixel values in second image'
            print 'Subtracting minimum value of image from pixel array'
            print 'Number of negative pixels: ', len(np.where(imageData2 < 0.0)[0])
            print 'Minimum pixel value: {} -> {}'.format(np.min(imageData2),0.0)  
            imageData2 -= np.min(imageData2)
        # Computing circular-aperture fluxes
        flux1,flux2 = np.zeros(np.shape(x)),np.zeros(np.shape(x))
        for i in range(len(x)):
            flux1[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData1)
            flux2[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData2)

        # Compute object-wise flux ratio
        fluxRatio = np.divide(flux1,flux2)
        
        #''' 
        # Print statistics:
        print ''
        print img1.category+' vs. '+img2.category
        print ''
        print 'Number of objects detected: ',len(fluxRatio)
        print 'Minimum flux values: ', np.min(flux1),' ',np.min(flux2)
        print 'Number of negative flux values: ', len(np.where(flux1 < 0.0)[0]),' ',len(np.where(flux2 < 0.0)[0])
        print 'Mean flux values: ',np.mean(flux1),' ',np.mean(flux2)
        print ''
        print 'Minimum flux ratio: ', np.min(fluxRatio)
        print 'Maximum flux ratio: ', np.max(fluxRatio)
        print 'Standard deviation of flux ratio: ', np.std(fluxRatio)
        print ''
        print 'Minimum values of images: ', np.min(imageData1),' ',np.min(imageData2)
        print 'Median values of images: ', np.median(imageData1),' ',np.median(imageData2)
        print 'Mean values of images: ', np.mean(imageData1),' ',np.mean(imageData2)
        print ' '
        #'''
        
        # Creating histogram of flux1 and flux2
        # plt.hist has range=(tup1,tup2) to set upper and lower bounds of histogram
        plt.hist(flux1,bins=70,normed=False,color='green',label=img1.category)
        plt.hist(flux2,bins=70,normed=False,color='blue',label=img2.category)
        plt.legend()
        plt.title('Histogram of Object Flux for Entire Image')
        plt.ylabel('Frequency (N)')
        plt.xlabel('Object flux')
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'flux1_flux2_hist.png'))
            plt.close()
        else:
            plt.show()

        # Creating histogram of flux1/flux2 (object-wise flux ratio)
        #plt.hist(fluxRatio,bins=70,range=(-25.0,25.0),color='green') # Range = (-25.0,25.0)
        #plt.hist(fluxRatio,bins=int(abs(np.max(fluxRatio)-np.min(fluxRatio))*10.0),range=(np.min(fluxRatio),np.max(fluxRatio)),color='green') 
        plt.hist(fluxRatio,bins=70,range=(np.min(fluxRatio),np.max(fluxRatio)),color='green')            
        plt.title('Histogram of Object-wise Flux Ratio')
        plt.ylabel('Frequency (N)')
        plt.xlabel('Object-wise flux ratio')
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'fluxRatio_hist.png'))
            plt.close()
        else:
            plt.show()
            
        # Creating histogram of FLUX_RADIUS for PHOT_FLUXFRAC = 0.75,0.8,0.9
        plt.hist(rads,bins=70,color='green')
        plt.title('Histogram of 9/10-light radii')
        plt.xlabel('9/10-light radius')
        plt.ylabel('Frequency (N)')
        plt.savefig(os.path.join(new_dir,'hist_flux_radius0.90.png'))
        plt.close()
        
        # Creating color plot of object-wise flux ratio
        cmap = matplotlib.cm.jet
        plt.scatter(x,y,s=5.0*rads,c=fluxRatio-np.median(fluxRatio),marker='o',vmin=np.min(fluxRatio),vmax=np.max(fluxRatio),alpha=0.85)
        plt.axis([0,1600,0,1600])  
        plt.colorbar()
        plt.title('Map of Object-wise Flux Ratio')
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'fluxRatio_colorplot.png'))
            plt.close()
        else:
            plt.show()
        # Creating histogram of object-wise flux ratio for 4x4 bins
        m,n = 4,4
        xBins,yBins,fluxRatioBins = sex_stats.binData(x,y,fluxRatio,M=m,N=n)
        fluxRatioBin_means = np.zeros([m,n])
        for i in range(m):
            for j in range(n):
                plt.subplot(m,n,(n*i + (j+1)))
                plt.hist(fluxRatioBins[i,j],bins=70)
                plt.axis([0.0,7.5,0.0,23.0])
                fluxRatioBin_means[i,j] = np.mean(fluxRatioBins[i,j])
        
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'fluxRatioBins_hist.png'))
            plt.close()
        else:
            plt.show()
            
        # Creating color grid for bin means
        cmap = matplotlib.cm.gray
        plt.pcolormesh(fluxRatioBin_means,cmap=cmap)
        plt.axis([0,m,0,n])
        plt.title('Object-wise Flux Ratio Bin Means')
        plt.colorbar()
        plt.xlabel('X Bin')
        plt.ylabel('Y Bin')
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'fluxRatioBin_means.png'))
            plt.close()
        else:
            plt.show()