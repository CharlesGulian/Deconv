# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:34:39 2016

@author: charlesgulian
"""
import os
curr_dir = os.getcwd()
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from astropy.io import fits

import pysex
import sex_stats
import sex_config
import fits_tools
from SciPyKDE import kde_scipy


# Image deconvolution project:
# Main script for data analysis, image comparison, photometric statistics, and more

class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID

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
comparisonImages = [coaddedImage,deconvolvedImage]

# ===============================================================================
# Image comparisons:
for img1 in comparisonImages:
    for img2 in comparisonImages:
        if img1 == img2:
            continue
        if img1.category == 'Deconvolved':
            continue
        
        print '\nDetection image: {}\nMeasurement image: {}\n'.format(img1.filename,img2.filename)
        
        # ===============================================================================
        # Make directory to save results and figures to different comparison categories (Good vs. Good, Good vs. Deconvolved, etc.)
        #image1,image2 = img1.filename,img2.filename # Get image filenames
        new_dir = os.path.join(curr_dir,'Figures','Oct12',img1.category+img1.ID+'_vs_'+img2.category+img2.ID)
        if not os.path.exists(new_dir):
            os.mkdir(os.path.join(new_dir))
            
        # ===============================================================================
        # Adjust image data for best comparison with SExtractor    
            
        imgs = [img1,img2]
        for i,img in enumerate(imgs):
            
            if img.category == 'Good':
                temp = img.filename.replace('.fits','_biasSub.fits')
                fits_tools.subtractBias(img.filename,new_image_file=temp,bias=1000.0)
                img.filename = temp
            
            if img.category == 'Coadded':
                temp = img.filename.replace('.fits','_masked.fits')
                fits_tools.maskImage(img.filename,maskFile,masked_image_file=temp)
                img.filename = temp
                
            if img.category == 'Deconvolved':
                temp = img.filename.replace('.fits','_masked.fits')
                fits_tools.maskImage(img.filename,maskFile,masked_image_file=temp)
                img.filename = temp
              
        # Re-align images
        header1 = fits.getheader(img1.filename)
        header2 = fits.getheader(img2.filename)
        
        x_offset1,y_offset1 = header1['x_offset'],header1['y_offset']
        x_offset2,y_offset2 = header2['x_offset'],header2['y_offset']
        
        del_x = x_offset2 - x_offset1
        del_y = y_offset2 - y_offset1
        
        # Shift second image to align with first image        
        shifted_image2_data = fits_tools.shift_image(fits.getdata(img2.filename),del_x,del_y)
        
        # New image file names:
        new_image1 = os.path.join(new_dir,img1.category+img1.ID+'_copy'+'.fits')
        new_image2 = os.path.join(new_dir,img2.category+img2.ID+'_shifted'+'.fits')
        
        fits.writeto(new_image1,fits.getdata(img1.filename),header=fits.getheader(img1.filename),clobber=True)
        fits.writeto(new_image2,shifted_image2_data,header=fits.getheader(img2.filename),clobber=True)
        
        img1.filename,img2.filename = new_image1,new_image2
        
        # =======================================================================
        # Write configuration files for SExtractor comparison
        
        # Create configuration object for first image
        fig1 = sex_config.configure(img1.filename+','+img1.filename,'default.sex','default.param',dual=True)
        # Do default configuration
        fig1.default_config()
        
        # Create configuration object for second image
        fig2 = sex_config.configure(img1.filename+','+img2.filename,'default.sex','default.param',dual=True)
        # Do default configuration
        fig2.default_config()
        
        # Reconfigure detection threshold parameters
        imgs = [img1,img2]
        figs = [fig1,fig2]
        for i,img in enumerate(imgs):
            if img.category == 'Deconvolved':
                figs[i].reconfigure('THRESH_TYPE','ABSOLUTE')
                figs[i].reconfigure('DETECT_THRESH',15.0)
                figs[i].reconfigure('BACK_TYPE','MANUAL')
                figs[i].reconfigure('BACK_VALUE',0.0)
                #print '\nSetting detection threshold = {0}'.format(figs[i].config_dict['DETECT_THRESH'])
                #print 'Turning off background estimation/subtraction\n'
            
            if img.category == 'Coadded':
                print '\nSubtracting median from co-added image: median = {0}\n'.format(np.median(fits_tools.getPixels(img.filename)))
                temp = (img.filename).replace('.fits','_medSub.fits')
                fits_tools.subtractMedian(img.filename,new_image_file=temp)
                img = Image(temp,'Coadded','')
                
                print '\nSetting DETECT_THRESH = 1.0 for co-added image\n'
                figs[i].reconfigure('THRESH_TYPE','ABSOLUTE')
                figs[i].reconfigure('DETECT_THRESH',15.0)
                figs[i].reconfigure('BACK_TYPE','MANUAL')
                figs[i].reconfigure('BACK_VALUE',0.0)
        
        # Write configuration files
        fig1.write_config_file(new_config_file='copy_compare1.sex',new_param_file='copy_compare1.param')
        fig2.write_config_file(new_config_file='copy_compare2.sex',new_param_file='copy_compare2.param')
        
        # =======================================================================
        # Compare images in SExtractor
        
        pysex.compare(img1.filename,img2.filename,'copy_compare1.sex','copy_compare2.sex')
        
        # =======================================================================
        # Retrieve data from output catalogs        
        
        # Get image tags
        img_tag1 = (os.path.split(img1.filename)[1])
        img_tag1 = img_tag1[0:len(img_tag1)-len('.fits')]
        img_tag2 = (os.path.split(img2.filename)[1])
        img_tag2 = img_tag2[0:len(img_tag2)-len('.fits')]
        
        outputCat1 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag1+'_compare.cat')
        print outputCat1
        if not os.path.exists(outputCat1):
            print 'Error: first output catalog path does not exist'
        
        outputCat2 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag2+'_compare.cat')
        if not os.path.exists(outputCat2):
            print 'Error: second output catalog path does not exist'
        
        # Create sex_stats.data objects:
        img1data = sex_stats.data(outputCat1)
        img2data = sex_stats.data(outputCat2)
        
        # Printing background detection threshold:
        threshold = img1data.get_data('THRESHOLD')
        print 'Background detection threshold: ', threshold
        
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
        
        imageData1,imageData2 = fits_tools.getPixels(img1.filename),fits_tools.getPixels(img2.filename)
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
            flux1[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData1) # When computing, use original (non-masked) images
            flux2[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData2)

        # Compute object-wise flux ratio
        fluxRatio = np.divide(flux1,flux2) # Divide flux 1 by flux 2
        log_fluxRatio = np.log(fluxRatio)        
        fluxRatio = log_fluxRatio # Note that from here on, "flux ratio" = ln(f1/f2)
        
        # =======================================================================
        """
        # Looking at flux ratio outliers in deconvolved vs. co-added
        outliers = {}
        outliers['high'],outliers['low'] = {},{}
        
        outliers['high']['inds'],outliers['high']['fluxRatio'],outliers['high']['x'],outliers['high']['y'],outliers['high']['flux1'],outliers['high']['flux2'] = [],[],[],[],[],[]
        outliers['high']['aperture_radius'] = []   
        outliers['low']['inds'],outliers['low']['fluxRatio'],outliers['low']['x'],outliers['low']['y'],outliers['low']['flux1'],outliers['low']['flux2'] = [],[],[],[],[],[]    
        outliers['low']['aperture_radius'] = []
              
        for i in range(len(fluxRatio)):
            if fluxRatio[i] >= np.median(fluxRatio): # Note that "high outliers" include any value above the median flux ratio value
                outliers['high']['inds'].append(i)
                outliers['high']['fluxRatio'].append(fluxRatio[i])
                outliers['high']['flux1'].append(flux1[i])
                outliers['high']['flux2'].append(flux2[i])
                outliers['high']['aperture_radius'].append(rads[i])
                outliers['high']['x'].append(x[i])
                outliers['high']['y'].append(y[i])
            else:
                outliers['low']['inds'].append(i)
                outliers['low']['fluxRatio'].append(fluxRatio[i])
                outliers['low']['flux1'].append(flux1[i])
                outliers['low']['flux2'].append(flux2[i])
                outliers['low']['aperture_radius'].append(rads[i])
                outliers['low']['x'].append(x[i])
                outliers['low']['y'].append(y[i])
        
        # Sort and find top outliers, bottom outliers
        temp,sorted_high_outlier_indices = zip(*sorted(zip(outliers['high']['fluxRatio'],range(len(outliers['high']['fluxRatio'])))))
        top_outlier_inds = sorted_high_outlier_indices[::-1][0:10]
        temp,sorted_low_outlier_indices = zip(*sorted(zip(outliers['low']['fluxRatio'],range(len(outliers['low']['fluxRatio'])))))
        bottom_outlier_inds = sorted_low_outlier_indices[0:10]
        
        # Create .reg files for top and bottom outliers:
        sex_stats.createRegFile_outliers(outliers['high'],top_outlier_inds,os.path.join(new_dir,'top_outliers.reg'))
        sex_stats.createRegFile_outliers(outliers['low'],bottom_outlier_inds,os.path.join(new_dir,'bottom_outliers.reg'))
        
        print 'High outliers: flux ratio | x | y'
        for i in top_outlier_inds:
            
            # 2) Remove outliers (???) (try np.delete(fluxRatio,fluxRatio[outliers['ind'][i]]) )? or else delete directly from the data dictionary/object
            # 3) Look @ top and bottom 5 outliers in deconvolved and co-added images; put .reg file in deconvolved vs. coadded folder
                # --> Use create_regfile type function for x,y,radius (circular)
            # 4) Commit to github repository?
            print outliers['high']['fluxRatio'][i],outliers['high']['x'][i],outliers['high']['y'][i]
            print ''
            #print outliers['high']['flux1'][i],outliers['high']['flux2'][i]
            #print outliers['high']['aperture_radius'][i]
        
        print 'Low outliers: flux | x | y'
        for i in bottom_outlier_inds:
            print outliers['low']['fluxRatio'][i],outliers['low']['x'][i],outliers['low']['y'][i]
            #print outliers['low']['flux1'][i],outliers['low']['flux2'][i]
            print ''
         
        '''
        # Delete high outliers
        np.delete(x,top_outlier_inds)
        np.delete(y,top_outlier_inds)
        np.delete(rads,top_outlier_inds)
        np.delete(fluxRatio,top_outlier_inds)
        #'''
        #"""    
        # =======================================================================
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
        print 'Median flux ratio: ', np.median(fluxRatio)
        print 'Mean flux ratio: ', np.mean(fluxRatio)
        print 'Standard deviation of flux ratio: ', np.std(fluxRatio)
        print ''
        print 'Minimum values of images: ', np.min(imageData1),' ',np.min(imageData2)
        print 'Median values of images: ', np.median(imageData1),' ',np.median(imageData2)
        print 'Mean values of images: ', np.mean(imageData1),' ',np.mean(imageData2)
        print ' '
        #'''
                
        # =======================================================================
        # Splitting up data into N different brightness regimes
                
        N = 4 # Number of brightness regimes (each corresponding to adjacent partitions of width b*std(flux1) of the flux1 distribution; final partition is all points greater than (N-1)*std(flux1))
        
        '''b = 0.05 # Controls width of partitions
        regimes = {}
        for i in range(N-1):
            regimes[i+1] = []
            inds1 = np.where(flux1 <= (i+1)*b*np.std(flux1))[0]
            inds2 = np.where(flux1 > (i)*b*np.std(flux1))[0]
            inds = np.intersect1d(inds1,inds2)
            regimes[i+1].extend(inds)
        regimes[N] = []
        inds = np.where(flux1 > (N-1)*b*np.std(flux1))[0]
        regimes[N].extend(inds)
        #'''
        
        inds_unsorted = np.arange(len(flux1))
        flux1_sorted,inds_sorted = zip(*sorted(zip(flux1,inds_unsorted)))
        L = len(flux1)
        bin_size = L/N
    
        regimes = {}
        regimes[0] = []
        inds = np.arange(len(flux1))
        regimes[0].extend(inds) # Oth regime corresponds to original (full) dataset
        for i in range(N-1):
            regimes[i+1] = []
            inds = np.array(inds_sorted[i*bin_size:(i+1)*bin_size-1])
            regimes[i+1].extend(inds)
        regimes[N] = []
        inds = np.array(inds_sorted[(N-1)*bin_size::])
        regimes[N].extend(inds)            
        
        # =======================================================================
        # Plotting
        
        # Creating histogram of flux1 and flux2
        # plt.hist has range=(tup1,tup2) to set upper and lower bounds of histogram
        histRange = (np.min(flux1),np.mean(flux1)+1.5*np.std(flux1))    
        plt.hist(flux1,bins=70,normed=False,range=histRange,color='green',label=img1.category)
        plt.hist(flux2,bins=70,normed=False,range=histRange,color='blue',label=img2.category)
        plt.legend()
        plt.title('Histogram of Object Flux for Entire Image')
        plt.ylabel('Frequency (N)')
        plt.xlabel('Object flux')
        SAVE = True
        if SAVE:
            plt.savefig(os.path.join(new_dir,'flux1_flux2_hist.png'))
            plt.close()
        else:
            plt.close()
            
        '''   
         # Creating histogram of FLUX_RADIUS for PHOT_FLUXFRAC = 0.75,0.8,0.9
        plt.hist(rads,bins=70,color='green')
        plt.title('Histogram of 9/10-light radii')
        plt.xlabel('9/10-light radius')
        plt.ylabel('Frequency (N)')
        SAVE = False
        if SAVE:
            plt.savefig(os.path.join(new_dir,'hist_flux_radius0.90.png'))
            plt.close()
        else:
            plt.close()
        #'''
        
        # =======================================================================
        # Plotting for different brightness regimes
        fluxRatio_og,x_og,y_og = fluxRatio,x,y
        for i in range(N+1):
            
            fluxRatio = fluxRatio_og[regimes[i]]
            x,y = x_og[regimes[i]],y_og[regimes[i]]
            
            print '\nNumber of objects in brightness regime: {}'.format(str(len(regimes[i])))
            print 'Median flux value of brightness regime: {}\n'.format(str(np.median(flux1[regimes[i]])))
            
            # Creating histogram of flux1/flux2 (object-wise flux ratio)
            #plt.hist(fluxRatio,bins=70,range=(-25.0,25.0),color='green') # Range = (-25.0,25.0)
            #plt.hist(fluxRatio,bins=int(abs(np.max(fluxRatio)-np.min(fluxRatio))*10.0),range=(np.min(fluxRatio),np.max(fluxRatio)),color='green') 
            histRange = (np.mean(fluxRatio)-2.5*np.std(fluxRatio),np.mean(fluxRatio)+2.5*np.std(fluxRatio))      
            plt.hist(fluxRatio,bins=70,range=histRange,color='green')           
            plt.title('Histogram of Object-wise Flux Ratio')
            plt.ylabel('Frequency (N)')
            plt.xlabel('Object-wise flux ratio')
            # Setting histogram axis scale
            hist_axis = [np.max([np.median(fluxRatio_og)-4.5*np.std(fluxRatio_og),np.min(fluxRatio_og)]),np.median(fluxRatio_og)+4.5*np.std(fluxRatio_og),-0.1,np.max(plt.hist(fluxRatio_og,bins=70,range=histRange,color='green')[0]+1.0)]
            plt.axis(hist_axis)
            print '\nHistogram x-axis range: {0},{1}\n'.format(str(hist_axis[0]),str(hist_axis[1]))
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatio_hist_{}.png'.format('regime'+str(i))))
                plt.close()
            else:
                plt.close()
            
            # Creating color plot of object-wise flux ratio
            cmap = cm.get_cmap('rainbow')
            #from scipy import stats
            colors = fluxRatio # - stats.mode(fluxRatio_og)[0][0] # How to find mode of continuous distribution? Which x-bin has highest count in distribution?
            vmin = np.max([np.median(fluxRatio_og)-2.5*np.std(fluxRatio_og),np.min(fluxRatio_og)])
            vmax = np.median(fluxRatio_og)+1.5*np.std(fluxRatio_og)
            plt.scatter(x,y,s=35.0,c=colors,marker='o',cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.85,linewidths=0.0)
            plt.axis([0,1600,0,1600])
            plt.colorbar()
            plt.title('Map of Object-wise Flux Ratio')
            plt.xlabel('X_IMAGE')
            plt.ylabel('Y_IMAGE')
            plt.figtext(.65,.85,'N = '+str(len(fluxRatio)))
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatio_colorplot_{}.png'.format('regime'+str(i))))
                plt.close()
            else:
                plt.close()
        
            #'''
            # Creating histogram of object-wise flux ratio for 4x4 bins
            m,n = 4,4
            xBins,yBins,fluxRatioBins = sex_stats.binData(x,y,fluxRatio,M=m,N=n)
            binSizes = []
            for k in fluxRatioBins.iteritems():
                binSizes.append(len(k[1]))
            maxBinSize = max(binSizes)
            histRange = (np.mean(fluxRatio)-2.5*np.std(fluxRatio),np.mean(fluxRatio)+2.5*np.std(fluxRatio))
            fluxRatioBin_means = np.zeros([m,n])
            for i in range(m):
                for j in range(n):
                    plt.subplot(m,n,(n*i + (j+1)))
                    bins = int(np.ceil((np.max(fluxRatioBins[i,j]) - np.min(fluxRatioBins[i,j]))*4))
                    plt.hist(fluxRatioBins[i,j],bins=bins)
                    fluxRatioBin_means[i,j] = np.mean(fluxRatioBins[i,j])
                    plt.axis([np.mean(fluxRatio)-4.5*np.std(fluxRatio),np.mean(fluxRatio)+4.5*np.std(fluxRatio),0,maxBinSize+1])
            SAVE = False
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatioBins_hist_{}.png'.format('regime'+str(i))))
                plt.close()
            else:
                plt.close()
        
            # Creating KDE of flux ratio
            from kde_danielsmith import kde as kde_ds        
            #fr_grid = np.linspace(np.min(fluxRatio)-0.5,np.max(fluxRatio)+0.5,1000)
            #fluxRatio_KDE = kde_scipy(fluxRatio,fr_grid,0.01) 
            bandwidth, fr_grid, fluxRatio_KDE = kde_ds(fluxRatio)
            #print 'Optimal bandwidth for KDE of flux ratio: ',bandwidth
            
            plt.plot(fr_grid,fluxRatio_KDE,'k',linewidth=2.5)
            histRange = (np.mean(fluxRatio_og)-2.5*np.std(fluxRatio_og),np.mean(fluxRatio_og)+2.5*np.std(fluxRatio_og))
            #plt.hist(fluxRatio,bins=140,range=histRange,color='r')           
            plt.title('KDE of Object-wise Flux Ratio')
            plt.ylabel('Frequency (N)')
            plt.xlabel('Object-wise flux ratio')
            kde_axis = [np.max([np.median(fluxRatio_og)-2.5*np.std(fluxRatio_og),np.min(fluxRatio_og)]),np.median(fluxRatio_og)+2.5*np.std(fluxRatio_og),-0.1,np.max(fluxRatio_KDE)+1.0]
            plt.axis(kde_axis)
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatio_KDE_{}.png'.format('regime'+str(i))))
                plt.close()
            else:
                plt.show()
        
            # Creating KDE of 4x4 flux ratio bins
            for i in range(m):
                for j in range(n):
                    plt.subplot(m,n,(n*i + (j+1)))
                    fluxRatioBin_KDE = kde_scipy(fluxRatioBins[i,j],fr_grid,16*bandwidth)
                    #bandwidth, fr_grid, fluxRatioBin_KDE = kde_ds(fluxRatioBins[i,j],guess=1.0)
                    plt.plot(fr_grid,fluxRatioBin_KDE,'k',linewidth=2.0)
                    plt.grid()
                    #plt.gca().get_xaxis().set_visible(False)
                    plt.gca().get_yaxis().set_visible(False)
                    fluxRatioBin_means[i,j] = np.mean(fluxRatioBins[i,j])
                    plt.axis(kde_axis)
                    #  plt.axes(frameon=False)
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatioBins_KDE_{}.png'.format('regime'+str(i))))
                plt.close()
            else:
                plt.show()
         
            # Creating color grid for bin means
            cmap = cm.get_cmap('gray_r')
            plt.pcolormesh(fluxRatioBin_means,cmap=cmap)
            plt.axis([0,m,0,n])
            plt.colorbar()
            plt.title('Object-wise Flux Ratio Bin Means')
            #plt.colorbar()
            plt.xlabel('X Bin')
            plt.ylabel('Y Bin')
            SAVE = False
            if SAVE:
                plt.savefig(os.path.join(new_dir,'fluxRatioBin_means_{}.png').format('regime'+str(i)))
                plt.close()
            else:
                plt.close()