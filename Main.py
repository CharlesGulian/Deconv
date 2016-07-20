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
import sex_config
import fits_tools

# Image deconvolution project:
# Main script for data analysis, image comparison, photometric statistics, and more

# Good image comparison
goodImage1 = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
goodImage2 = 'AstroImages/Good/fpC-7006-x5226-y115_stitched_alignCropped.fits'
goodImage3 = 'AstroImages/Good/fpC-4868-x4211-y138_stitched_alignCropped.fits'
goodImage4 = 'AstroImages/Good/fpC-6383-x5176-y121_stitched_alignCropped.fits'

goodImgs = [goodImage1,goodImage2,goodImage3,goodImage4]
goodImgs = [goodImage1,goodImage3]

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
        
        # Printing image file names
        print testImage1+' '+testImage2
        
        medSub = False
        if medSub:
            testImage1_medSub = testImage1.replace('.fits','_medSub.fits')
            testImage2_medSub = testImage2.replace('.fits','_medSub.fits')
            
            fits_tools.subtractMedian(testImage1,new_image_file=testImage1_medSub)
            fits_tools.subtractMedian(testImage2,new_image_file=testImage2_medSub)
            
            testImage1 = testImage1_medSub
            testImage2 = testImage2_medSub
        
        # Write configuration files for SExtractor comparison
        
        # Writing configuration file for first image
        fig1 = sex_config.configure(testImage1+','+testImage1,'default.sex','default.param',dual=True)
        # Do default configuration        
        fig1.default_config()
        
        # Writing onfiguration file for second image
        fig2 = sex_config.configure(testImage1+','+testImage2,'default.sex','default.param',dual=True)
        # Do default configuration
        fig2.default_config()
        
        varParam = 'Filter'
        PhotAper = np.linspace(5.0,25.0,5)
        BackSize = np.linspace(150,350,5)
        Filter = ['Y','N']
        varParamRange = Filter
        for k in range(len(varParamRange)):
            # Adjust Kron factor
            fig1.reconfigure('FILTER',Filter[k])
            # Change check image name
            #temp = fig1.config_dict['CHECKIMAGE_NAME'].replace('Results/CheckImages','Figures/Jul19/'+varParam).replace('.fits','_{}{}.fits'.format(varParam,varParamRange[k]))
            #fig1.reconfigure('CHECKIMAGE_NAME',temp)
            # Write new configuration file for first image
            fig1.write_config_file(new_config_file='copy_compare1.sex',new_param_file='copy_compare1.param')
            
            # Adjust Kron factor
            fig2.reconfigure('FILTER',Filter[k])
            # Change check image name
            #temp = fig2.config_dict['CHECKIMAGE_NAME'].replace('Results/CheckImages','Figures/Jul19/'+varParam).replace('.fits','_{}{}.fits'.format(varParam,varParamRange[k]))
            #fig2.reconfigure('CHECKIMAGE_NAME',temp)
            # Write new configuration file for second image
            fig2.write_config_file(new_config_file='copy_compare2.sex',new_param_file='copy_compare2.param')
            
            # Compare images in SExtractor
            output = pysex.compare(testImage1,testImage2,'copy_compare1.sex','copy_compare2.sex')
            
            # Get image tags
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
                
            # Create sex_stats.data objects:
            img1data = sex_stats.data(outputCat1)
            img2data = sex_stats.data(outputCat2)
            
            # Create .reg files from output catalogs
            CREATE_regFiles = False
            if CREATE_regFiles:
                img1data.create_regFile()
                img2data.create_regFile()
        
        
        
            #-----------------------------------------------------------------------------#
            # Flux ratio analysis:
            
            flux1,flux2 = img1data.get_data('FLUX_APER'),img2data.get_data('FLUX_APER')
            #mag1,mag2 = img1data.get_data('MAG_AUTO'),img2data.get_data('MAG_AUTO')
            #flux1,flux2 = mag1,mag2
            x,y = img1data.get_data('X_IMAGE'),img1data.get_data('Y_IMAGE')
            
            flux1,flux2 = np.array(flux1),np.array(flux2)
            
            #'''
            fluxAvg = 0.5*(flux1+flux2)
            fluxRatio = np.divide(flux1,flux2)
            fluxRatio_mean = np.mean(fluxRatio)
            fluxRatio_std = np.std(fluxRatio)
            #fluxRatio_meanSubtracted = fluxRatio - fluxRatio_mean # (NOT MEAN SUBTRACTED)
            #'''
            
            #''' 
            print ''
            print '{} = {}'.format(varParam,varParamRange[k])
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
            print 'Minimum values of images: ', np.min(fits_tools.getPixels(testImage1)),' ',np.min(fits_tools.getPixels(testImage2))
            print 'Median values of images: ', np.median(fits_tools.getPixels(testImage1)),' ',np.median(fits_tools.getPixels(testImage2))
            print 'Mean values of images: ', np.mean(fits_tools.getPixels(testImage1)),' ',np.mean(fits_tools.getPixels(testImage2))
            print ' '
            #'''
            
            # Creating histogram of flux1 and flux2
            # plt.hist has range=(tup1,tup2) to set upper and lower bounds of histogram
            plt.hist(flux1,bins=70,normed=False,color='green')
            plt.hist(flux2,bins=70,normed=False,color='blue')
            plt.title('Histogram of Object Flux for Entire Image')
            plt.ylabel('Frequency (N)')
            plt.xlabel('Object flux')
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(curr_dir,'Figures','Jul19',varParam,'flux1_flux2_hist_{}{}.png'.format(varParam,varParamRange[k])))
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
                plt.savefig(os.path.join(curr_dir,'Figures','Jul19',varParam,'fluxRatio_hist_{}{}.png'.format(varParam,varParamRange[k])))
                plt.close()
            else:
                plt.show()
            
            # Creating color plot of object-wise flux ratio
            cmap = matplotlib.cm.jet
            plt.scatter(x,y,s=25.0*img1data.get_data('A_IMAGE'),c=fluxRatio-np.median(fluxRatio),marker='o',vmin=np.min(fluxRatio),vmax=np.max(fluxRatio),alpha=0.85)
            plt.axis([0,1600,0,1600])  
            plt.colorbar()
            plt.title('Map of Object-wise Flux Ratio')
            plt.xlabel('X_IMAGE')
            plt.ylabel('Y_IMAGE')
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(curr_dir,'Figures','Jul19',varParam,'fluxRatio_map_{}{}.png'.format(varParam,varParamRange[k])))
                plt.close()
            else:
                plt.show()
            # Creating histogram of object-wise flux ratio for 4x4 bins
            m,n = 4,4
            xBins,yBins,fluxRatioBins = sex_stats.binData(x,y,fluxRatio,M=m,N=n)
            for i in range(m):
                for j in range(n):
                    plt.subplot(m,n,(n*i + (j+1)))
                    #plt.hist(fluxRatioBins[i,j],bins=int(abs(np.max(fluxRatio)-np.min(fluxRatio))*10.0))
                    plt.hist(fluxRatioBins[i,j],bins=70)
                    plt.axis([np.min(fluxRatioBins[i,j]),np.max(fluxRatioBins[i,j]),0.0,15.0])
            
            SAVE = True
            if SAVE:
                plt.savefig(os.path.join(curr_dir,'Figures','Jul19',varParam,'fluxRatioBins_hist_{}{}.png'.format(varParam,varParamRange[k])))
                plt.close()
            else:
                plt.show()
            
        '''
        fluxRatioBin_Avgs = np.zeros([m,n])
        emptyBins = []
        for i in range(m):
            for j in range(n):
                # Clipping data in bins:
                fluxRatioBins_sigmaClipped = []
                fluxRatioBins_excess = []
                for k in range(len(fluxRatioBins[i,j])):
                    if np.abs((fluxRatioBins[i,j])[k]) <= maxSig[s]*np.std(fluxRatioBins[i,j]):
                        fluxRatioBins_sigmaClipped.append(fluxRatioBins[i,j][k])
                    else:
                        fluxRatioBins_excess.append()
                if len(fluxRatioBins_sigmaClipped) == 0:
                    emptyBins.append('{},{}'.format(str(i),str(j)))
                fluxRatioBins[i,j] = fluxRatioBins_sigmaClipped
                fluxRatioBin_Avgs[i,j] = np.mean(fluxRatioBins_sigmaClipped)
                    
        # Masking NaNs in fluxRatioBin_Avgs:

        fluxRatioBin_Avgs_Masked = np.ma.array(fluxRatioBin_Avgs,mask=np.isnan(fluxRatioBin_Avgs))
        cmap = matplotlib.cm.gray
        cmap.set_bad('r',1.)
        
        #print np.nanmean(fluxRatioBin_Avgs)-2.0,' ',np.nanmean(fluxRatioBin_Avgs)+2.0
        plt.pcolormesh(fluxRatioBin_Avgs_Masked,cmap=cmap,vmin=np.nanmean(fluxRatioBin_Avgs)-2.0,vmax=np.nanmean(fluxRatioBin_Avgs)+2.0)
        plt.colorbar()
        plt.xlabel('X Bin')
        plt.ylabel('Y Bin')
        plt.title('Flux Ratio Bin Averages: {} x {}'.format(m,n))
        if not os.path.exists(os.path.join(curr_dir,'Figures','Jul14Imgs','ObjBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10]))):
            os.mkdir(os.path.join(curr_dir,'Figures','Jul14Imgs','ObjBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10])))
        plt.savefig(os.path.join(curr_dir,'Figures','Jul14Imgs','ObjBin','{}_{}'.format(img_tag1[0:10],img_tag2[0:10]),'fluxRatioBin_Avgs_sigmaClip{}.png'.format(str(maxSig[s])[0:4])))
        plt.close()
        
        plot = False # Warning: do not change to true unless length of maxSig small
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
        #'''
        break
    break


    

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