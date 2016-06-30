# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:36:10 2016

@author: charlesgulian
"""
import os
import matplotlib.pyplot as plt
import numpy as np
#from astropy.io import fits

dir_name = os.getcwd()

n = open(dir_name + '/img_name_compare.txt','r')
img_tag1 = (n.readline()).strip()
img_tag2 = (n.readline()).strip()

    
#filename = dir_name+'/AstroImages/'img_tag+'.fits'
#hdulist = fits.open(filename)  # open a FITS file
#fits_hdr = hdulist[0].header
#hdulist.close()

f1 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag1+'_compare.cat','r')
lines = f1.readlines()
cat_hdr_size = len(lines[len(lines)-1].split()) # Catalog header size; number of output columns in default.param
f1.seek(0) # Reset pointer to start of file
lines = None

params = []
PRINT = True
if PRINT:
    for i in range(cat_hdr_size):
        param = (f1.readline()).split()[2] # Get parameter name
        print param
        params.append(param)
else:
    for i in range(cat_hdr_size):
        param = (f1.readline()).split()[2] # Get parameter name
        params.append(param)

data1 = []
for line in f1:
    line = line.strip()
    column = line.split()
    source = {}
    for i in range(len(params)):
        source[params[i]] = float(column[i])
    data1.append(source)
    
f2 = open(dir_name+'/Results/'+img_tag1+'_'+img_tag2+'_compare.cat','r')
data2 = []
for i in range(cat_hdr_size):
    f2.readline()
for line in f2:
    line = line.strip()
    column = line.split()
    source = {}
    for i in range(len(params)):
        source[params[i]] = float(column[i])
    data2.append(source)

if len(data1) != len(data2):
    print "Warning: data arrays are of different lengths"

flux1,flux2 = [],[]
x,y = [],[]
for i in range(len(data1)):
    Obj1,Obj2 = data1[i],data2[i]
    flux1.append(Obj1['FLUX_BEST'])
    flux2.append(Obj2['FLUX_BEST'])
    x.append(Obj1['X_IMAGE'])
    y.append(Obj1['Y_IMAGE'])

x,y = np.array(x),np.array(y)

flux1,flux2 = np.array(flux1),np.array(flux2)
fluxAvg = 0.5*(flux1+flux2)
fluxRatio = np.divide(flux1,flux2)
fluxRatio_mean = np.mean(fluxRatio)
fluxRatio_std = np.std(fluxRatio)
fluxRatio_meanSubtracted = fluxRatio - fluxRatio_mean

maxSig = np.linspace(0.05,1.0,10) # Sigma



for j in range(len(maxSig)): 
    
    #''' # Exclude clipped data points, copy to fluxRatio_excess
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
        
    ''' # Include clipped data points, set equal maxSig[j]
    fluxRatio_meanSubtracted_sigmaClipped = fluxRatio_meanSubtracted
    for i in range(len(fluxRatio)):
            if np.abs(fluxRatio_meanSubtracted_sigmaClipped[i]) >= maxSig[j]*fluxRatio_std: # maxSig[j] = 1.0 later
                fluxRatio_meanSubtracted_sigmaClipped[i] = maxSig[j]*fluxRatio_std*np.sign(fluxRatio_meanSubtracted_sigmaClipped[i])
            else:
                next
    #'''
    
    plot = True
    if plot:
        # Plotting source-wise flux ratio
        #plt.scatter(x, y, s=25*np.log10(0.1*fluxAvg), c=fluxRatio_meanSubtracted_sigmaClipped, vmin=-1.5*maxSig[j], vmax=1.5*maxSig[j], alpha=0.75)
        plt.scatter(x_clip, y_clip, s=25*np.log10(0.1*np.array(fluxAvg_clip)), c=fluxRatio_meanSubtracted_sigmaClipped, vmin=-1.5*maxSig[j]*fluxRatio_std, vmax=1.5*maxSig[j]*fluxRatio_std, alpha=0.75)        
        plt.axis([0,1600,0,1600])
        plt.colorbar()
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')
        plt.title('Flux Ratio Color Map: sigma cutoff = '+str(maxSig[j])[0:4])
        plt.savefig((dir_name+'/Figures/{}_{}_maxSig{}_fluxRatio.png'.format(img_tag1, img_tag2, str(maxSig[j])[0:4])))
        plt.close()
    
    hist = False
    if hist:
        #plt.subplot(211)
        #plt.hist(fluxRatio,bins=20,color='green')
        #plt.title('Histogram of Flux Ratio')
        #plt.subplot(212)
        plt.hist(fluxRatio_meanSubtracted_sigmaClipped,bins=20,color='green')
        plt.title('Histogram of Flux Ratio')
        plt.ylabel(('Mean subtracted + clipped @ {} sigma').format(str(maxSig[j])[0:4]))
        plt.xlabel('Flux ratio')
        plt.savefig((dir_name+'/Figures/Hist_{}_{}_maxSig{}_fluxRatio.png'.format(img_tag1, img_tag2, str(maxSig[j])[0:4])))
        plt.close()
