# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 20:31:36 2016
@author: charlesgulian
"""
# Library of statistics functions for photometric statistics and data analysis

import os
os.chdir('/Users/annepstein/Work/Deconv')
curr_dir = os.getcwd()
import numpy as np
import pandas as pd
from astropy.io import fits

class data:
    
    def __init__(self,outputCatName):
        
        self.outputCatName = outputCatName
        with open(outputCatName) as f:
            f = open(outputCatName)
            lines = f.readlines()
            #print lines
            cat_hdr_size = len(lines[len(lines)-1].split()) # Catalog header size; number of output columns in default.param
            f.seek(0) # Reset pointer to start of file
            lines = None
            
            params = []
            PRINT = False
            if PRINT:
                for i in range(cat_hdr_size):
                    param = (f.readline()).split()[2] # Get parameter name
                    print param
                    params.append(param)
            else:
                for i in range(cat_hdr_size):
                    param = (f.readline()).split()[2] # Get parameter name
                    params.append(param)
            self.Params = params
            
            _Data = []
            for line in f:
                line = line.strip()
                column = line.split()
                source = {}
                for i in range(len(params)):
                    source[params[i]] = float(column[i])
                _Data.append(source)
            self.Data = _Data
    
    def get_data(self,paramName):
        
        paramData = []
        for i in range(len(self.Data)):
            Obj = self.Data[i] # Get data for i'th object in catalog
            paramData.append(Obj[paramName])
        return np.array(paramData)
    
    def create_regFile(self,regFileName=None):
        # Creates a DS9 .reg file for SExtractor output catalog using appropriate columns
    
        if regFileName == None:
            regFileName = self.outputCatName.replace('.cat','.reg').replace('Results/','Results/regFiles/')
        
        with open(regFileName,'w') as g:
            # Initializing .reg file
            g.write('# Region file format: DS9 version 4.1\n')
            g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            g.write('image\n')
            for i in range(len(self.Data)):
                Obj = self.Data[i]
                Line = 'ellipse('+str(Obj['X_IMAGE'])+','+str(Obj['Y_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['A_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['B_IMAGE'])+','+str(Obj['THETA_IMAGE'])+')\n'
                g.write(Line)

def binData(x1,x2,y,M=3,N=3,ImgD1=[0,1600],ImgD2=[0,1600]):
    '''
    - Bins y-data along x1 and x2 into MxN bins (default MxN = 3x3)
    - Bin sizes/intervals are measured from image dimensions in pixels, NOT from
        min/max values of x1 and x2
        * For this option, specify ImgD1=[min(x1),max(x1)], ImgD2=[min(x2),max(x2)]
    - Default image dimension = 1600x1600 pixels
    '''
    
    if len(x1) != len(x2):
        raise Exception('Error: x1 and x2 must have same size')
        return
    
    x1Dict,x2Dict,yDict = {},{},{}
    for i in range(M):
        for j in range(N):
            x1Dict[i,j] = []
            x2Dict[i,j] = []
            yDict[i,j] = []
    
    x1bin_length= float(ImgD1[1]-ImgD1[0])/float(M)
    x2bin_length = float(ImgD2[1]-ImgD2[0])/float(N)
    
    for i in range(M):
        x1_inds_lower = np.where(x1 >= i*x1bin_length)
        x1_inds_upper = np.where(x1 < (i+1)*x1bin_length)
        x1_inds = np.intersect1d(x1_inds_lower,x1_inds_upper)
        for j in range(N):
            x2_inds_lower = np.where(x2 >= j*x2bin_length)
            x2_inds_upper = np.where(x2 < (j+1)*x2bin_length)
            x2_inds = np.intersect1d(x2_inds_lower,x2_inds_upper)
            
            inds = np.intersect1d(x1_inds,x2_inds) # These indices correspond to points in the ith,jth bin
            
            x1Dict[i,j].extend(x1[inds])
            x2Dict[i,j].extend(x2[inds])
            yDict[i,j].extend(y[inds])
            yDict[i,j] = np.array(yDict[i,j])
            
    return x1Dict,x2Dict,yDict


def binImage(pixelArray,M=3,N=3):
    # Bins pixels along image axes into MxN bins (default MxN = 3x3)  
    
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
    
    
def createRegFile(outputCatName,regFileName=None):
    if regFileName == None:
        regFileName = outputCatName.replace('.cat','.reg')
    
    temp = data(outputCatName)
    outputData = temp.Data
    with open(os.path.join(curr_dir,'Results',regFileName),'w') as g:
        # Initializing .reg file
        g.write('# Region file format: DS9 version 4.1\n')
        g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        g.write('image\n')
        for i in range(len(outputData)):
            Obj = outputData[i]
            Line = 'ellipse('+str(Obj['X_IMAGE'])+','+str(Obj['Y_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['A_IMAGE'])+','+str(Obj['KRON_RADIUS']*Obj['B_IMAGE'])+','+str(Obj['THETA_IMAGE'])+')\n'
            g.write(Line)
            
def createRegFile_outliers(outlierDict,outlierInds,regFileName):
    # Create a .reg file to locate outliers
    with open(regFileName,'w') as g:
        g.write('# Region file format: DS9 version 4.1\n')
        g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        g.write('image\n')
        for i in outlierInds:
            Line1 = 'circle('+str(outlierDict['x'][i])+','+str(outlierDict['y'][i])+','+str(outlierDict['aperture_radius'][i])+')\n'
            Line2 = 'circle('+str(outlierDict['x'][i])+','+str(outlierDict['y'][i])+','+str(10.0*outlierDict['aperture_radius'][i])+')\n'
            #Line3 = '# Flux Ratio = '+str(outlierDict['fluxRatio'][i]) 
            g.write(Line1)
            g.write(Line2)
            #g.write(Line3)

def chiSquareNormalized(Obs,Exp,numParams):
    ''' Compute the "reduced" or "normalized" chi-square value given a set of observed values, 
    a set of expected values, and the number of fitted parameters in the model of interest
    - This program assumes all observations are from the same population, and thus have the same
    variance
    - Assumes that entries correspond between the Obs and Exp sets
    '''
    
    if len(Obs) != len(Exp):
        print 'Error: sizes of Obs and Exp must match'
        
    df = float(len(Obs) - numParams) # Degrees of freedom = N - n
    var = np.var(Obs) # Sample variance of observations
    
    chiSquareNorm = np.sum(np.square(Obs - Exp))*(1/var)*(1/df) # Compute normalized chi-square
    return chiSquareNorm
        


def linReg3D(x1,x2,y):
    # Linear least-squares model fitting for 3D data
    A = np.column_stack((np.ones(np.shape(x1)),x1,x2))
    coeffs,residuals,rank,s = np.linalg.lstsq(A,y)
    return [coeffs,residuals,rank,s]
    

def rSquaredAdjusted(Obs,Exp,numParams):
    # Computes adjusted R^2 value from linear model residuals and data
    y,residuals = Obs, Obs - Exp
    N,p = y.size,numParams # Number of data points, number of 
    rSq = 1.0 - np.sum(np.square(residuals))/np.sum(np.square(y-y.mean()))
    rSqAdj = 1.0 - ((1.0 - rSq)*(float(N - 1.0)))/(float(((N - p) - 1.0)))
    return rSqAdj

#testImage = 'AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits'
#imgDict = binImage(getPixelValues(testImage),M=4,N=4)

#testCat = '/Users/annepstein/Work/Deconv/Results/G_S82_outerItr2_test7_Robust16_1.25FilteredUpdate_0.01PSFTresh_1.0Blur_noSR_150pxMask_wobgEst_blurPSF0.2_psfSize200_noSR_G_S82_outerItr2_test7_Robust16_1.25FilteredUpdate_0.01PSFTresh_1.0Blur_noSR_150pxMask_wobgEst_blurPSF0.2_psfSize200_noSR_compare.cat'
#testData = data(testCat)