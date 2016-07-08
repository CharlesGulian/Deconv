# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 20:31:36 2016

@author: charlesgulian
"""
# Library of statistics functions for photometric statistics and data analysis

import os
os.chdir('/Users/annepstein/Work/Deconv')
import numpy as np
import pandas as pd

class data:
    
    def __init__(self,outputCatName):
        f = open(outputCatName)
        lines = f.readlines()
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


    
    