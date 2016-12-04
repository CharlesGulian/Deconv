# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 19:03:41 2016

@author: charlesgulian
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress as linreg

   
def analyze(x,y,output_file):
    
    # =======================================================================
    # Splitting up data into N different brightness regimes
            
    N = 5 # Number of brightness regimes (each corresponding to adjacent partitions of width b*std(flux1) of the flux1 distribution; final partition is all points greater than (N-1)*std(flux1))
    
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
    
    inds_unsorted = np.arange(len(x))
    x_sorted,inds_sorted = zip(*sorted(zip(x,inds_unsorted)))
    L = len(x)
    bin_size = L/N

    regimes = {}
    regimes[0] = []
    inds = np.arange(len(x))
    regimes[0].extend(inds) # Oth regime corresponds to original (full) dataset
    for i in range(N-1):
        regimes[i+1] = []
        inds = np.array(inds_sorted[i*bin_size:(i+1)*bin_size-1])
        regimes[i+1].extend(inds)
    regimes[N] = []
    inds = np.array(inds_sorted[(N-1)*bin_size::])
    regimes[N].extend(inds) 

    # =======================================================================
    # Testing null hypothesis that data has slope m = 1.0
    residual = y - x
    y = residual
    
    # =======================================================================
    # Writing output file
    with open(output_file,'w') as g:
        g.write('# Linear regression analysis results\n')
        g.write('# Columns: Slope | Intercept | r-value | p-value \n')
        g.write('\n')
        for i in range(N+1):
            
            # =======================================================================
            # Slice ith brightness regime of data
            xr,yr = x[regimes[i]],y[regimes[i]]

            # =======================================================================
            # Do linear regression analysis with data in x,y
            slope, intercept, r_value, p_value, std_err = linreg(xr,yr)
            
            # =======================================================================
            # Write data to output file
            xl,xh = xr[0],xr[len(xr)-1] # Bounds of brightness regime
            g.write('# ({0}) Results for data in x-axis interval {1}:{2} (n={3})\n'.format(i+1,xl,xh,len(regimes[i])))
            g.write('{0} {1} {2} {3}\n'.format(slope,intercept,r_value,p_value))
            if p_value < 0.05:
                g.write('# Statistically significant deviation from slope of 1.0\n')
            g.write('\n')
    

'''
noise1 = (0.1*(np.random.rand(400)-0.5))
noise2 = 0.2*np.random.randn(400)
artifact1 = np.zeros(400)
artifact1[0:56] = -0.6*np.linspace(0.0,1.0,56)+0.9

for i in range(5):
    
    curr_dir = os.getcwd()
    output_file = os.path.join(curr_dir,'Results','LinRegAnalysis','test{0}.txt'.format(i+1))
    
    # Artificial data
    m = 1.0 + 0.025*i # Slope of data
    
    x = np.linspace(4.0,16.0,400) + noise1
    y = m*x + noise2 + 0.5
    
    xsort,ysort = zip(*sorted(zip(x,y)))
    
    x,y = np.array(xsort),np.array(ysort)
    y += artifact1
    
    SHOW = False
    if SHOW:
        plt.plot(np.linspace(0.0,18,400),np.linspace(0.0,18,400),'b--')
        plt.scatter(x,y,c='m',linewidth=0.1)
        plt.axis([0.0,19.0,0.0,19.0])
        plt.show()

    analyze(x,y,output_file)
#'''