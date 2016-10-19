# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 23:59:22 2016

@author: charlesgulian
"""

# Kernel Density Estimation (KDE) with SciPy
# (This code is borrowed from https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/)

from scipy.stats import gaussian_kde

def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)
    
    
# Note that "x_grid" is actually a 1-D vector along which the curve is plotted;
# x is the data whose "density" we wish to estimate