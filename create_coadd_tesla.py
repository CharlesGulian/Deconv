# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 16:29:44 2016

@author: charlesgulian
"""

# Co-adding program

# Load a list of images; create a co-added image whose pixels equal the median
# pixel of the corresponding pixels of each image in the list

import os
os.chdir('/home/cgulian2/Deconv')
curr_dir = os.getcwd()
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


# ===============================================================================
# Getting image file paths

# Image directory:
image_dir = '/home/DATA/STRIPE82_330-360_AlignCropped/test7'
image_files = glob.glob(os.path.join(image_dir,'*alignCropped.fits'))
image_files.remove('/home/DATA/STRIPE82_330-360_AlignCropped/test7/fpC-4927-x4127-y118_stitched_alignCropped.fits')
#print type(image_files)
#print image_files[0:5]

# ===============================================================================
# Generating co-add

# Do not load more than 0.5 GB of data at a time on Tesla machine
max_dir_size = 0.5 # GB
print 'Specified maximum allowed memory usage: {0} GB'.format(max_dir_size)
# Find size of directory
image_file_sizes = []
for i in range(len(image_files)):
    file_size = os.path.getsize(image_files[i]) # Bytes
    # Convert to gigabytes:
    temp = float(file_size)/float((2**30))
    image_file_sizes.append(temp)
directory_size = sum(image_file_sizes)
print '\nComputing co-added image from {0} frames'.format(len(image_files))
print 'Cumulative size of image files: {0} GB\n'.format(directory_size)

# ===============================================================================
# Computing optimal allocation of memory for computing co-add
image_dimensions = [1600,1600]
if directory_size > max_dir_size:
    print 'Cannot load all images simultaneously; must compute co-add iteratively'
    num_pix = image_dimensions[0]*image_dimensions[1]
    
    # Computing maximum NxN pixel sub-region of images that can be loaded simultaneously to compute co-add
    # * Note that this calculation assumes the file size of each image corresponds directly to # pixels in image
    
    proportion = max_dir_size/directory_size
    # At any time, we can load this proportion of the pixels in each image
    # It is more secure to load slightly less than this number of pixels at any time
    # To make computations easier, we will load the nearest square factor of 1600x1600 (num_pix)
    M = 1.
    while 1./(M**2) >= proportion:
        M += 1.
    N = np.ceil(float(image_dimensions[0])/M)
    num_pix_maxload = N*N
    if num_pix_maxload > proportion*num_pix:
        print 'Error: maximum memory usage exceeded'
    
else:
    M = 1.
    N = image_dimensions[0]

# ===============================================================================
# Creating co-add

# From here, write a loop to iteratively load and save an NxN portion of the image
# Take median of all pixels iteratively and save bins of co-add
# Then stitch together to create co-add

# Define empty array for new co-added image
coadd_image = np.zeros(image_dimensions)

# Define type of co-add (mean vs. median)
MEDIAN = True
MEAN = not(MEDIAN)       
if MEDIAN:
	op = np.median
       	ext = 'median'
elif MEAN:
	op = np.mean
	ext = 'mean'

# Iterate through each of MxM regions of image:
temp = int(M)
M = temp
for i in range(M):
    for j in range(M):
        
        # Create dictionaries for storing data
        index_dict = {}
        imageBin_dict = {}
        
        # Define indices by which to slice the (i,j)th bin of image
        xBinSize,yBinSize = N,N
        # Save indices
        indices = [int(np.floor(float(i)*xBinSize)),int(np.ceil(float(i+1)*xBinSize)),int(np.floor(float(j)*yBinSize)),int(np.ceil(float(j+1)*yBinSize))]
        index_dict[i,j] = indices
        
        for p in range(len(image_files)):
            # Select image
            image_file = image_files[p]
            # Get data
            image = fits.getdata(image_file)
            # Slice out the (i,j)th bin
            imageBin = image[indices[0]:indices[1],indices[2]:indices[3]]
            # Save the data from this bin
            imageBin_dict[p] = imageBin
            # Delete image data to save memory
            del image
        
        # Define empty array for (i,j)th bin
        coadd_bin = np.zeros([indices[1]-indices[0],indices[3]-indices[2]])
        # Create pixel-wise co-add
        for v in range(indices[1]-indices[0]):
            for w in range(indices[3]-indices[2]):
                # Create vector for (v,w)th pixel of each image
                pixel_vector = []
                for p in range(len(image_files)):
                    pixel = imageBin_dict[p][v,w]
                    pixel_vector.append(pixel)
		coadd_pixel = op(pixel_vector)
		coadd_bin[v,w] = coadd_pixel

        coadd_image[indices[0]:indices[1],indices[2]:indices[3]] = coadd_bin

print '\nCo-added image complete'
print 'Co-added image dimensions: ',np.shape(coadd_image)
coadd_image_file = os.path.join(curr_dir,'AstroImages','Coadd','custom_coadd_{0}.fits'.format(ext))
fits.writeto(coadd_image_file,coadd_image,clobber=True)
