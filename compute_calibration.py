# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:34:07 2017

@author: charlesgulian
"""

from astropy.io import fits

# Load all image filenames into dictionary
dataDict =  dict.fromkeys([line.strip() for line in open('/Users/annepstein/Work/Deconv/sample_image_filenames.txt','r').readlines()])
sanity=[]
calibrationData = '/Users/annepstein/Work/Deconv/CamCol2_calibration_data.csv'
with open(calibrationData) as g:
    _ = g.readline()
    for line in g.readlines():
        bits = line.split(',')
        curr_file = 'fpC-%06d-r2-%04d.fit' % (int(bits[4]), int(bits[7]))
        if curr_file in dataDict:
            dataDict[curr_file] = {'aa_r':float(bits[0]),
                                   'kk_r':float(bits[1]),
                                   'airmass_r':float(bits[2]),
                                   'sky_r':float(bits[3])}
            sanity.append(curr_file)
            
assert len(sanity) == len(dataDict)

# Write new header keywords to image headers
for image_file in dataDict.iterkeys():
    header = fits.getheader(image_file)
    for keyword in dataDict[image_file].iterkeys():
        header[keyword] = dataDict[image_file][keyword]
    


#'aa_r	kk_r	airmass_r sky_r	run	rerun	 camcol	field'

#Sample image name: fpC-006504-r2-0250.fit

