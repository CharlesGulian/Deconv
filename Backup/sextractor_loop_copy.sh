#!/bin/bash
# Iteratively runs sextractor on all .FITS files in current directory; creates .cat catalog file and .FITS check image of the same name as original .FITS file (check image name appended with _checkimg)


for f in *.fits
do
	cp default.sex config.sex
	# Create a temporary configuration file for each image (this technique can be used to customize any part of config file)
	FNAME=$(echo $f | sed 's/.fits//')
	# Copy name of .FITS file to variable FNAME
	sed -i '' "s/\[\[catalog\]\]/${FNAME}/" config.sex
	# Rewrite .cat filename to match original .FITS filename
	sed -i '' "s/\[\[check\]\]/${FNAME}_checkimg/" config.sex
	# Rewrite check image filename to match original .FITS filename with "_checkimg" appended
	sex $f -c config.sex
	# Run SExtractor on current image
done
