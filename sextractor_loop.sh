#!/bin/bash
# Iteratively runs sextractor on all .FITS files in current directory; creates .cat catalog file and .FITS check image of the same name as original .FITS file (check image name appended with _checkimg)

echo > img_name.txt
# Clear current contents of img_name
for f in AstroImages/*.fits
do
	cp default.sex config.sex

	FNAME=$(echo $f | sed 's/.fits//' | sed 's/AstroImages\///')
	# Copy name of .FITS file to variable FNAME
	
	echo $FNAME >> img_name.txt
	# Append filename to img_name.txt

	sed -i '' "s/\[\[check\]\]/Results\/${FNAME}_checkimg/" config.sex
	# For some reason doing this in the call to sex $f below doesn't work

	MG_ZRO=$(python get_magzero.py)
	# Get MAG_ZEROPOINT parameter from .fits file header (Python script)

	sex $f -c config.sex  -CATALOG_NAME Results/$FNAME.cat -MAG_ZEROPOINT $MG_ZRO
	# Run SExtractor on current image
done

./create_regfile.py
