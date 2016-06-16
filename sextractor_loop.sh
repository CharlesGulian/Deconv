#!/bin/bash
# Iteratively runs sextractor on all .FITS files in current directory; creates .cat catalog file and .FITS check image of the same name as original .FITS file (check image name appended with _checkimg)

echo > img_name.txt
# Clear current contents of img_name
for f in AstroImages/*.fits
do
	FNAME=$(echo $f | sed 's/.fits//' | sed 's/AstroImages\///')
	# Copy name of .FITS file to variable FNAME
	
	echo $FNAME >> img_name.txt
	# Append filename to img_name.txt

	MG_ZRO=$(python get_magzero.py)
	# Get MAG_ZEROPOINT parameter from .fits file header (Python script)

	sex $f -c default.sex  -CATALOG_NAME Results/$FNAME.cat -MAG_ZEROPOINT $MG_ZRO -CHECKIMAGE_TYPE APERTURES -CHECKIMAGE_NAME Results/$FNAME_checkimg.fits
	# Run SExtractor on current image with preferences
done

./create_regfile.py
