#!/bin/bash
# Iteratively runs sextractor on all .FITS files in current directory; creates .cat catalog file and .FITS check image of the same name as original .FITS file (check image name appended with _checkimg)

echo > img_name.txt
# Clear current contents of img_name
for f in AstroImages/*.fits
do
	cp default.sex config.sex
	# Create a temporary configuration file for each image (this technique can be used to customize any part of config file)

	FNAME=$(echo $f | sed 's/.fits//' | sed 's/AstroImages\///')
	# Copy name of .FITS file to variable FNAME
	
	echo $FNAME >> img_name.txt
	# Append filename to img_name.txt

	sed -i '' "s/\[\[catalog\]\]/Results\/${FNAME}/" config.sex
	# Rewrite .cat filename to match original .FITS filename

	sed -i '' "s/\[\[check\]\]/Results\/${FNAME}_checkimg/" config.sex
	# Rewrite check image filename to match original .FITS filename with "_checkimg" appended

	MG_ZRO=$(python FITSheader.py)
	sed -i '' "s/\[\[MG_ZERO\]\]/${MG_ZRO}/" config.sex

	sex $f -c config.sex
	# Run SExtractor on current image
done

./create_regfile.py
