#!/bin/bash
# Iteratively runs sextractor on all .FITS files in current directory; creates .cat catalog file and .FITS check image of the same name as original .FITS file (check image name appended with _checkimg)

#echo > img_name.txt
# Clear current contents of img_name
for f in AstroImages/*.fits
do
	# Copy default versions of configuration file and output catalog parameter file
	cp default.sex copy.sex
	cp default.param copy.param

	# Rewrite configuration and output catalog parameters in Python script do_config.py
	./do_config.py $f copy.sex copy.param
	
	######
	#echo $FNAME >> img_name.txt
	# Append filename to img_name.txt	
	######

	sex $f -c copy.sex
	# Run SExtractor
	#exit
done

#./create_regfile.py
