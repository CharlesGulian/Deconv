#!/bin/bash
# Compares two .fits images using SExtractor in dual-image mode; first run detects from first image, measures from first image;
# second run detects from first image, measures from second image; output catalogs of both runs are saved


NPARAM=2
if [ $# -eq "$NPARAM" ]
then	

	# Make copies of default configuration file and output catalog parameter file, reconfigure in Python script do_config.py
	cp default.sex copy_compare.sex
	cp default.param copy_compare.param
	./do_config.py $1,$1 copy_compare.sex copy_compare.param -d	
	
	# Get .fits image tags
	FNAME1=$(echo $1 | sed 's/.fits//' | sed 's/AstroImages\///')
	FNAME2=$(echo $2 | sed 's/.fits//' | sed 's/AstroImages\///')	
	
	# Copy .fits image tags to .txt file for analysis in Python
	echo $FNAME1 > img_name_compare.txt
	echo $FNAME2 >> img_name_compare.txt
	
	sex -c copy_compare.sex $1,$1
	# Detect from first image, measure from first image
	
	# Make copies of default configuration file and output catalog parameter file, reconfigure in Python script do_config.py
        cp default.sex copy_compare.sex
        cp default.param copy_compare.param
	./do_config.py $1,$2 copy_compare.sex copy_compare.param -d

	sex -c copy_compare.sex $1,$2
	# Detect from first image, measure from second image
else
	echo "Input two .FITS files to be compared"
	exit
fi
