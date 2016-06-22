#!/bin/bash

NPARAM=2
if [ $# -eq "$NPARAM" ]
then	
	# ~~~First comparison: first image with first image~~~

	# Copy and modify configuration file, get .fits image tags
	cp default.sex config.sex
	FNAME1=$(echo $1 | sed 's/.fits//' | sed 's/AstroImages\///')
	FNAME2=$(echo $2 | sed 's/.fits//' | sed 's/AstroImages\///')	

	# Copy .fits image tags to .txt file for analysis in Python
	echo $FNAME1 > img_name_compare.txt
	echo $FNAME2 >> img_name_compare.txt
	
	sed -i '' "s/\[\[check\]\]/Results\/${FNAME1}_compare_checkimg/" config.sex
	sed -i '' "s/\[\[catalog\]\]/Results\/${FNAME1}_${FNAME1}_comparison/" config.sex	

	MG_ZRO=$(python get_magzero.py)
	# Get MAG_ZEROPOINT parameter from .fits file header (Python script)	

	sex -c config.sex $1,$1 -MAG_ZEROPOINT $MG_ZERO
	# SExtract first image, compare with first image

	# ~~~Second comparison: first image with second image~~~

	# Copy and modify configuration file  
        cp default.sex config.sex
        sed -i '' "s/\[\[check\]\]/Results\/${FNAME2}_compare_checkimg/" config.sex
	sed -i '' "s/\[\[catalog\]\]/Results\/${FNAME1}_${FNAME2}_comparison/" config.sex	

        MG_ZRO=$(python get_magzero.py)
        # Get MAG_ZEROPOINT parameter from .fits file header (Python script)

	sex -c config.sex $1,$2 -MAG_ZEROPOINT $MG_ZERO
        # SExtract first image, compare with second image

else
	echo "Input two .FITS files to be compared"
	exit
fi
