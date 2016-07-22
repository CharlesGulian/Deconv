# Deconvolution Project Script Tutorial

This tutorial demonstrates the functions of the following scripts in this repository:

- MainCompare.py
- sex_stats.py
- sex_config.py
- fits_tools.py
- pysex.py

### MainCompare.py

* This is the "main comparison" script, i.e. where the actual image comparisons get done, and where the resulting data is collected and analyzed
* It should be run via the command line using ```python MainCompare.py```
* MainCompare.py directly or indirectly makes use of all of the other scripts listed in this tutorial; as such, most demonstrations of other scripts will pertain to MainCompare.py
* The structure of MainCompare.py is as folows:

* First, we define the ```Image``` class, which is used to save the filenames of different categories of .fits images (e.g. "Good," "Bad", "Co-added," and "Deconvolved" images)
``` python
class Image:
    
    def __init__(self,filename,category,ID):
        self.filename = filename
        self.category = category
        self.ID = ID
```
* Here is an example of an instance of this class
``` python
goodImage1 = Image('AstroImages/Good/fpC-6484-x4078-y134_stitched_alignCropped.fits','Good','1')
coaddedImage = Image('AstroImages/Coadd/fpC-206-x4684-y126_stitched_alignCropped-COADD.fits','Coadded','')
```
* The ```Image``` class instances are primarily used to save and categorize data analysis results/figures

* Next, we define a list of ```Image``` class instances and iterate through it in a nested loop such that every image is paired for comparison with every other image exactly once
* We perform trivial data cleaning operations; bias subtraction, median subtraction
    * The ```Image.category``` method can be used to selectively perform operations on certain categories of images; for example, we only subtract a bias from "Good" images
* In the next part of the comparison loop, we use sex_config.py to write custom SExtractor configuration files for each image

## sex_config.py

* Well,
