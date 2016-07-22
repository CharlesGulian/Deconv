# Deconvolution Project Script Tutorial

This tutorial demonstrates the functions of the following scripts in this repository:

- MainCompare.py
- sex_config.py
- sex_stats.py
- pysex.py
- fits_tools.py


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

* sex_config.py solely contains one class: ```configure```:
``` python
class configure:
    
    default_params = []
    def __init__(self,image_file,config_file,param_file,dual=False,default_params=default_params):
            
            self.config_file,self.param_file = config_file,param_file
            self.dual = dual
            
            if self.dual:
                print 'Configuring SExtractor for dual image mode'
                self.image_file1,self.image_file2 = (image_file).split(',')
            else:
                self.image_file = image_file
                    
            self.config_dict = {}
            self.param_dict = {}
            ...
            ...
```
* The class is initiated as follows
``` python
fig = sex_config.configure(image_file.fits,config_file.sex,param_file.param)
```
or, in dual image mode, we write the image file as two image files separated by a comma
``` python
fig = sex_config.configure(image_file1.fits+','+image_file2.fits,config_file.sex,param_file.param,dual=True)
```
* In MainCompare.py, we create two instances of ```configure``` (```fig1``` and ```fig2```) in dual image mode, such that we compare image_file1 with itself, and then with image_file2
* The ```configure``` class has a method called ```reconfigure```, which takes the name of a configuration parameter and its value as input (e.g. ```fig.reconfigure('BACK_TYPE','MANUAL')```
    * Internally, the ```configure``` class stores and edits the value/state of every configuration parameter and output catalog parameter in two separate dictionaries
* To save space, most of the basic configuring is handled in the ```configure.default()``` method, but parameter values can still be changed from MainCompare.py (perhaps in an iterative fashion) 

* After executing ```fig.default()``` in MainCompare.py, we use the the ```configure.write_configfile``` method to write new configuration and output parameter files

## pysex.py

* pysex.py is a simple SExtractor wrapper that uses Python's ```subprocess``` module to call SExtractor via the command line
* The ```compare()``` function compares two images using their respective configuration files, as in MainCompare.py:
``` python
 pysex.compare(image1,image2,'copy_compare1.sex','copy_compare2.sex')
```
 

