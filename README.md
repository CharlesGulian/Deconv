# Deconvolution Project: Python Script Tutorial

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

### sex_config.py

* sex_config.py contains a single class, ```configure```:
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

### pysex.py

* pysex.py is a simple SExtractor wrapper that uses Python's ```subprocess``` module to call SExtractor via the command line
``` python
import os
import subprocess

sex_path = '/usr/local/bin/sex'
curr_dir = os.getcwd()

# All image filenames must be absolute paths or relative to current working directory

def _call(args):
    try:
        subprocess.call(args,shell=False)
    except subprocess.CalledProcessError:
        print 'Error: sex_call CalledProcessError'
        pass # Handle errors in the called executable
    except OSError:
        print 'Error: OSError'
        pass # Executable not found

def call_sex(imgName,config_file=None,args_ext=[]):
    # Call SExtractor from Python script
    
    args = [sex_path,imgName]
    
    if not os.path.exists(imgName):
        if os.path.exists(os.path.abspath(imgName)):
            imgName = os.path.abspath(imgName)
        else:
            imgName = os.path.join(curr_dir,imgName)
        args[1] = imgName
        
    if config_file != None:
        if not os.path.exists(config_file):
            config_file = os.path.join(curr_dir,config_file)
        args.extend(['-c',config_file])
        
    args.extend(args_ext)
    _call(args)
...
...
```
*  ```compare()``` runs SExtractor in dual image mode for two images using their respective configuration files, as in MainCompare.py:
``` python
 pysex.compare(image1,image2,'copy_compare1.sex','copy_compare2.sex')
```
* This generates two SExtractor output catalogs containing information on the detected sources' positions, shapes, and brightnesses

### sex_stats.py

* sex_stats.py contains a class called ```data```, which is used to extract and organize data from SExtractor output catalogs
* It is initialized as follows
```
imageData = sex_stats.data(output_catalog.cat)
```
* The ```data``` class has a method called ```get_data()```, which takes a SExtractor output parameter name as input and returns the corresponding column of the output catalog as a NumPy array; e.g.
``` python
flux_radius_array = imageData.get_data('FLUX_RADIUS')
```
* sex_stats.py also contains functions to bin image/array data (```binImage()```) and scattered 3-D data (```binData()```), as well as a graveyard for the statistics functions I wrote

* In MainCompare.py, we use the ```sex_stats.data``` class to retrieve and analyze information about sources detected via SExtractor
``` python
 # Get image "tags" (filenames without filetype extensions)
img_tag1 = os.path.split(image1)[1]
img_tag1 = img_tag1[0:len(img_tag1)-len('.fits')]
img_tag2 = os.path.split(image2)[1]
img_tag2 = img_tag2[0:len(img_tag2)-len('.fits')]

# Get first output catalog filename
outputCat1 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag1+'_compare.cat')
if not os.path.exists(outputCat1):
    print 'Error: first output catalog path does not exist'
# Get second output catalog filename
outputCat2 = os.path.join(os.getcwd(),'Results',img_tag1+'_'+img_tag2+'_compare.cat')
if not os.path.exists(outputCat2):
    print 'Error: second output catalog path does not exist'

# Create sex_stats.data objects:
img1data = sex_stats.data(outputCat1)
img2data = sex_stats.data(outputCat2)
```

### fits_tools.py

* fits_tools.py is a collection of functions that I've written to manipulate .fits images
* fits_tools.py comprises the following functions
    * Binning images: ```binImage()```
    * Retrieving .fits image data: ```getPixels()```
        * Note that ```astropy.io.fits.get_data()``` does the same thing as ```getPixels()```
    * Masking a .fits image with a .fits binary mask: ```maskImage()```
    * Bias and median subtraction: ```subtractBias()``` and ```subtractMedian()```
    * Computing the circular-aperture flux of an object: ```computeObjectFlux()```
* E.g.
``` python
# Assign .fits image path to image_file variable   
image_file = 'image_file.fits'
# Subtract bias of 1000.0 from image_file data, write to new image file
subtractBias('image_file.fits',new_image_file='image_file_biasSubtracted.fits',bias=1000.0)
# Reassign new (bias-subtracted) image to image_file variable
image_file = 'image_file_biasSubtracted.fits'

# Assign .fits binary mask path to mask_file variable
mask_file = 'mask.fits'
# Apply .fits binary mask to image_file data, write to new .fits file
maskImage(image_file,mask_file,new_image_file='image_file_masked.fits')
# Reassign new (masked) image to image_file variable
image_file = 'image_file_masked.fits'

imageData = getPixels(image_file.fits)      # Retrieve image data array
imageBinDict = binImage(imageData,M=4,N=4)  # Bin image into 4x4 bins
lowerRight_data = imageBinDict[3,3]                 # Get data from bin in lower right corner of image

# Assuming a sex_stats.data instance exists for image_file's SExtractor output catalog (outputData):
x,y = outputData.get_data('X_IMAGE'),outputData.get_data('Y_IMAGE')
radii = outputData.get_data('FLUX_RADIUS')
# Compute the circular aperture flux for the 26th object in the output catalog
i = 25
flux = computeObjectFlux(x[i],y[i],radii[i],imageData)
```

### MainCompare.py (continued)

* The rest of MainCompare.py focuses on analyzing the flux ratios of detected sources
    * Note that the order/number of detected objects in the output catalogs will be the same for both comparison images, as the first image is used as the "detection image" for both calls to ```pysex.compare()```
* We retrieve data on the (x,y) coordinates and half-light radii of all detected sources and therefrom compute the sources' circular-aperture fluxes with ```fits_tools.computeObjectFlux()```
``` python
# Computing circular-aperture fluxes
flux1,flux2 = np.zeros(np.shape(x)),np.zeros(np.shape(x))
for i in range(len(x)):
    flux1[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData1)
    flux2[i] = fits_tools.computeObjectFlux(x[i],y[i],rads[i],imageData2)

# Compute object-wise flux ratio
fluxRatio = np.divide(flux1,flux2)
```
* We compute and print several statistics from the image data and output catalog data, and then create/save a number of figures and histogram to analyze the flux ratio data
    * These include a histogram of object fluxes from each image, a histogram of flux ratios, a figure containing MxN histograms of flux ratios (one corresponding to each bin), a color-grid representing the mean flux ratio of each bin, a colored scatter plot of flux ratios (one for each object), and a histogram of half-light radii
