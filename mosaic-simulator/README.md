# WFSC Ancillary Tools for Commissioning

## 1. Mosaic Simulator
The Mosaic Simulator can be used to generate images with WCS header information consistent with the different observations of 
OTE-01 Initial Image Mosaic (APT program 1134):

![alt text](https://grit.stsci.edu/wfsc/tools/blob/master/mosaic-simulator/OTE-01-Observations.png "Initial Image Mosaic Observations")

The necessary datafiles can be downloaded from central store: ***REMOVED***/lajoie/InitialMosaicData/


The sequence of events as the tool runs is as follows:
1. Randomly determine a star position and, given the uncertainty on the boresight and the segment deployment, draw 18 random 
positions for the mirror segments around that target position;
2. Using the pointings information from APT Program 1134 (see attached txt file), read the dither/offset of each visit and loop 
and determine if a segment falls on any of the NIRCam detectors at that particular position;
3. If a segment is "captured" by NIRCam, stich in PSF image, else stich in a dark (see attached FITS 
files). 
4. WCS information is also incorporated to the image header information at that point, and the image is saved in ObsX folder as 
testXXXX.fits.

The tool has the logic to determine is overlapping tiles capture the same segment(s), and if so, to ensure it appears on all 
overlapping images. This is important because of the way QUIP displays overlapping images.

The user should specify in the script which observation number is desired. Values range from 1 to 7 inclusively (see image 
above).



## 2. Mosaic Congrid
This tool is a function taken out of the Mosaic Simulator and downsamples the input images into 100x100 images, which can then be 
fed to QUIP. The WCS header information needs to be modified in order for the scaling to be correct; in particular the CDELTX and 
CDX_Y keywords.
