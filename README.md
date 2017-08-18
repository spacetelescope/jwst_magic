# WFSC Tools for use during Ground Testing

## 1. FIX_ISIM
The fix_isim tool is used on ground calibrated data (e.g. ISIM DHAS) to produce WAS-compatible header information as well as 
DMS-equivalent image orientation. 

The script uses command-line arguments and expect at least one argument: the path to the directory containing the FITS images to 
fix. If no other arguments are provided, the tool will process all the files in that folders, and output duplicated updated 
files with suffix _FIXED.  Here are the arguments that can be provided to the tool:

* dir: full/relative path to the input images to be fixed.
* -s : [OPTIONAL] string to match when processing the files (e.g. "-s NRCA" will process all the files containing NRCA). Default 
behavior is to process all the files in the directory.
* -o : [OPTIONAL] full/relative path to save the output images. Default value is the input directory.

Example usage: python fix_isim.py /Users/Lajoie/Desktop/OTIS_Images/ -s NRC310 -o ./Fixed_NRC310

For more information: python fix_isim.py --help


## 2. LVDT TRACKING TOOL
The LVDT Tracking tool is used to monitor the difference between the expected MCS LVDT value and the returned ADU LVDT value. THe 
input to the tool is the main directory where all your unzipped archive packages live. Before running the notebook, you should 
also verify that the MSDB and WFC Record file names match what's in your archive packages.

Here is a screenshot of a typical output for the tool:

![alt text](https://grit.stsci.edu/wfsc/tools/blob/jsc-testing/LVDT_Tool_preview.png "LVDT Tool Preview")
