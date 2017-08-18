# WFSC Tools for use during Ground Testing

## FIX_ISIM
The fix_isim tool is used on ground calibrated data (e.g. ISIM DHAS) to produce WAS-compatible header information as well as 
DMS-equivalent image orientation. 

The script uses command-line arguments and expect at least one argument: the path to the directory containing the FITS images to 
fix. If no other arguments are provided, the tool will process all the files in that folders, and output duplicated updated 
files with suffix _FIXED.  Here are the arguments that can be provided to the 
tool:

* dir: full/relative path to the input images to be fixed
* -s : string to match when processing the files (e.g. "-s NRCA" will process all the files starting with NRCA)
* -o : full/relative path to save the output images. Default is the input directory.

Example usage: python fix_isim.py /Users/Lajoie/Desktop/OTIS_Images/ -s NRC310 -o ./Fixed_NRC310

For more information: python fix_isim.py --help

## LVDT TRACKING TOOL
