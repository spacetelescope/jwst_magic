FGS Commissioning Tools
=======================

![GUI](notebooks/FGSCommTools_GUI.png)

The FGS Commissioning Tools package provides convenient access to  numerous ancillary tools that will be used, as the name suggests, with the JWST FGS during OTE Commissioning. The package allows for user interaction with commissioning data and creates files that are needed for the operation of the flight software and the execution of visits.

These tools comprise of four main components that can be run individually
or together:

### 1. NIRCam to FGS image conversion (``convert_image``)
This tool can take in a simulated (or real) NIRCam image and will convert
it to an FGS (guider 1 or guider 2) image. In addition to rotating the image,
adjusting the pixel scale and image size, this tool corrects bad pixels and
normalizes the image to a specific magnitude of star.


### 2. Star Selection Tool (``star_selector``)
This tool will take the FGS image either created with the first tool, or
an FGS image that it is passed by the user, and allow the user to choose
the guide star and reference stars using a GUI.


### 3. Flight Software File Writer (``fsw_file_writer``)
This module requires an FGS image and a file that includes a list of the
coordinates of the guide star and reference stars, along with their count
rates. This tool will create all files necessary to test this image different
flight software simulators (FGS DHAS, FGSES, etc.) These include all the
files necessary to run the ID, ACQ1, ACQ2, and TRK steps in these simulators.


### 4. Segment Guiding Tool (``segment_guiding``)
Used to facilitate guiding on unstacked segment arrays during commissioning. When
provided 1) the commanded RA and Dec of a guide star and 2) the V2/V3 (or x/y)
positions of all segments in an array, the segment guiding tool calculates the
effective RA and Dec of all segments on the sky.


Installation notes
------------------
This package was developed in a python 3.5 environment. Python 2.7 is not yet supported.

The following supplemental packages are required, and are available for easily installation with `pip`:
* `pysiaf`
* `photutils`

##### To install:

1) Clone the gitlab repository to your local machine

	git clone git@grit.stsci.edu:wfsc/tools.git

2) Install the `jwst_fgs_commissioning_tools` package:

	cd tools/fgs-commissioning

	pip install -e .a



Running the Tools
-----------------
These tools are best run in the IPython terminal. Simply launch the GUI with the following steps:

    In[1]: import jwst_fgs_commissioning_tools

    In[2]: jwst_fgs_commissioning_tools.run_tool_GUI()


Tutorial
-----------------
See the [`FGS Commissioning Tool Tutorial.ipynb`](notebooks/FGS Commissioning Tool Tutorial.ipynb) notebook for examples on how to run the tools on a modular level from within a Jupyter notebook.


###### Developed by Keira Brooks, Lauren Chambers, Kathryn St. Laurent and collaborators, 2016-2018.
