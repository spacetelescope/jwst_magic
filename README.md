FGS Commissioning Tools
=======================

*Python 3*

These tools comprise of four main components that can be run individually
or together.

1. NIRCam to FGS image conversion -
This tool can take in a simulated (or real) NIRCam image and will convert
it to an FGS, guider 1 or guider 2 image. In addition to rotating the image, adjusting the pixel scale and
image size, this tool corrects bad pixels and normalizes the image to
a specific magnitude of star.


2. Star Selection Tool -
This tool will take the FGS image either created with the first tool, or
an FGS image that it is passed by the user, and allow the user to choose
the guide star and reference stars using a simple GUI.


3. File creation -
This part requires an FGS image and a file that includes a list of the
coordinates of the guid star and reference stars along with their count
rates. This tool will create all files necessary to test this image different
flight software simulators (FGS DHAS, FGSES, etc.) These include all the
files necessary to run the ID, ACQ1, ACQ2, and TRK steps in these simulators.


4. Segment Guiding Tool -
*In progress*

Installation notes
------------------

This package was developed in a python 3.5 environment. Python 2.7 is not yet supported.

Clone the repository

``git clone git@grit.stsci.edu:JWST-FGS/Commissioning-tools.git``

Install the JWST FGS Commissioning Tools:

``cd Commissioning_tools/``

``python setup.py install --user``



Running the Tools
-----------------
These tools are best run in IPython or Jupyter Notebook.

* Jupyter Notebook
See the `FGS Commisioning Tool Tutorial.ipynb` notebook for examples on how
to run the first three tools.

* IPython
In IPython and execute the following steps:

```
In[1]: import run_fgs_commissioning_tool

In[2]: input_image = ‘full/path/to/image.fits’

In[3]: guider = 1

In[4]: run_fgs_commissioning_tool.run_all(input_image, guider)
```
