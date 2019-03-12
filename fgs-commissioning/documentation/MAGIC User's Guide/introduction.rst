Introduction
============

The Multi-Application Guiding Interface for Commissioning (MAGIC) is a Python package that provides convenient access to a set of tools that will be used, as the name suggests, with the JWST FGS during OTE Commissioning. The package allows for user interaction with commissioning data and creates files that are needed for the operation of the flight software and the execution of visits.

These tools comprise of four main components that can be run individually
or together:

NIRCam to FGS image conversion (Image Converter)
------------------------------------------------
This tool can take in a simulated (or real) NIRCam image and will convert
it to an FGS (guider 1 or guider 2) image. In addition to rotating the image,
adjusting the pixel scale and image size, this tool corrects bad pixels and
normalizes the image to a specific magnitude of star.

Star Selection Tool (Star Selector)
-----------------------------------
This tool will take the FGS image either created with the first tool, or
an FGS image that it is passed by the user, and allow the user to choose
the guide star and reference stars using a GUI.

Flight Software File Writer (Flight Software (FSW) File Writer)
---------------------------------------------------------------
This module requires an FGS image and a file that includes a list of the
coordinates of the guide star and reference stars, along with their count
rates. This tool will create all files necessary to test this image different
flight software simulators (FGS DHAS, FGSES, etc.) These include all the
files necessary to run the ID, ACQ1, ACQ2, and TRK steps in these simulators.

Segment Guiding Tool (Segment Guiding)
--------------------------------------
Used to facilitate guiding on unstacked segment arrays during commissioning. When
provided 1) the commanded RA and Dec of a guide star and 2) the V2/V3 (or x/y)
positions of all segments in an array, the segment guiding tool calculates the
effective RA and Dec of all segments on the sky.

----------------------------------------------------

The following documents will walk you through how to use each of the parts of MAGIC and verify those results. 
Note: A powerful feature of MAGIC is that is retains information about the image you are working on as long as the window stays open and the file is loaded, so we recommend that while working with MAGIC, you keep the GUI open. If it crashes (which can happen) don’t worry, just following section II again and use the same input values and MAGIC will find the dataset that you are using. 

----------------------------------------------------

“We’re not *really* going to use magic?” – Ronald Weasley

