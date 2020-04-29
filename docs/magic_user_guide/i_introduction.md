I. Introduction
===============

<p align="center">
  “We’re not <b>really</b> going to use magic?” – Ronald Weasley
</p>


The Multi-Application Guiding Interface for Commissioning (MAGIC) is a Python package that provides convenient access to a set of tools that will be used, as the name suggests, with the JWST FGS during OTE Commissioning. The package allows for user interaction with commissioning data and creates files that are needed for the operation of the flight software and the execution of visits. The documents here will walk you through how to use each of the parts of MAGIC and verify those results.

A powerful feature of MAGIC is that is retains information about the image you are working on as long as the window stays open and the file is loaded, so we recommend that while working with MAGIC, you keep the GUI open. If it crashes (which can happen) don’t worry, just go to [Section II](ii_setting_up.md), input the same values that you have been using for that data, and MAGIC will find the dataset that you are working with.

These tools comprise of four main components that can be run individually
or together:

NIRCam or FGS science to FGS raw image conversion (Image Converter)
------------------------------------------------
This module can take in a simulated or real NIRCam (all shortwave and long wave detectors) or FGS (all detectors) image and will convert it to a pseudo-FGS raw (guider 1 or guider 2) image. This pseudo-FGS raw image is not a simulated image, but an expectation of the raw CDS image that the flight software will see. This module works with the [`jwst-fgs-countrate`](https://github.com/spacetelescope/jwst-fgs-countrate) module to use the count rates in the input image to renormalize the image to the expected count rates for a different object. This module appropriately rotates the image from the science frame, adjusts the pixel scale and image size if converting from a NIRCam image, corrects bad pixels, and normalizes the image to the magnitude and count rates of a specified guide star.

[Section III: Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md)

Star Selection Tool (Star Selector)
-----------------------------------
This module takes in a raw FGS image and allows the user to choose the guide and reference star PSFs using a GUI.

[Section IV: Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

Flight Software File Writer (Flight Software (FSW) File Writer)
---------------------------------------------------------------
This module requires an FGS image and a file that includes a list of the image coordinates and count rates of the selected guide and reference star PSFs. The output is all files necessary to test this image with different flight software simulators (FGS DHAS, FGSES, etc.) includes all the files necessary to run the ID, ACQ1, ACQ2, and TRK steps in these simulators.

[Section IV: Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

Segment Guiding Tool (Segment Guiding)
--------------------------------------
Used to facilitate guiding on unstacked segment PSF arrays (Segment Override), and stacked, but un-phased segment PSFs (Photometry Override) during commissioning.

**Segment Override** - When provided with:
1. the commanded RA and Dec of a guide star, and
2. the V2/V3 (or x/y) positions of all segments in an array

the segment guiding tool calculates the effective RA and Dec of all segments on the sky. This information can then be used to override the guide star catalog with new RA, Dec, and count rate values for the requested guide and reference star PSFs

[Section VII: Writing the Segment Override File (SOF)](vii_write_sof.md)

**Photometry Override** - When provided with a count rate factor and a threshold factor, the tool will provide the necessary information to override the guide star catalog with a new count rate for a specific guide star

[Section VIII: Writing the Photometry Override File (POF)](viii_write_pof.md)

----------------------------------------------------

#### Next

###### Section II: [Setting Up MAGIC](ii_setting_up.md)

###### Section III: [Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md)

###### Section IV: [Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

###### Section V: [Testing Selections in DHAS](v_testing_in_dhas.md)

###### Section VI: [Contingency: Re-selecting Stars and Re-running DHAS](vi_contingency_reselect_stars.md)

###### Section VII: [Writing the Segment Override File (SOF)](vii_write_sof.md)

###### Section VIII: [Writing the Photometry Override File (POF)](viii_write_pof.md)

###### Appendix A: [Installing the JWST MAGIC Package](appendix_a_installing_magic.md)

###### Appendix B: [Setting Up DHAS](appendix_b_opening_dhas.md)

###### Appendix C: [Using APT to Get Guide Star RA & Dec](appendix_c_apt.md)

###### Appendix D: [Mirror State Procedures](appendix_d_mirror_states.md)
