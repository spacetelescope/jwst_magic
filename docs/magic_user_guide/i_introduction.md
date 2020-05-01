I. Introduction
===============

Guiding During the JWST Optical Telescope Element commissioning
---------------------------------------------------------------
### Getting a raw FGS image that doesn't yet exist
The main usage for MAGIC is to override the guide star catalog for the cases of un-stacked and un-phased PSFs during the commissioning of the JWST Optical Telescope Element (OTE). This means giving the FGS flight software new count rates and/or locations of the guide and reference stars when we have 18 PSFs per guide star, or one PSF in which the segment PSFs are not phased (and therefore the count rate will be reduced).

To do this, we use the last NIRCam image taken (with hopefully the same mirror configuration), and convert it to what we expect from the FGS detector for guiding modes (ID, ACQ, and TRK) with the guide star. We call this our pseudo-FGS image. We put this image in the raw frame in order to simulate what the FGS flight software will be seeing when it tries to ID the guide star.

The target in the NIRCam detector and the FGS detector will rarely, if ever, be the same, so it is not only necessary to update the pixelscale and correct for any detector-specific differences, but we also have to account for a different count rate in each PSF. To do this, we can "renormalize" the input image to the count rate of the guide star (in the FGS detector), rather than the count rate of the target star (in the NIRCam detector), since we know the count rate of star that is now spread over multiple PSFs. We assume that the impact due to background stars will be insignificant due of the isolation criteria placed on the guide star candidates.  

#### A Note on Filters, jitter, and Exposure Time
The full frame FGS images (called calibration images) and NIRCam full frame images have similar exposure times of ~10s, however the NIRCam image is taken with a filter, and FGS images are not since there are no filters on the FGS. This difference manifests itself in what appears to be a blurring or smearing of each PSF. An important note is that when guiding, the FGS calibration image is not, in fact, used, instead strips (36 strips of 2048 X 64 pixels, with some overlap between strips) are used for the the Identification mode and these strips are taken with a much shorter exposure time (~1/3s). Therefore, we make the rough assumption that additional jitter due to longer exposure times in the full frame NIRCam image, is roughly equivalent to smearing due to additional wavelengths observed in the FGS ID, ACQ, and TRK images.


### Determining what to guide on
After creating our pseudo-FGS raw image, we can then hand-select the PSFs that we feel will be most successful as guide and reference stars. We will particularly want to do this when we have un-stacked segments since we will have up to 18 PSFs per star in the detector. Since the isolation criteria will hopefully mean that the guide star will be the only "bright" object in the detector, the reference stars will not be catalog stars, but other segment PSFs, and the "guide star" will be the "best" segment PSF candidate for guiding. When the segments are stacked, but not phased, the catalog reference stars can be used since the catalog RA and Dec are valid, but the count rate will not be what is in the guide star catalog.

### Testing
Once the guide and reference stars have been selected (and with that selection process we calculate the 3x3 count rate around the peak pixel in the selected PSF), simulated images of the ID, ACQ, and TRK modes can easily be created along with any files to help run these images through flight software simulator. This step will help verify that we can expect (or not expect) to successfully guide given these count rates and/or locations of stars in the detector.

### Overriding the guide star catalog
The final step is to use the values calculated for the count rate and/or locations of the guide and reference stars, to create a guide star catalog override command.


The Tools
----------
The Multi-Application Guiding Interface for Commissioning (MAGIC) is a Python package that provides convenient access to the set of tools that accomplishes the goals laid out above and will be used, as the name suggests, with the JWST FGS during OTE Commissioning. The package allows for user interaction with commissioning data and creates the necessary files needed for the operation of the flight software and the execution of visits. The documents here will walk you through how to use each of the parts of MAGIC and verify those results.

A powerful feature of MAGIC is that is retains information about the image you are working on as long as the window stays open and the file is loaded, so we recommend that while working with MAGIC, you keep the GUI open. If it crashes (which can happen), don’t worry, just go to [Section II](ii_setting_up.md), input the same values that you have been using for that data, and MAGIC will find the dataset that you are working with.

The components discussed above directly map to the four main components of MAGIC that can be run individually or together:

NIRCam or FGS science to FGS raw image conversion (Image Converter)
------------------------------------------------
This module can take in a (simulated or real) NIRCam image from any of the short- and longwave detectors, or full-frame calibration FGS image from either detector, and convert it to a pseudo-FGS raw (guider 1 or guider 2) image. This pseudo-FGS raw image that is created is not a simulated image, but an expectation of the raw image that the FGS flight software will see. This module works with the [`jwst-fgs-countrate`](https://github.com/spacetelescope/jwst-fgs-countrate) module to use the count rates in the input image to renormalize the image to the expected count rates for a different object. This module appropriately rotates the image from the science frame, adjusts the pixel scale and image size if converting from a NIRCam image, corrects bad pixels, and normalizes the image to the magnitude and count rates of a specified guide star.

  [Section III: Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md)

Star Selection Tool (Star Selector)
-----------------------------------
This module takes in a raw FGS image and allows the user to choose the guide and reference star PSFs using a GUI.

  [Section IV: Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

Flight Software File Writer (Flight Software (FSW) File Writer)
---------------------------------------------------------------
This module requires an FGS raw image (for example, the image that comes out of the Image Converter tool, and accompanying file (for example, the "guiding_selections" file that comes out of the Star Selection Tool) that includes a list of the guide and reference PSF image coordinates and their associated measured count rates. The output is all files necessary to test this image with different flight software simulators (FGS DHAS, FGSES, etc.) includes all the files necessary to run the ID, ACQ1, ACQ2, and TRK steps in these simulators.

  [Section IV: Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

Segment Guiding Tool (Segment Guiding)
--------------------------------------
Used to facilitate guiding on unstacked segment PSF arrays (Segment Override), and stacked, but un-phased segment PSFs (Photometry Override) during commissioning.

**Segment Override**

When provided with:
1. the commanded RA and Dec of a guide star, and
2. the V2/V3 (or x/y) positions of all segments in an array

the segment guiding tool calculates the effective RA and Dec of all segments on the sky. This information can then be used to override the guide star catalog with new RA, Dec, and count rate values for the requested guide and reference star PSFs

  [Section VII: Writing the Segment Override File (SOF)](vii_write_sof.md)

**Photometry Override**

When provided with a count rate factor and a threshold factor, the tool will provide the necessary information to override the guide star catalog with a new count rate for a specific guide star

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
