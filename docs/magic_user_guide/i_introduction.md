I. Introduction
===============

Guiding During the JWST Optical Telescope Element commissioning
---------------------------------------------------------------
### Why we need to interfere with nominal guiding during JWST Optical Telescope Element commissioning
The commissioning of the JWST Optical Telescope Element (OTE) will take the PSFs that come off every primary mirror segment and stack and phase them. At this point we expect 18 PSFs for every star in the field of view (FOV) until we can identify them all, and over the course of several months, stack and phase those PSFs.

The nominal guiding scenario retrieves the location, count rate, and count rate threshold of a guide star from the guide star catalog (GSC) and matches that to what is seen on the FGS detector using the FGS flight software. During the JWST OTE commissioning, we will need to start guiding so we can utilize the resulting stability to get better measurements during this stacking and phasing process. However, with 18 PSFs for each star in the FOV with vastly different count rates and locations from the values in the GSC, the nominal guiding strategy will not work. This resulted in the need for a tool, or set of tools, that will allow us to implement guiding during the JWST OTE commissioning by overriding the GSC with adjusted locations, count rates, and count rate thresholds for the guide and reference stars. These tools will have to work for the case of un-stacked and un-phased PSFs, as well as stacked, but still un-phased PSFs.

### The challenges
One of the biggest challenges of guiding during OTE commissioning, is anticipating what the FGS flight software will encounter when it attempts to identify and then closed-loop guide on a non-nominal PSF. In addition to there being 18x more PSFs in the FOV, we expect additional aberrations in the PSFs as a result of small changes in the actuators controlling each segment. Essentially, each PSF created by a specific mirror segment (segment PSFs) will likely have additional and different aberrations that will be corrected over the course this phase of JWST commissioning, but will also add to the challenge of closed-loop guiding. Instead of using models and simulated images to anticipate the PSFs, particularly when there are so many unknowns, the most dependable source of information about the locations and count rates of the PSFs would be the last image taken by the telescope. This gives a sense of the PSF configuration (locations and shapes of the PSFs) which we can then use to anticipate what the FGS flight software will see when it goes to identify the guide star. By intercepting the flow of information at this point, we can create a modified mini-GSC that can be used to override the main GSC.

During OTE commissioning, the main instrument used is NIRCam. This means that the images we are using for input will have the following differences from what we would expect from an FGS image: FOV, pixelscale, wavelength coverage, detector effects, etc. Additionally, the guiding modes for FGS (ID, ACQ, TRK, and FG) have varying image sizes and exposure times, none of which match the NIRCam full frame images. And, while the NIRCam images we get during commissioning have been processed, since we have to anticipate what the FGS flight software will see, we need to have an idea what the raw frame image looks like.

Lastly, the target star in the NIRCam detector in the input image will likely never be used as the guide star in the FGS detector, so it is not only necessary to update the pixelscale, move from the science frame to FGS raw frame, and correct for any detector effects that are specific to each instrument and detector, but also account for a different count rate in each PSF.

While overcoming these challenges seems like a lot of work to override the GSC, it is still likely more effective than trying to simulate what we will see.


#### Moving from a NIRCam image to an FGS image: A Note on Filters, jitter, and Exposure Time
The NIRCam full frame shortwave detector images have an exposure times of ~10s. These are the images that will be used to determine the necessary locations and count rate values that will be used to override the GSC. For the identification (ID) mode of guiding, 36 strips of 2048 X 64 pixels, with some overlap between strips are taken with a much shorter exposure time of ~1/3s. This would mean that the amount of PSF smearing due to jitter will be much more in the NIRCam images and therefore would need to be accounted for. That being said, the NIRCam image is taken with a filter, while FGS images are taken with the full wavelength coverage since there are no filters on the FGS instrument. This wavelength coverage difference manifests itself in what appears to be a blurring or smearing of each PSF. Therefore, we make the rough assumption that additional jitter due to longer exposure times in the full frame NIRCam image, is roughly equivalent to smearing due to additional wavelengths observed in the FGS ID, ACQ, and TRK images.


### Determining what to guide on
When the segments are un-stacked, since there will be multiple PSFs per star in the FOV, the catalog guide and reference stars' commanded positions will not make sense. Instead, it will be necessary to select a guide and reference star from the guide star's segment PSFs. An isolation criteria will be imposed on guide star candidates during commissioning so that the guide star will be the only "bright" candidate in the FOV. When the segments are stacked, but not phased, the catalog locations for the guide and reference stars can be used since the catalog RA and Dec values are valid.

### Testing
Once the guide and reference stars have been selected, simulated images of the ID, ACQ, and TRK guiding modes can easily be created that can then be run through a flight software simulator like the FGS DHAS. This will allow us to test the success of the current configuration.

### Overriding the guide star catalog
Once the configuration has been verified, the values calculated for the count rate, count rate threshold, and/or locations of the guide and reference stars, can be used to create a command to override the GSC.


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
