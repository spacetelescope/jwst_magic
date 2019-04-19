(To create a PDF of this user's guide, go [here](magic_pdf_disclaimer.md))


![MAGIC Logo](../../magic_logo.png)


# JWST MAGIC: User's Guide

How to use the Multi-Application Guiding Interface for Commissioning (MAGIC) for simulated and real JWST data. 

#### Authors:
* Lauren Chambers (lchambers@stsci.edu) 
* Keira Brooks (kbrooks@stsci.edu)
* Sherie Holfeltz

--------


## Table of Contents

#### I. [Introduction](i_introduction.md) 

   &nbsp;&nbsp;&nbsp;&nbsp; Description of the different tools available in MAGIC
   

#### II. [Setting Up](ii_setting_up.md)

   &nbsp;&nbsp;&nbsp;&nbsp; Pull latest version of the MAGIC repository

#### III. [Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md) 
         
   &nbsp;&nbsp;&nbsp;&nbsp; Load an image that you will use for running the MAGIC tools

#### IV. [Selecting Guide & Reference Stars for an Input Image and Writing Files](iv_select_stars_and_write_files.md) 

   &nbsp;&nbsp;&nbsp;&nbsp; Select the guide and reference stars to be used for the image loaded in Section III

#### V. [Testing Selections in DHAS](v_testing_in_dhas.md) 

   &nbsp;&nbsp;&nbsp;&nbsp; Run ID, ACQ, and TRK images through the DHAS

   1. [Testing ID](v_testing_in_dhas.md#testing-id-in-dhas)

   2. [Testing ACQ](v_testing_in_dhas.md#testing-acq-in-dhas)

   3. [Testing TRK](v_testing_in_dhas.md#testing-trk-in-dhas)
   
#### VI. [*Contingency*: Re-selecting Stars and Re-running DHAS](vi_contingency_reselect_stars.md) 
 
   &nbsp;&nbsp;&nbsp;&nbsp; If the DHAS failed on selected guide and reference stars, re-run the tool with new selections

#### VII. [Writing the Segment Override File](vii_write_sof.md) 

   &nbsp;&nbsp;&nbsp;&nbsp; Create the Segment Override File for situations with more than one PSF for a guide star

#### VIII. [Writing the Photometry Override File](viii_write_pof.md) 

   &nbsp;&nbsp;&nbsp;&nbsp; Create the Photometry Override File for situations with one, un-phased, PSF for a guide star


### Appendices

#### A. [Installing the JWST MAGIC Package](appendix_a_installing_magic.md)

#### B. [Setting up DHAS](appendix_b_opening_dhas.md)

#### C. [Using APT to Get Guide Star RA & Dec](appendix_c_apt.md)

#### D. [Mirror State Procedures](appendix_d_mirror_states.md)
