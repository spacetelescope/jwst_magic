IV.	Selecting Guide & Reference Stars for an Input Image and Writing Out Files
==============================================================================
One of the main features of MAGIC is that it allows the user to determine the guide and reference stars for a specific scene. While during normal opterations this is determined by the Guide Star Selection System (GSSS), during commissioning the MAGIC user will determine which PSFs will be used for guiding and as reference stars. You can turn on this feature by selecting the “Star Selector” check box. 

1. Set star selection parameters:

    [figure 7]

    Figure 7 - Star Selection section of the Main GUI

   a. Ensure the “Star Selector” box is checked.
   b. Inspect the input image and **indicate if the PSFs are non-standard**. (*A*) This flag alters the PSF-finding algorithm in the star selector tool to widen the smoothing filter for diffuse images in early commissioning stages when the telescope is unphased. If you are unsure if the PSFs are phased, consult the “Guiding Method” row in the Guider Commissioning Summary Table (https://innerspace.stsci.edu/display/INSTEL/Guider+Commissioning+Summary+Table) on Innerspace.
   c. If desired, **load pre-selected guide and reference stars from a file** (*B*) by selecting the “Load from File” option and selecting the desired input file. This file must include X/Y pixel coordinates and count rates in the form of a filepath to a regfile.txt or .incat file. Providing this will bypass using the Star Selection GUI to click-to-select the guide and reference stars. 

2. Set file writer parameters:

    [figure 8]
  
    Figure 8 - Flight Software file writer section for the Main GUI

   a. Ensure the “Flight Software (FSW) File Writer” box is checked.
   b. Check that all of the *necessary FGS steps* are selected. 
     a. For general guiding, this includes all of the operational steps: ID, ACQ, and TRK. (These are the default selections.) 
     b. For calibration observations, add the CAL step.
   c. If you want to shift your image so that the selected guide star is moved to the center of the image, ensure the “Place the guide star at the ID attitude” box is checked. Designate whether the guiding field is crowded enough that the alternate ID attitude at (Ideal X, Ideal Y) = (-45, 0) should be used (“Crowded field”). Otherwise, leave the “Nominal” button selected such that the star is placed at (Ideal X, Ideal Y) = (0, 0).

3. Run the tool

  [figure a - run button]

4. Monitor the terminal window from which you launched the GUI to notice any possible errors that are raised. 

  Note:	The output that appears in the command line is also written to::
     
     /data/jwst/wss/guiding/MAGIC_logs/

5. When the Star Selection GUI appears: 

    [figure 9]
  
    Figure 9 - Star Selection GUI window

   a. Inspect the PSFs in the image by moving your cursor over different PSFs. Examine the profile plot to see the distribution of light.
   b. Select, by clicking, which PSFs will be the guide star and the reference stars. The first star selected will be the guide star, while any subsequent stars will be reference stars. *See Appendix D to choose the guide and reference stars based on the mirror state.*
   c. If you want to change your selections while in the tool, use the **“Make Guide Star”** (*A*) button to change the guide star, use the **“Delete”** button (*B*) to remove individual selections, and use the **“Clear Selections”** button (*C*) to start over.
   d. When you are happy with your selections, click “Done”  
   
   The output files will be located in the specified out directory.

