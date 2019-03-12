III.	Determining and Loading the Input Image
=============================================

MAGIC takes in any FGS image or a NIRCam image that was taken with the CLEAR filter (the NIRCam weak lens (WL) filter will cause MAGIC to crash). If you want MAGIC to convert this image into a raw dectector FGS image, you can indicate this by checking the “Convert Image” check box. In most cases MAGIC can determine from the header information which instrument and detector the image is from, however if it can’t figure this out from header information, you will have to help it out by telling it which instrument and detector your input image comes from. You can also renormalize and/or add background images to your seed image.

1. From your astroconda environment, start an IPython session and launch the main GUI:

       $ ipython
    
       In [1]: import jwst_magic
       In [2]: jwst_magic.run_tool_GUI()

 [figure 1]

 Figure 1 - Main GUI for the JWST MAGIC Tool


2. Set general input parameters:

    [figure 2]
   
    Figure 2 - General Input section of the Main GUI
   
   a. **Load the input image** (*A*) and a preview of the image and the full path to the image will appear in the Image Preview box at right.
   b. **Specify the guider** (*B*) that the final image should simulate. If this is not known, check the APT file (see Appendix B for more information about using APT).
 
3. If you are running MAGIC on the SOGS network to generate files for commissioning:
    
    [figure 3]
    
    Figure 3 – Commissioning naming parameter section of the main GUI 

   a. **Check the “Commissioning” radio button** (*C*) to set the naming method.
   b. **Select the practice name** (*D*) corresponding to the current activity.
   c. **Select the CAR or step name** (*E*) of the activity you are generating an override file for.
   d. **Select the observation number** (*F*) of the activity you are generating an override file for.
   
   
   Considering these parameters all together, the output files will be saved in the ``/data/jwst/jwss/guiding/{practice}/{car}/out/for_obs{obs}/`` directory, with the root ``for_obs{obs}_G{guider}``.


4. If you are running MAGIC off of SOGS, or to generate test data:

    [figure 4]
    
    Figure 4 – Manual naming parameter section of the main GUI
   
   a. **Check the “Manual” radio button** (*G*) to set the naming method.
   b. **Specify a root name** (*H*) If different than the default name that was created when the input image was uploaded. The root will be used to to create the output directory where all created files will reside, out/{root}.
   c. **Change the out directory** (*I*) Choose the location to where the files were copied in Part II. An out/ directory will be created in this location, and this is where all the files will be saved.

   
   Considering these parameters all together, the output files will be saved in the ``{out}/out/{root}/`` directory, with names of the format ``{root}_G{guider}``
   
5. Set image conversion parameters: (Note: The steps labelled “optional” below will create higher-fidelity simulations, but are not necessary when using MAGIC to generate FSW input or segment override files.)

    [figure 5] 
   
    Figure 5 - Image Converter section of the Main GUI
   
   a. (Optional) **Simulate the effects of coarse pointing** (*A*)  by specifying the jitter rate of the observatory. A jitter rate of 0.7 arcsec/sec creates images that are similar to ITM simulations in coarse point. Otherwise, ensure the “Add jitter rate” box is unchecked.
   b. **Check that the instrument** (*B*) **and NIRCam detector** (*C*) used to take the input image are set to the correct values; change them if not. (If the NIRCam detector is not defined, the tool will attempt to parse it from the input FITS header.) The FGS-formatted image will be saved to ``out/{root}/FGS_imgs/{input_image}_G{guider}.fits``
   c. (Optional) **Specify the magnitude or counts for the normalization** (*D*) of the final image. Otherwise, ensure the “Normalize to” box is unchecked.
   d. (Optional) **Add background stars** to the final image.
       
      i. Click “Add Background Stars”. (*E*) The background stars dialog box will appear:
      
      [figure6]
         
      Figure 6 - Background stars dialog window
          
      ii. Select which method you wish to use to add stars to the image: randomly, with a user-defined table, or with a Guide Star Catalog (GSC) 2.4.1 query.
          
          1. To add stars randomly:
             
             a. Select the **“Add Stars Randomly”** (*A*) checkbox.
             b. Input the number of stars you want to add to the image
             c. Specify the magnitude range that these additional stars will lie between (relative to the magnitude of the guide star)
          
          2. To add stars individually:
             
             a. Select the **“Define Stars to Add”** (*B*)  checkbox.
             b. If you wish to load star locations and brightness from a file, indicate the location of that file.
             c. Otherwise, enter into the table the X position in pixels, the Y position in pixels, and the countrate in J Magnitude of each star you wish to add. Click the “Add Another Star” button to add another row to the table, or the “Delete Star” button to remove a row.
          
          3. To add stars using a web query from the Guide Star Catalog:
             
             a. Select the **“Query Stars from Guide Star Catalog 2.4.1”** (*C*) checkbox.
             b. Enter the RA and Dec of the guide star, being sure to specify if the RA units as either hours or degrees.
             c. Enter the position angle (roll angle) of the observatory.
             d. Click the “Query GSC” button to add the stars that are visible in the FOV of the selected guider.
      iii. Click “Done” to save and apply these selections, or click “Cancel” to close the window without updating the background star selections.
      iv. Verify that the indicator shows that thcorrect number of background stars have been added.





   


