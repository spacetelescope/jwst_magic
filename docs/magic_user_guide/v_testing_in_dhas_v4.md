#### Previous

###### Section I: [Introduction](i_introduction.md)

###### Section II: [Setting Up MAGIC](ii_setting_up.md)

###### Section III: [Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md)

###### Section IV: [Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

-----------------------------------------

V.	Testing Selections in DHAS
==============================
**Note: This documentation is for DHAS *Version 3.99* and later.**

After selecting the guide and reference stars, we have to determine whether or not this selection will be successful with the flight software (FSW). To do this, we can use the files created in the last section to run our images through the FSW simulator, the FGS Data Handling Analysis System (DHAS) (see figure below). (If you do not have the DHAS set up, see [Appendix B](appendix_b_opening_dhas.md))

*Don’t close the MAGIC GUI while you are working in DHAS; you will need to use that open MAGIC window later.*

![DHAS GUI](./figs/figure10_dhas_gui_v4.png)

Note: If you are having issues with how the DHAS GUI is displayed, check out our [DHAS Known Issues page](https://jwstitarwiki.stsci.edu/display/WFSCOWG/Known+Issues+and+Work+Arounds#KnownIssuesandWorkArounds-DHASKnownIssues) for a possible solution.

This section includes testing ID, ACQ, and TRK images with the DHAS. Make sure you know which steps you are testing before continuing, and then use the links below:

1. [Testing ID](#testing-id-in-dhas)

2. [Testing ACQ](#testing-acq-in-dhas)

3. [Testing TRK](#testing-trk-in-dhas)

-----

## Testing ID in DHAS

1. Load the ``{root}_IDstrips.fits`` file for the current data:

   1. Click the green **Load .FITS Images Files** button

      ![DHAS Load FITS button](./figs/figure_b_dhas_load_fits_v4.png)

   2. Using the *Current Directory* drop-down menu (or by typing the path directly into the text box), navigate to the ``out/{root}`` directory where the images have been saved

   3. Go the ``dhas`` subdirectory

   4. Check the **Show All Files** box

   5. Select the ``{root}_IDstrips.fits`` file

   6. Click **Add →**

   7. Click **Done**


2. Add the Bad Pixel Map

   1. Check the **Apply Bad Pixel Map?** box highlighted in yellow, in the Data Analysis section #ADD_FIGURE. *If the box is checked with a path next to it, a bad pixel map has already applied for the guider used for the last simulation; continue to step 3.*

   2. In the finder window that opens, choose the most recent bad pixel map for the corresponding guider (check the “rev” number on files with a .pix extension). For DHAS version 4, this file is `G1_F11_CV3_Rev_64.pix` for guider 1 and `G2_F14_CV3_Rev_64.pix` for guider 2.

   3. Click **Open**

   4. When the finder box disappears, the **Image Size Mismatch…** dialog box, will appear; this concerns the alignment of the subarrays with bad pixel map and the default values may be accepted by clicking **OK**.


3. Run the ID simulator

   1. Click the corresponding G1 or G2 button depending on if this is a guider 1 vs. guider 2 image at the top of the page

      ![Choosing the Guider in DHAS](./figs/figure_c_dhas_guider_v4.png)

   2. Click thepurple **ID Simulation** button in the Guide Simulations & Telemetry section

      ![Selecting the ID simulator in DHAS](./figs/figure_d_dhas_mode_v4.png)


   3. This will open the idSim_Config window.


      ![idSim_Config Window](./figs/figure_dhas_idSim_Config.png)

   4. To run in single file mode with a .prc file:

      1. Go to last line in the window and click. This will open a finder dialog box. Navigate to where the ``{root}_IDstrips.fits`` you selected is located and choose the associated ``{root}_ID.prc`` file and click Open.

      2. You should see the path to this .prc file in the box where you clicked.

      3. Click the "Write to Star File" button above the box where the path to the .prc file is. This will convert your .prc file into a .star file.

      4. Check the "Star Settings" section to make sure that the selections look correct

      5. If you want to alter your guide & reference star selections, do so here by toggling the **GS** (guide star) and **REF#** (reference star) buttons

    5. To run in batch mode (note that you will need .star files for each configuration)

      1. Click on the box for "Enable Batch Mode"

      2. Click "Select Files" and choose the relevant files ???????

    6. Click the “Run” button.

4. Wait for the simulator to run

5. DHAS results: Commanded (CMD) reference stars are denoted by yellow triangles (![Yellow Triangle](./figs/dhas_commanded_ref_star.png)). The commanded (CMD) guide star is denoted by a yellow cross/plus sign (![Yellow Plus](./figs/dhas_commanded_guide_star.png)). The reference stars that the DHAS has found are denoted by blue x’s (![Blue X](./figs/dhas_found_ref_star.png)) and the guide star is denoted by a blue asterisk (![Blue X](./figs/dhas_found_guide_star.png)). See the figure below for an example of a successful DHAS run.

    ![Example of a successful DHAS finding of guide and reference stars](./figs/figure11_dhas_success.png)

6. Inspect the DHAS results

   1. Do the stars that DHAS found (in blue) match the stars you commanded it to find (in yellow)? If not, ID has FAILED. In the example below, note that DHAS labeled this as a success, even though you can tell in the plot that it failed.

      ![Example of a failed DHAS run](./figs/figure_i_dhas_id_telem_failure.png)

   2. For a more detailed DHAS diagnosis:

      1. Does DHAS think it successfully found the guide star? (Does **Status** equal SUCCESS?)

      2. Did DHAS find all of the stars? (How many **# Candidates** were found (if less than 18 for GA or image array steps, are some of the PSFs difuse?)?)

      3. If necessary, click the **Export** button to more closely examine DHAS’s output plot and/or save the image as a .PNG

        ![Example of a failed DHAS run](./figs/figure_j_dhas_id_telem_success2.png)

If DHAS ID fails, we need to try a different orientation of guide and reference stars until we find a successful one. Continue to [Section VI](vi_contingency_reselect_stars.md).

If DHAS ID succeeds, continue on to test ACQ.

## Testing ACQ in DHAS

1. Load the ``{root}_ACQ.fits`` file you just created:

   1. Click the green **Load .FITS Images Files** button

     ![DHAS Load FITS button](./figs/figure_b_dhas_load_fits_v4.png)

   2. Using the *Current Directory* drop-down menu (or by typing the path directly into the textbox), navigate to the ``out/{root}`` directory where the images have been saved

   3. Go the ``dhas`` subdirectory

   4. Check the **Show All Files** box

   5. Select the ``{root}_ACQ1.fits`` and ``{root}_ACQ2.fits`` files (preserving order, ACQ1 must be the top filename in the box on the right)

   6.	Click **Add →**

   7. Click **Done**

2. Add the Bad Pixel Map *(if the box is checked with a path next to it, a bad pixel map has already applied for the guider used for the last simulation)*. If the box is not already checked, following the steps in step 2 for ID above.

3. Run the ACQ simulator

   1. Click the corresponding G1 or G2 button depending on if this is a guider 1 vs. guider 2 image at the top of the page

      ![Choosing the Guider in DHAS](./figs/figure_c_dhas_guider_v4.png)

   2. Click the purple **ACQ Simulation** button in the Guide Simulations & Telemetry section

      ![Selecting the ACQ simulator in DHAS](./figs/figure_f_dhas_acq_v4.png)

   3. A finder window will appear; select the appropriate `ACQ.prc` file and click **Open** (or double click on the `.prc` file). (Note that Acquisition still uses .prc files)

   4. The **acq_star_catalog_page** dialog box will appear. *Record the Row and Column of the ACQ2 window - you will need this for track*. When you are happy with the values, click **DONE**.

       ![DHAS ACQ Star Catalog](./figs/figure_q_acq_star_catalog.png)

   5. Adjust ACQ parameters:

      **NOTE: You will NOT have to change these parameters if the "Commissioning Parameters" button in the Guide Simulations & Telemetry section is checked. If you button is checked when you run the simulations, skip to step 4.**

      1. First, the **ACQ Mode Setup – ACQ1 Settings** dialog box will appear. You may adjust ACQ1 and ACQ2 table parameters here if you wish.

         For simulations of commissioning data pre-MIMF, change the following values:

         | Parameter Name            | Old Value | New Value |
         |---------------------------|-----------|-----------|
         | `Update default settings?`| 0         | 1         |

         ![DHAS ACQ Mode Setup - ACQ1 Settings Window](./figs/figure_m_dhas_acq1_params.png)

         When you are happy with the values, click **OK**. If you did not change `Update default settings?` to `1`, continue on to step 4.

      2. Second, the **ACQ Mode Setup - ACQ1 Scanning Parameters** dialog box will appear. You may adjust ACQ1 table parameters here if you wish.

         No parameters on this window need to be changed for pre-MIMF commissioning simulations.

         When you are happy with the values, click **OK**.

      3. Third, the **ACQ Mode Setup - ACQ1 Centroid Parameters** dialog box will appear. You may adjust ACQ1 table parameters here if you wish.

         No parameters on this window need to be changed for pre-MIMF commissioning simulations.

         When you are happy with the values, click **OK**.

      4. Fourth, the **ACQ Mode Setup - ACQ2 Parameters** dialog box will appear (see screenshot below). You may adjust ACQ2 table parameters here if you wish.

         For simulations of commissioning data pre-MIMF, make sure the following values are as stated or you will need to update them:

         | Parameter Name            |   Value   |
         |---------------------------|-----------|
         | `softwareSubWindSize`     | 28        |
         | `perimeterSubWindSize`    | 28        |
         | `maxPSFWidth`             | 28        |
         | `maxPSFHeight`            | 28        |
         | `badDeltaSignalThreshold` | 8         |
         | `badDeltaNoiseThreshold`  | 8         |

         When you are happy with the values, click **OK**.

         ![DHAS ACQ Mode Setup - ACQ2 Parameters Window](./figs/figure_p_ACQ_parameters.png)

4.	Inspect the DHAS Results:

    The *ACQ Centroid Plot* and *ACQ Centroid Telemetry* reports will pop up when the simulation is complete.

    In the *ACQ Centroid Telemetry* report table, verify that 1) there are 8 acqSim entries, and 2) the last `IFGS_ACQ_STAT` telemetry value is `SUCCESS`. If either of these is not the case, ACQ has FAILED.

    View the *ACQ Centroid Plot* (see figure below). The acquired guide star locations are denoted with `X` and `+`, and the commanded guide star position is denoted with `*`. Check that 1) the acquired guide star locations are not more than ~10 pixels away from the commanded guide star position and 2) there is not a lot of spread (i.e. more than a couple pixels) within the acquired guide star locations. If either of these is not the case, consult with the FGS SI team about the success of ACQ.

    ![ACQ Centroid Plot from DHAS](./figs/figure_n_ACQ_centroid_plot.png)

If DHAS ACQ fails, we need to try a different orientation of guide and reference stars until we find a successful one. Continue to [Section VI](vi_contingency_reselect_stars.md).

If DHAS ACQ succeeds, continue on to test TRK.

## Testing TRK in DHAS

1. Load the ``{root}_TRK.fits`` file for the current data:

   1. Click the green **Load .FITS Images Files** button

     ![DHAS Load FITS button](./figs/figure_b_dhas_load_fits_v4.png)

   2. Using the *Current Directory* drop-down menu (or by typing the path directly into the text box), navigate to the ``out/{root}`` directory where the images have been saved

   3. Go the ``dhas`` subdirectory

   4. Check the **Show All Files** box

   5. Select the ``{root}_TRK.fits`` file

   6.	Click **Add →**

   7. Click **Done**

2. Run the TRK simulator

   1. Click the corresponding G1 or G2 button depending on if this is a guider 1 vs. guider 2 image at the top of the page

      ![Choosing the Guider in DHAS](./figs/figure_c_dhas_guider_v4.png)

   2. Click the **Stand-Alone** button in the Guider Simulations & Telemetry section, and then the purple **TRACK Simulations** button

      ![Selecting the TRK simulator in DHAS](./figs/figure_g_dhas_trk_v4.png)

   3. The **fd_TFG_sim_config** dialog box will appear (*numbers in the cells highlighted will only need to be updated if the table loads have not applied*):

      ![TRK fd_TFG_sim_config window in DHAS](./figs/figure_l_dhas_trk_window.png)

      Click **fixed window location** and the radio button next to 1024, 1024 in the *Window Location* pane on the right of the dialog box; and replace the values for the row and column with the values you recorded for ACQ2 window above.

      ![TRK Window Location selection with the DHAS](./figs/figure_h_dhas_trk_window_loc.png)

      If the **Commissioning Parameters** check box has not been checked, you can adjust some of the table settings here. For simulations of commissioning data pre-MIMF, make sure that the table parameters listed in the table are as stated (these are the numbers included in the useGAtables table update):

      | Parameter Name            |   Value   |
      |---------------------------|-----------|
      | `perimeterSubWindSize`    | 28        |
      | `softwareSubWindSize`     | 28        |
      | `badDeltaSignalThreshold` | 8         |
      | `badDeltaNoiseThreshold`  | 8         |
      | `maxPSFWidth`             | 28        |
      | `maxPSFHeight`            | 28        |

      *Note that the `subWinUpdateCriteria` and `swWinJump` parameters should ideally also be changed, to 28 and 28 respectively, but those changes would need to be made directly to the `fd_TRK_table.m` file in the DHAS source code. For now, proceed without changing those parameters.*

      When you are happy with the values, click the green **Run** button.

   4. A finder window will appear; select the appropriate bad pixel map file (G1_F11_CV3_Rev_64.pix for guider 1 and G2_F14_CV3_Rev_64.pix for guider 2) and click **Open**.

3. It may take a few seconds for the thermometer-like status bar (see below) to appear and/or to begin to register progress in the simulation. If the status bar ever stalls after initial progress has been made, the simulation is not likely to recover. __*This is considered this a FAILED simulation*__.
       ![TRK Status Bar in DHAS](./figs/figure_k_dhas_trk_status_bar.png)

4. Inspect the DHAS Results:

   Several diagnostic plots will automatically pop up. The same plots are written to the same directory where you got the TRK files as a *.ps* file.

   View the **test_ACQ_G1_TRK: Centroid Location** plot (see example below). Ensure that the NEA (noise equivalent angle) value in pixels, listed underneath the plot title, is below the required threshold.

    | CAR                         | NEA [mas] | NEA [pixels] |
    |-----------------------------|-----------|--------------|
    | OTE-03: SM Focus Sweep      | 1000      | 14.5         |
    | OTE-04: Segment ID          | 1000      | 14.5         |
    | OTE-06: Segment-Image Array | 1000      | 14.5         |
    | OTE-07: Global Alignment    | 100       | 1.45         |
    | OTE-09: Image Stacking      | 50        | 0.72         |
    | OTE-12: Coarse Phasing      | 50        | 0.72         |
    | OTE-18: Fine Phasing        | 4         | 0.058        |
    | OTE-26: WF Monitoring       | 4         | 0.058        |
    | OTE-29: MIMF                | 4         | 0.058        |

   (FGS plate scale: 0.069 arcsec/pixel, 69 mas/pixel)

    *(This table is a work in progress and will eventually hold values for additional CARs)*  


   ![TRK Centroid Location & NEA plot from DHAS](./figs/figure_o_TRK_centroid_plot.png)

   If the NEA is above the threshold, TRK has FAILED.  

If DHAS TRK fails, we need to try a different orientation of guide and reference stars until we find a successful one. Continue to [Section VI](vi_contingency_reselect_stars.md).

If DHAS TRK succeeds, continue on to [Section VII](vii_write_sof.md) to create a segment override file (SOF) or [Section VIII](viii_write_pof.md) to create a photometry override file (POF).

---------------------------------

#### Next

###### Section VI: [Contingency: Re-selecting Stars and Re-running DHAS](vi_contingency_reselect_stars.md)

###### Section VII: [Writing the Segment Override File (SOF)](vii_write_sof.md)

###### Section VIII:  [Writing the Photometry Override File (POF)](viii_write_pof.md)

###### Appendix A: [Installing the JWST MAGIC Package](appendix_a_installing_magic.md)

###### Appendix B: [Setting Up DHAS](appendix_b_opening_dhas.md)

###### Appendix C: [Mirror State Procedures](appendix_c_mirror_states.md)
