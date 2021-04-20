#### Previous

###### Section I: [Introduction](i_introduction.md)

###### Section II: [Setting Up MAGIC](ii_setting_up.md)

###### Section III: [Determining and Loading the Input Image](iii_determining_and_loading_the_input_image.md)

###### Section IV: [Selecting Guide & Reference Stars for an Input Image and Writing Out Files](iv_select_stars_and_write_files.md)

###### Section V: DHAS Versions 3.0 and Earlier:[Testing Selections in DHAS Version 3.0 and Ealier](v_testing_in_dhas.md), DHAS Versions 3.99 and Laster:[Testing Selections in DHAS Version 3.99 and Laster](v_testing_in_dhas_v4.md)

###### Section VI: [Contingency: Re-selecting Stars and Re-running DHAS](vi_contingency_reselect_stars.md)

###### Section VII: [Writing the Segment Override File (SOF)](vii_write_sof.md)

-----------------------------------------

VIII.	Writing the Photometry Override File (POF)
================================================
In the case of MIMF where we only need to change the photometry of the guide star (the RA and Dec and expected count rates are taken from the APT file), we need to make an override for Planning & Scheduling but this is for the photometry and will have no information about the segments (because the PSFs from the individual segments are stacked).


Creating a photometry override file through MAGIC:
--------------------------------------------------
1. Load the file for this observation, select the guider, the set the out directory and root.

2. As in section VII, in the main GUI, select the **Segment Guiding** check box. All other options in the interface will be disabled. Note: You do not need to run the other parts of MAGIC when creating a photometry override file.

3. Select the **Create photometry override file** radio button (*B* in figure below)

   ![Segment Guiding Section of the Main GUI](./figs/figure13_segment_guiding.png)

4. Run the tool.

   ![Run MAGIC](./figs/figure_a_run.png)

5. When the Segment Guiding Dialog Box appears (shown in Figure 17), define the segment guiding parameters, including:

   ![Photometry Override Dialog Box](./figs/figure17_photometry_override_dialog.png)

   1. **Program Number** – of the current APT program; three to five digits. Input *only* the Program Number if *all* observations and visits in the program will use the same file.

   2. **Observation Number(s)** (optional) - of the observation(s) that will be executed. To write a file for only one observation, simply write one number (e.g. "13"). To write a file for multiple observations, write a comma-separated list of numbers, and use hyphens to denote ranges (e.g. "1, 3, 5" or "3-7" or "1, 3, 5-7, 9"). Input *only* the Program Number and Observation Number if *all* visits in a given observation will use the same file.

   3. **Visit Number** (optional) – of the visit that will be executed (this is usually 1, but will be different when mosaics, etc. are taken)

   4. **Count rate factor** - A multiplication factor to be applied to the computed count rate of each guide star and reference object of the visit. This factor is used for cases such as MIMF and CP when the segments are stacked but unphased, and so the brightness of the guide star is dimmed. This factor should be greater than 0.0 and less than 1.0.

   5. **Count rate uncertainty factor** - A multiplication factor to be applied to the computed count rate of each guide star and reference object of the visit that will be set at the count rate uncertainties. This factor should be greater than or equal to 0.01 and less than 1.0.

   Note: If you used the APT query functionality in [Section III](iii_determining_and_loading_the_input_image.md), the Program ID, Observation Number, Visit Number, and RA and DEC of the guide star should be pre-populated.

   *See [Appendix C](appendix_c_mirror_states.md) for information about the count rate factor based on the mirror state.*

6. Click **OK**

---------------------------------

#### Next

###### Appendix A: [Installing the JWST MAGIC Package](appendix_a_installing_magic.md)

###### Appendix B: [Setting Up DHAS](appendix_b_opening_dhas.md)

###### Appendix C: [Mirror State Procedures](appendix_c_mirror_states.md)
