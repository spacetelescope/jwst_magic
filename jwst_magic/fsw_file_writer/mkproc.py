"""Create CECIL proc files needed to simulate ID and ACQ with the DHAS.

Authors
-------
    - Keira Brooks
    - Lauren Chambers

Use
---
    This module can be used as such:
    ::
        from jwst_magic.fsw_file_writer import mkproc
        mkproc.Mkproc(guider, root, xarr, yarr, count_rates, step)

    Required arguments:
        ``guider`` - guider number (1 or 2)
        ``root`` - name used to create the output directory, {out_dir}/out/{root}
        ``xarr`` - X coordinates of guide and reference stars (pixels)
        ``yarr`` - Y coordinates of guide and reference stars (pixels)
        ``count_rates`` - count rates of guide and reference stars
        ``step`` - name of the step to create files for

    Optional arguments:
        ``thresh_factor`` - factor by which to multiply the countrates
            to determine the threshold count rate
        ``out_dir`` - where output files will be saved. If not provided,
            the image(s) will be saved within the repository at
            jwst_magic/
        ``acq1_imgsize`` - dimension of ACQ1 images
        ``acq2_imgsize`` - dimension of ACQ2 images

Notes
-----
    This code is adapted from mkIDproc.pro/mkACQproc.pro from Sherie Holfeltz.
"""

# Standard Library Imports
import os
import logging
import shutil

import numpy as np

# Local Imports
from jwst_magic.utils import coordinate_transforms

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ directory

# Start logger
LOGGER = logging.getLogger(__name__)

class Mkproc(object):
    """Makes CECIL proc files for FGS guider 1 and 2
    """

    def __init__(self, guider, root, xarr, yarr, count_rates, step, threshold=None,
                 out_dir=None, dhas_dir='dhas', ground_system_dir='ground_system',
                 acq1_imgsize=None, acq2_imgsize=None):
        """ Initialize the class and create CECIL proc files for guider 1 and 2.

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : array
            X coordinates of guide and reference stars (pixels)
        yarr : array
            Y coordinates of guide and reference stars (pixels)
        count_rates : list
            Count rates of guide and reference stars
        step : str
            Name of the step to create files for
        threshold : float or list, optional
            Absolute threshold value(s) to go into the prc file. Must be a list
            of the same length as the number of PSFs in the config
        out_dir : str, optional
            Where output files will be saved. If not provided, the
            image(s) will be saved within the repository at
            jwst_magic/
        dhas_dir : str
            Name of dhas directory. Either 'dhas' or 'dhas_shifted'. Only written to for ACQ prc
        ground_system_dir : str
            Name of ground_system directory. Either 'ground_system' or 'ground_system_shifted'
        acq1_imgsize : int
            Array size of the ACQ1 window
        acq2_imgsize : int
            Array size of the ACQ2 window
        """

        # Create output directory if does not exist
        if out_dir is None:
            self.out_dir = os.path.join(OUT_PATH, 'out', root)
        else:
            self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        self.dhas_dir, self.ground_system_dir = dhas_dir, ground_system_dir
        for dir in [dhas_dir, ground_system_dir]:
            if not os.path.exists(os.path.join(self.out_dir, dir)):
                os.makedirs(os.path.join(self.out_dir, dir))

        # Find templates. If template path not given, script will assume that a
        # 'templates' directory that includes are necessary prc templates lives
        # in the same directory as this script
        # TODO: be smarter
        data_path_central_store = '/itar/jwst/tel/share/wf_guiding/magic_data_files'
        template_path = os.path.join(data_path_central_store, 'templates')
        if not os.path.exists(template_path):
            template_path = os.path.join(PACKAGE_PATH, 'data', 'templates')
            LOGGER.warning('Make Proc: Cannot access central store. No prc can be made unless access to central store can be restored.')

        self.find_templates(guider, step=step, template_path=template_path)

        # Confirm the threshold is a list
        if not isinstance(threshold, (list, np.ndarray)):
            threshold = [threshold]

        # Depending on the 'step' create the correct CECIL proc files.
        if step == 'ID':
            self.create_id_proc_file(guider, root, xarr, yarr, count_rates, threshold=threshold)
        elif step == 'ACQ':
            self.create_acq_proc_file(guider, root, xarr, yarr, count_rates, threshold=threshold,
                                      acq1_imgsize=acq1_imgsize, acq2_imgsize=acq2_imgsize)

    def find_templates(self, guider, step, template_path):
        """Open the different templates used to make the proc file.

        Parameters
        ----------
        guider : int
            Guider number (1 or 20)
        step : str
            Name of step ('ID' or 'ACQ')
        template_path : str
            Path to prc templates
        """
        guider = str(guider)
        self.guider = 'GUIDER{}'.format(guider)
        self.template_hdr = os.path.join(template_path, 'g{}{}templateHDR.prc'.format(guider, step))
        self.template_a = os.path.join(template_path, 'g{}{}templateA.prc'.format(guider, step))
        self.template_b = os.path.join(template_path, 'g{}{}templateB.prc'.format(guider, step))
        self.template_c = os.path.join(template_path, 'g{}{}templateC.prc'.format(guider, step))
        self.template_d = os.path.join(template_path, 'g{}{}templateD.prc'.format(guider, step))
        if step == 'ID':
            self.template_e = os.path.join(template_path, 'g{}{}templateE.prc'.format(guider, step))
            self.template_f = os.path.join(template_path, 'g{}{}templateF.prc'.format(guider, step))

    def create_id_proc_file(self, guider, root, xarr, yarr, count_rates, threshold=None):
        """Creates the CECIL proc file for the identification (ID) step.
        Writes to {out_dir}/out/{root}/dhas/{root}_G{guider}_ID.prc

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : array
            X coordinates of guide and reference stars (pixels)
        yarr : array
            Y coordinates of guide and reference stars (pixels)
        count_rates : list
            Count rates of guide and reference stars
        threshold : float, optional
            Absolute threshold value.
        """
        eol = '\n'
        nref = len(xarr) - 1
        ground_system_filename = os.path.join(self.out_dir, self.ground_system_dir,
                                              '{0}_G{1}_ID.prc'.format(root, guider))

        with open(ground_system_filename, 'w') as file_out:
            self.write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ID'.format(root, guider))
            file_out.write(eol)

            self.write_from_template(self.template_a, file_out)

            # Convert real pixel to DHAS ideal angle
            xangle, yangle = coordinate_transforms.raw2dhas(xarr, yarr, guider)

            file_out.write('@IFGS_GUIDESTAR {0}, DFT, {1:12.4f}, {2:12.4f}, \
                            {3:12.4f}, {4:12.4f}'.format(self.guider,
                                                         xangle[0],
                                                         yangle[0],
                                                         float(count_rates[0]),
                                                         float(threshold[0])))
            file_out.write(eol)

            self.write_from_template(self.template_b, file_out)

            if nref >= 1:  # ref stars > 0
                file_out.write('@IFGS_REFCOUNT DETECTOR={0}, REFSTARS={1}'.format(self.guider,
                                                                                  nref))
                file_out.write(eol)

                self.write_from_template(self.template_c, file_out)

                for istar in range(1, nref):
                    file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, \
                                   {3:12.6f}, {4:12.4f}, {5:12.4f}'.format(self.guider,
                                                                           int(istar),
                                                                           xangle[istar],
                                                                           yangle[istar],
                                                                           float(count_rates[istar]),
                                                                           float(threshold[istar])))
                    file_out.write(eol)
                    self.write_from_template(self.template_d, file_out)

                # The last reference star ends with a different template so it
                # written outside of the for loop
                file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, {3:12.6f}, \
                               {4:12.4f},{5:12.4f}'.format(self.guider,
                                                           int(nref),
                                                           xangle[nref],
                                                           yangle[nref],
                                                           float(count_rates[nref]),
                                                           float(threshold[nref])))
                file_out.write(eol)

                self.write_from_template(self.template_e, file_out)

            self.write_from_template(self.template_f, file_out)

        file_out.close()
        LOGGER.info("Successfully wrote: {}".format(ground_system_filename))

    def create_acq_proc_file(self, guider, root, xarr, yarr, count_rates,
                             acq1_imgsize, acq2_imgsize, threshold=None):
        """Creates the CECIL proc file for the acquisition (ACQ) steps.
        Writes to {out_dir}/out/{root}/dhas/{root}_G{guider}_ACQ.prc

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : array
            X coordinates of guide and reference stars (pixels)
        yarr : array
            Y coordinates of guide and reference stars (pixels)
        count_rates : list
            Count rates of guide and reference stars
        acq1_imgsize : int
            Array size of the ACQ1 window
        acq2_imgsize : int
            Array size of the ACQ2 window
        threshold : float, optional
            Absolute threshold value.
        """
        eol = '\n'

        if len(xarr) != 1:
            xarr, yarr = xarr[0], yarr[0]

        # Corner coordinates & guide star count rate
        gs_xangle, gs_yangle = coordinate_transforms.raw2dhas(xarr, yarr, guider)
        a1_xangle, a1_yangle = coordinate_transforms.raw2dhas(int(xarr - acq1_imgsize / 2),
                                                              int(yarr - acq1_imgsize / 2), guider)
        a2_xangle, a2_yangle = coordinate_transforms.raw2dhas(int(xarr - acq2_imgsize / 2),
                                                              int(yarr - acq2_imgsize / 2), guider)

        # Get threshold from countrate (not from STC file)
        if len(count_rates) != 1:
            count_rates = count_rates[0]
        threshold = threshold[0]

        dhas_filename = os.path.join(self.out_dir, self.dhas_dir,
                                     '{0}_G{1}_ACQ.prc'.format(root, guider))

        with open(dhas_filename, 'w') as file_out:
            self.write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ACQ'.format(root, guider))
            file_out.write(eol)

            self.write_from_template(self.template_a, file_out)

            # Write guide star coordinates
            file_out.write('@IFGS_GUIDESTAR {0}, 2, {1:12.4f}, {2:12.4f}, \
                           {3:12.4f}, {4:12.4f}'.format(self.guider, float(gs_xangle),
                                                        float(gs_yangle), float(count_rates),
                                                        float(threshold)))

            self.write_from_template(self.template_b, file_out)

            # Write ACQ1 box specs
            file_out.write('@IFGS_CONFIG {0}, SWADDRESS=spaceWireAddr1, SLOT=1, NINTS=1, \
            NGROUPS=groupNum1, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=128, NCOLS=128, \
            ROWCORNER={1:12.4f},COLCORNER={2:12.4f}'.format(self.guider,
                                                            float(a1_xangle),
                                                            float(a1_yangle)))
            file_out.write(eol)

            self.write_from_template(self.template_c, file_out)

            # Write ACQ2 box specs
            file_out.write('@IFGS_CONFIG {0}, spaceWireAddr2, SLOT=2, NINTS=1, \
            NGROUPS=groupNum2, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=32, NCOLS=32, \
            ROWCORNER={1:12.4f}, COLCORNER={2:12.4f}'.format(self.guider,
                                                             float(a2_xangle),
                                                             float(a2_yangle)))
            file_out.write(eol)

            self.write_from_template(self.template_d, file_out)
        file_out.close()
        LOGGER.info("Successfully wrote: {}".format(dhas_filename))
        shutil.copy2(dhas_filename,
                     os.path.join(self.out_dir, self.ground_system_dir))
        LOGGER.info("Successfully wrote: {}".format(os.path.join(
            self.out_dir, self.ground_system_dir, '{0}_G{1}_ACQ.prc'.format(root,
                                                                            guider)
        )
        ))

    @staticmethod
    def write_from_template(template, file_out):
        """Write the lines from the specified template to the specified file

        Parameters
        ----------
        template : str
            Filepath to prc template
        file_out : str
            Output file
        """
        with open(template, 'r') as temp:
            for line in temp:
                file_out.write(line)
        temp.close()
