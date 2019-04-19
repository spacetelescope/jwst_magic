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
        mkproc.Mkproc(guider, root, xarr, yarr, counts, step)

    Required arguments:
        ``guider`` - guider number (1 or 2)
        ``root`` - name used to create the output directory, {out_dir}/out/{root}
        ``xarr`` - X coordinates of guide and reference stars (pixels)
        ``yarr`` - Y coordinates of guide and reference stars (pixels)
        ``counts`` - count rates of guide and reference stars
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
import shutil

# Local Imports
from jwst_magic.utils import coordinate_transforms

# Paths
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_PATH = os.path.split(__location__)[0]
OUT_PATH = os.path.split(PACKAGE_PATH)[0]  # Location of out/ directory


class Mkproc(object):
    """Makes CECIL proc files for FGS guider 1 and 2
    """

    def __init__(self, guider, root, xarr, yarr, counts, step, thresh_factor=0.5,
                 out_dir=None, acq1_imgsize=None, acq2_imgsize=None):
        """ Initialize the class and create CECIL proc files for guider 1 and 2.

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : list
            X coordinates of guide and reference stars (pixels)
        yarr : list
            Y coordinates of guide and reference stars (pixels)
        counts : list
            Count rates of guide and reference stars
        step : str
            Name of the step to create files for
        thresh_factor : float, optional
            Factor by which to multiply the countrates to determine
            the threshold count rate
        out_dir : str, optional
            Where output files will be saved. If not provided, the
            image(s) will be saved within the repository at
            jwst_magic/
        acq1_imgsize : int, optional
            Dimension of ACQ1 images
        acq2_imgsize : int, optional
            Dimension of ACQ2 images
        """

        # Create output directory if does not exist
        if out_dir is None:
            self.out_dir = os.path.join(OUT_PATH, 'out', root)
        else:
            self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # Find templates. If template path not given, script will assume that a
        # 'templates' directory that includes are necessary prc templates lives
        # in the same directory as this script
        template_path = os.path.join(PACKAGE_PATH, 'data', 'templates')

        self.find_templates(guider, step=step, template_path=template_path)

        # Depending on the 'step' create the correct CECIL proc files.
        if step == 'ID':
            self.create_id_proc_file(guider, root, xarr, yarr, counts,
                                     thresh_factor=thresh_factor)
        elif step == 'ACQ':
            self.create_acq_proc_file(guider, root, xarr, yarr, counts,
                                      acq1_imgsize, acq2_imgsize)

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

    def create_id_proc_file(self, guider, root, xarr, yarr, counts, thresh_factor):
        """Creates the CECIL proc file for the identification (ID) step.
        Writes to {out_dir}/out/{root}/dhas/{root}_G{guider}_ID.prc

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : list
            X coordinates of guide and reference stars (pixels)
        yarr : list
            Y coordinates of guide and reference stars (pixels)
        counts : list
            Count rates of guide and reference stars
        thresh_factor : float
            Factor by which to multiply the countrates to determine
            the threshold count rate
        """
        eol = '\n'
        nref = len(xarr) - 1
        dhas_filename = os.path.join(self.out_dir, 'dhas',
                                     '{0}_G{1}_ID.prc'.format(root, guider))

        with open(dhas_filename, 'w') as file_out:
            self.write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ID'.format(root, guider))
            file_out.write(eol)

            self.write_from_template(self.template_a, file_out)

            # Convert real pixel to DHAS ideal angle
            xangle, yangle = coordinate_transforms.Raw2DHAS(xarr, yarr, guider)

            thresh = thresh_factor * counts

            file_out.write('@IFGS_GUIDESTAR {0}, DFT, {1:12.4f}, {2:12.4f}, \
                            {3:12d}, {4:8d}'.format(self.guider,
                                                    xangle[0],
                                                    yangle[0],
                                                    int(counts[0]),
                                                    int(thresh[0])))
            file_out.write(eol)

            self.write_from_template(self.template_b, file_out)

            if nref >= 1:  # ref stars > 0
                file_out.write('@IFGS_REFCOUNT DETECTOR={0}, REFSTARS={1}'.format(self.guider,
                                                                                  nref))
                file_out.write(eol)

                self.write_from_template(self.template_c, file_out)

                for istar in range(1, nref):
                    file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, \
                                   {3:12.6f}, {4:8d}, {5:8d}'.format(self.guider,
                                                                     int(istar),
                                                                     xangle[istar],
                                                                     yangle[istar],
                                                                     int(counts[istar]),
                                                                     int(thresh[istar])))
                    file_out.write(eol)
                    self.write_from_template(self.template_d, file_out)

                # The last reference star ends with a different template so it
                # written outside of the for loop
                file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, {3:12.6f}, \
                               {4:8d},{5:8d}'.format(self.guider,
                                                     int(nref),
                                                     xangle[nref],
                                                     yangle[nref],
                                                     int(counts[nref]),
                                                     int(thresh[nref])))
                file_out.write(eol)

                self.write_from_template(self.template_e, file_out)

            self.write_from_template(self.template_f, file_out)

        file_out.close()
        print("Successfully wrote: {}".format(dhas_filename))
        shutil.copy2(dhas_filename,
                     os.path.join(self.out_dir, 'ground_system'))
        print("Successfully wrote: {}".format(os.path.join(self.out_dir,
                                                           'ground_system',
                                                           '{0}_G{1}_ID.prc'.
                                                           format(root, guider))))

    def create_acq_proc_file(self, guider, root, xarr, yarr, counts,
                             acq1_imgsize, acq2_imgsize):
        """Creates the CECIL proc file for the acquisition (ACQ) steps.
        Writes to {out_dir}/out/{root}/dhas/{root}_G{guider}_ACQ.prc

        Parameters
        ----------
        guider : int
            Guider number (1 or 2)
        root : str
            Name used to create the output directory, {out_dir}/out/{root}
        xarr : list
            X coordinates of guide and reference stars (pixels)
        yarr : list
            Y coordinates of guide and reference stars (pixels)
        counts : list
            Count rates of guide and reference stars
        acq1_imgsize : int
            Dimension of ACQ1 images
        acq2_imgsize : int
            Dimension of ACQ2 images
        """
        eol = '\n'

        ### THIS IS BAD -- CHANGE THIS --- KJB 1/25/2018
        # Corner coordinates & guide star counts
        xarr1, yarr1 = coordinate_transforms.Raw2DHAS(xarr - acq1_imgsize // 2,
                                                      yarr - acq1_imgsize // 2,
                                                      guider)
        xarr2, yarr2 = coordinate_transforms.Raw2DHAS(xarr - acq2_imgsize // 2,
                                                      yarr - acq2_imgsize // 2,
                                                      guider)

        # Get threshold from countrate (not from STC file)
        threshgs = 0.50 * counts

        dhas_filename = os.path.join(self.out_dir, 'dhas',
                                     '{0}_G{1}_ACQ.prc'.format(root, guider))

        with open(dhas_filename, 'w') as file_out:
            self.write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ACQ'.format(root, guider))
            file_out.write(eol)

            self.write_from_template(self.template_a, file_out)

            # Write guide star coordinates
            file_out.write('@IFGS_GUIDESTAR {0}, 2, {1:12.4f}, {2:12.4f}, \
                           {3:12.4f}, {4:8d}'.format(self.guider, float(xarr1),
                                                     float(yarr1), float(counts),
                                                     int(threshgs)))

            self.write_from_template(self.template_b, file_out)

            # Write ACQ1 box specs
            file_out.write('@IFGS_CONFIG {0}, SWADDRESS=spaceWireAddr1, SLOT=1, NINTS=1, \
            NGROUPS=groupNum1, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=128, NCOLS=128, \
            ROWCORNER={1:12.4f},COLCORNER={2:12.4f}'.format(self.guider,
                                                            float(xarr1),
                                                            float(yarr1)))
            file_out.write(eol)

            self.write_from_template(self.template_c, file_out)

            # Write ACQ2 box specs
            file_out.write('@IFGS_CONFIG {0}, spaceWireAddr2, SLOT=2, NINTS=1, \
            NGROUPS=groupNum2, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=32, NCOLS=32, \
            ROWCORNER={1:12.4f}, COLCORNER={2:12.4f}'.format(self.guider,
                                                             float(xarr2),
                                                             float(yarr2)))
            file_out.write(eol)

            self.write_from_template(self.template_d, file_out)
        file_out.close()
        print("Successfully wrote: {}".format(dhas_filename))
        shutil.copy2(dhas_filename,
                     os.path.join(self.out_dir, 'ground_system'))
        print("Successfully wrote: {}".format(os.path.join(
            self.out_dir, 'ground_system', '{0}_G{1}_ACQ.prc'.format(root,
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
