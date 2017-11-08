'''This code is adapted from mkIDproc.pro/mkACQproc.pro from Sherie Hofeltz.'''

#STDLIB
import os
import shutil

# Third Party
import numpy as np

# LOCAL
import coordinate_transforms

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
class Mkproc(object):
    '''
    Makes CECIL proc files for FGS guider 1 and 2
    '''
    def __init__(self, guider, root, xarr, yarr, counts, step, thresh_factor=0.5,
                 out_dir=None):
        '''
        Parameters
        ==========
        guider: int
            Guider "1" or "2"
        root: str
            Root for data
        reg_file: str
            <location>/<root>_regfile.txt (Include y, x, and counts for
            each star)
        step: str
            'ID' or 'ACQ'
        thresh_factor: float
            Threshhold factor to be multiplied to the counts. Default = 0.5
        '''


        # Create output directory if does not exist
        if out_dir is None:
            self.out_dir = os.path.join(LOCAL_PATH, 'out', root)
        else:
            self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # Find templates. If template path not given, script will assume that a
        #'templates' directory that includes are necessary prc templates lives
        #in the same directory as this script
        template_path = os.path.join(LOCAL_PATH, 'templates')

        self.find_templates(guider, step=step, template_path=template_path)

        # Depending on the 'step' create the correct CECIL proc files.
        if step == 'ID':
            self.create_id_proc_file(root, guider, xarr, yarr, counts,
                                     thresh_factor=thresh_factor)
        elif step == 'ACQ':
            self.create_acq_proc_file(root, guider)

    def find_templates(self, guider, step, template_path):
        '''
        Pull out the different templates used to make the proc file.

        Parameters
        ==========
        guider: int
            Guider '1' or '2'

        step: str
            'ID' or 'ACQ'

        template_path: str
            Path to the templates
        '''
        self.guider = 'GUIDER{}'.format(guider)
        self.template_hdr = os.path.join(template_path, 'g{}{}templateHDR.prc'.format(guider, step))
        self.template_a = os.path.join(template_path, 'g{}{}templateA.prc'.format(guider, step))
        self.template_b = os.path.join(template_path, 'g{}{}templateB.prc'.format(guider, step))
        self.template_c = os.path.join(template_path, 'g{}{}templateC.prc'.format(guider, step))
        self.template_d = os.path.join(template_path, 'g{}{}templateD.prc'.format(guider, step))
        if step == 'ID':
            self.template_e = os.path.join(template_path, 'g{}{}templateE.prc'.format(guider, step))
            self.template_f = os.path.join(template_path, 'g{}{}templateF.prc'.format(guider, step))

    def create_id_proc_file(self, root, guider, xarr, yarr, counts, thresh_factor):
        '''
        Create the CECIL proc file for the identification phase
        '''
        eol = '\n'
        nref = len(xarr)-1
        print('Number of reference stars: {}'.format(nref))
        #sgid='_G{}'.format(guider)

        with open(os.path.join(self.out_dir, 'dhas',
                               '{0}_G{1}_ID.prc'.format(root, guider)), 'w') as file_out:
            write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ID'.format(root, guider))
            file_out.write(eol)

            write_from_template(self.template_a, file_out)

            #y,x,c = (np.loadtxt(reg_file, delimiter=' ', skiprows=1)).T  # star coords & gs counts

            # Convert real pixel to DHAS ideal angle
            xangle, yangle = coordinate_transforms.rptoia(xarr, yarr, guider)
            xangle, yangle = coordinate_transforms.iatoDHAS(xangle, yangle, guider)

            thresh = thresh_factor * counts

            file_out.write('@IFGS_GUIDESTAR {0}, DFT, {1:12.4f}, {2:12.4f}, \
                            {3:12d}, {4:8d}'.format(self.guider,
                                                    xangle[0],
                                                    yangle[0],
                                                    int(counts[0]),
                                                    int(thresh[0])))
            file_out.write(eol)

            write_from_template(self.template_b, file_out)

            if nref >= 1: # ref stars > 0
                file_out.write('@IFGS_REFCOUNT DETECTOR={0}, REFSTARS={1}'.format(self.guider,
                                                                                  nref))
                file_out.write(eol)

                write_from_template(self.template_c, file_out)

                for istar in range(1, nref):
                    file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, \
                                   {3:12.6f}, {4:8d}, {5:8d}'.format(self.guider,
                                                                     int(istar),
                                                                     xangle[istar],
                                                                     yangle[istar],
                                                                     int(counts[istar]),
                                                                     int(thresh[istar])))
                    file_out.write(eol)
                    write_from_template(self.template_d, file_out)

                #The last reference star ends with a different template so it
                #written outside of the for loop
                file_out.write('@IFGS_REFSTAR {0}, {1:5d}, {2:12.6f}, {3:12.6f}, \
                               {4:8d},{5:8d}'.format(self.guider,
                                                     int(nref),
                                                     xangle[nref],
                                                     yangle[nref],
                                                     int(counts[nref]),
                                                     int(thresh[nref])))
                file_out.write(eol)

                write_from_template(self.template_e, file_out)

            write_from_template(self.template_f, file_out)

        file_out.close()
        print("Successfully wrote: {0}_G{1}_ID.prc".format(root, guider))
        shutil.copy2(os.path.join(self.out_dir, 'dhas', '{0}_G{1}_ID.prc'.format(root,
                                                                                 guider)),
                     os.path.join(self.out_dir, 'ground_system'))


    def create_acq_proc_file(self, root, guider):
        '''
        Create the CECIL proc file for the acquisition phase
        '''
        eol = '\n'

        #corner coords & gs counts
        ind1, xarr1, yarr1, counts1 = np.loadtxt(os.path.join(self.out_dir, 'stsci',
                                                              '{}_G{}_ACQ1.stc'.format(root,
                                                                                       guider))).T
        threshgs = 0.50 * counts1

        ind2, xarr2, yarr2, counts2 = np.loadtxt(os.path.join(self.out_dir, 'stsci',
                                                              '{}_G{}_ACQ2.stc'.format(root,
                                                                                       guider))).T

        with open(os.path.join(self.out_dir, 'dhas',
                               '{0}_G{1}_ACQ.prc'.format(root, guider)), 'w') as file_out:
            write_from_template(self.template_hdr, file_out)

            file_out.write('PROC {0}_G{1}_ACQ'.format(root, guider))
            file_out.write(eol)

            write_from_template(self.template_a, file_out)


            file_out.write('@IFGS_GUIDESTAR {0}, 2, {1:12.4f}, {2:12.4f}, \
                           {3:12.4f}, {4:8d}'.format(self.guider, float(xarr1),
                                                     float(yarr1), float(counts1),
                                                     int(threshgs)))
            write_from_template(self.template_b, file_out)
            file_out.write('@IFGS_CONFIG {0}, SWADDRESS=spaceWireAddr1, SLOT=1, NINTS=1, \
            NGROUPS=groupNum1, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=128, NCOLS=128, \
            ROWCORNER={1:12.4f},COLCORNER={2:12.4f}'.format(self.guider, xarr1, yarr1))
            file_out.write(eol)

            write_from_template(self.template_c, file_out)


            file_out.write('@IFGS_CONFIG {0}, spaceWireAddr2, SLOT=2, NINTS=1, \
            NGROUPS=groupNum2, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=32, NCOLS=32, \
            ROWCORNER={1:12.4f}, COLCORNER={2:12.4f}'.format(self.guider, xarr2, yarr2))
            file_out.write(eol)

            write_from_template(self.template_d, file_out)
        file_out.close()
        print("Successfully wrote: {0}_G{1}_ACQ.prc".format(root, guider))
        shutil.copy2(os.path.join(self.out_dir, 'dhas', '{0}_G{1}_ACQ.prc'.format(root,
                                                                                  guider)),
                     os.path.join(self.out_dir, 'ground_system'))

def write_from_template(template, file_out):
    '''
    Write the lines from the specifies template to the specified file
    '''
    with open(template, 'r') as temp:
        for line in temp:
            file_out.write(line)
    temp.close()
