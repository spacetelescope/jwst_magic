#This code is adapted from mkIDproc.pro from Sherie Hofeltz.

import numpy as np
import os

#local
import grptoia

class mkproc(object):
    '''
    Makes CECIL proc files for FGS guider 1 and 2
    '''
    def __init__(self,guider,root,x,y,c,step,thresh_factor=0.5,out_dir=None,
                 template_path=None):
        '''
        Parameters
        ==========
        guider: int
            Guider "1" or "2"
        root: str
            Root for data
        reg_file: str
            <location>/<root>_psf_count_rates.txt (Include y, x, and counts for
            each star)
        step: str
            'ID' or 'ACQ'
        thresh_factor: float
            Default = 0.5
        template_path: str
            Location oftemplate files
        '''
        local_path = os.path.dirname(os.path.realpath(__file__))

        # Create output directory if does not exist
        if out_dir is None:
            self.out_dir = os.path.join(local_path,'out',root)
        else:
            self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # Find templates. If template path not given, script will assume that a
        #'templates' directory that includes are necessary prc templates lives
        #in the same directory as this script
        if template_path is None:
            template_path = os.path.join(local_path,'templates')

        self.find_templates(guider,step=step,template_path=template_path)


        # Depending on the 'step' create the correct CECIL proc files.
        if step == 'ID':
            self.create_ID_proc_file(root,guider,x,y,c,thresh_factor=thresh_factor)
        elif step == 'ACQ':
            self.create_ACQ_proc_file(root,guider)

    def find_templates(self,guider, step, template_path):
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
        self.guider='GUIDER{}'.format(guider)
        self.templateHDR= os.path.join(template_path,'g{}{}templateHDR.prc'.format(guider,step))
        self.templateA = os.path.join(template_path,'g{}{}templateA.prc'.format(guider,step))
        self.templateB = os.path.join(template_path,'g{}{}templateB.prc'.format(guider,step))
        self.templateC = os.path.join(template_path,'g{}{}templateC.prc'.format(guider,step))
        self.templateD = os.path.join(template_path,'g{}{}templateD.prc'.format(guider,step))
        if step == 'ID':
            self.templateE = os.path.join(template_path,'g{}{}templateE.prc'.format(guider,step))
            self.templateF = os.path.join(template_path,'g{}{}templateF.prc'.format(guider,step))

    def write_from_template(self,template,file_out):
        with open(template, 'r') as t:
            for line in t:
                file_out.write(line)
        t.close()

    def create_ID_proc_file(self,root,guider,x,y,c,thresh_factor,nref=9):
        eol='\n'

        sgid='_G{}'.format(guider)

        with open(os.path.join(self.out_dir,'{0}{1}_ID.prc'.format(root,sgid)), 'w') as file_out:
            self.write_from_template(self.templateHDR,file_out)

            file_out.write('PROC {0}{1}_ID'.format(root,sgid))
            file_out.write(eol)

            self.write_from_template(self.templateA,file_out)

            #y,x,c = (np.loadtxt(reg_file, delimiter=' ', skiprows=1)).T  # star coords & gs counts
            xAngle, yAngle = convert_pix_to_ideal(guider, x, y)

            threshfactor=thresh_factor
            thresh=threshfactor*c
            nstar=len(xAngle)

            #nref=nstar-1
            file_out.write('@IFGS_GUIDESTAR {0}, DFT,{1:12.4f},{2:12.4f},{3:12d},{4:8d}'.format(self.guider,xAngle[0],yAngle[0],int(c[0]),int(thresh[0])))
            file_out.write(eol)

            self.write_from_template(self.templateB,file_out)

            if nref >= 1: # ref stars > 0
                file_out.write('@IFGS_REFCOUNT DETECTOR={0}, REFSTARS={1}'.format(self.guider,nref))
                file_out.write(eol)

                self.write_from_template(self.templateC,file_out)


                for istar in range(1,nref+1):
                    file_out.write('@IFGS_REFSTAR {0},{1:5d},{2:12.6f},{3:12.6f},{4:8d},{5:8d}'.format(self.guider,int(istar-1),xAngle[istar],yAngle[istar],int(c[istar]),int(thresh[istar])))
                    file_out.write(eol)
                    self.write_from_template(self.templateD,file_out)

                istar=nstar-1
                file_out.write('@IFGS_REFSTAR {0},{1:5d},{2:12.6f},{3:12.6f},{4:8d},{5:8d}'.format(self.guider,int(istar-1),xAngle[istar],yAngle[istar],int(c[istar]),int(thresh[istar])))
                file_out.write(eol)

                self.write_from_template(self.templateE,file_out)

            self.write_from_template(self.templateF,file_out)

        file_out.close()
        print("Successfully wrote: {0}{1}_ID.prc".format(root,sgid))


    def create_ACQ_proc_file(self,root,guider):
        eol='\n'

        sgid='_G{}'.format(guider)

        i1,x1,y1,c1 = np.loadtxt(os.path.join(self.out_dir,'{}{}_ACQ1.stc'.format(root,sgid))).T #corner coords & gs counts
        threshgs=0.50*c1

        i2,x2,y2,c2 = np.loadtxt(os.path.join(self.out_dir,'{}{}_ACQ2.stc'.format(root,sgid))).T

        with open(os.path.join(self.out_dir,'{0}{1}_ACQ.prc'.format(root,sgid)), 'w') as file_out:
            self.write_from_template(self.templateHDR,file_out)

            file_out.write('PROC {0}{1}_ACQ'.format(root,sgid))
            file_out.write(eol)

            self.write_from_template(self.templateA,file_out)


            file_out.write('@IFGS_GUIDESTAR {}, 2,{:12.4f},{:12.4f},{:12.4f},{:8d}'.format(self.guider,float(x1),float(y1),float(c1),int(threshgs)))
            self.write_from_template(self.templateB,file_out)
            file_out.write('@IFGS_CONFIG {0}, SWADDRESS=spaceWireAddr1, SLOT=1, NINTS=1, \
            NGROUPS=groupNum1, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=128, NCOLS=128, \
            ROWCORNER={1:12.4f},COLCORNER={2:12.4f}'.format(self.guider,x1,y1))
            file_out.write(eol)

            self.write_from_template(self.templateC,file_out)


            file_out.write('@IFGS_CONFIG {0}, spaceWireAddr2, SLOT=2, NINTS=1, \
            NGROUPS=groupNum2, NFRAMES=1, NSAMPLES=1, GROUPGAP=1, NROWS=32, NCOLS=32, \
            ROWCORNER={1:12.4f}, COLCORNER={2:12.4f}'.format(self.guider,x2,y2))
            file_out.write(eol)

            self.write_from_template(self.templateD,file_out)
        file_out.close()
        print("Successfully wrote: {0}{1}_ACQ.prc".format(root,sgid))

def convert_pix_to_ideal(guider,x,y):
        if guider == 1:
            xAngle, yAngle = grptoia.g1RPtoIA(x,y) #this is correct: x then y
        elif guider == 2:
            xAngle, yAngle = grptoia.g2RPtoIA(x,y)
        else:
            print("How many guiders do you think FGS has? Expecting '1' or '2'")

        return xAngle, yAngle
