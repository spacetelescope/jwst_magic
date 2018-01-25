#!/usr/bin/env python

"""
SegmentGuidingTK.py

Segment Guiding tool using tkinter

Created by Colin Cox on 2017-05-01.
Modified by Lauren Chambers January 2018
"""

import glob

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pysiaf
from pysiaf.utils import rotations
from tkinter import Tk, StringVar, Radiobutton, Label, Entry, Button
import numpy as np
from astropy.io import ascii as asc
from astropy.io import fits

# Open the SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')

class SegmentForm:
    def __init__(self, segment_infile, GUI=True, GS_params_dict=None):

        # Will the tool be run through the GUI?
        self.GUI = GUI

        if self.GUI:
            # Check that you are using the right backend (or python will crash!)
            if matplotlib.get_backend() != 'TkAgg':
                errmsg = 'Cannot run GUI with current matplotlib backend (' + \
                         matplotlib.get_backend() + '). Please restart python ' \
                         'and load with TkAgg as matplotlib backend, or elect ' \
                         'to run the tool without the GUI.'
                raise EnvironmentError(errmsg)
            # Initialize the GUI
            self.initUI()

        else:
            # If not running through the GUI, define necessary attributes
            if not GS_params_dict:
                raise ValueError('If running the tool outside of the GUI, must '
                                 'supply a dictionary of the required parameters '
                                 'to the GS_params_dict argument.')

            self.V2Boff = GS_params_dict['V2Boff']
            self.V3Boff = GS_params_dict['V3Boff']
            self.fgsNum = GS_params_dict['fgsNum']
            self.RA = GS_params_dict['RA']
            self.Dec = GS_params_dict['Dec']
            self.PA = GS_params_dict['PA']
            self.segNum = GS_params_dict['segNum']

        # Get aperture parameters from FGS SIAF
        self.FGSsetup()

        # Parse the input file type
        self.parse_infile(segment_infile)

    def initUI(self):
        self.root = Tk()
        self.root.title('Segmented Guide Stars')
        self.root['bg'] = 'Snow2'

        # Choose FGS using Radio Buttons
        self.fgsNum = StringVar()
        r1 = Radiobutton(self.root, text='FGS1', variable=self.fgsNum,
                         value='1', command=self.FGSsetup)
        r1.grid(row=1)
        r2 = Radiobutton(self.root, text='FGS2', variable=self.fgsNum,
                         value='2', command=self.FGSsetup)
        r2.grid(row=1, column=1)
        r1.select()  # Default to FGS1

        # Boxes to receive FGS lookup
        Label(self.root, text='V2Ref').grid(row=2, column=1)
        Label(self.root, text='V3Ref').grid(row=2, column=2)
        Label(self.root, text='Angle').grid(row=2, column=3)
        self.fgsV2 = StringVar()
        self.fgsV3 = StringVar()
        self.fgsAngle = StringVar()
        self.fgsParity = StringVar()
        self.fgstitle = Label(self.root, text='FGS' + self.fgsNum.get())
        self.fgstitle.grid(row=3)
        self.LfgsV2 = Label(self.root, text=self.fgsV2.get())
        self.LfgsV2.grid(row=3, column=1)
        self.LfgsV3 = Label(self.root, text=self.fgsV3.get())
        self.LfgsV3.grid(row=3, column=2)
        self.LfgsA = Label(self.root, text=self.fgsAngle.get())
        self.LfgsA.grid(row=3, column=3)

        # Choose segment
        Label(self.root, text='Segment number').grid(row=4)
        self.segNum = StringVar()
        self.seg = Entry(self.root, textvariable=self.segNum)
        self.seg.grid(row=4, column=1)
        Label(self.root, text='  Segment 0 means use segment centroid').grid(row=4, column=2)
        self.segNum.trace('w', self.ChosenSeg)

        self.chSegdV2 = StringVar()
        self.chSegdV3 = StringVar()
        Label(self.root, text='Segment Offset').grid(row=5)
        self.chsegdv2 = Label(self.root, text='TBD') # Chosen segment
        self.chsegdv2.grid(row=5, column=1)
        self.chsegdv3 = Label(self.root, text='TBD')
        self.chsegdv3.grid(row=5, column=2)

        # Ideal
        Label.idl = Label(self.root, text='Ideal Coordinates').grid(row=6)
        self.segxidl = Label(self.root, text='')
        self.segxidl.grid(row=6, column=1)
        self.segyidl = Label(self.root, text='')
        self.segyidl.grid(row=6, column=2)

        # Boresight Offset
        self.V2Boff = StringVar()
        self.V3Boff = StringVar()
        Label(self.root, text='Boresight Offset').grid(row=7)
        self.EV2Boff = Entry(self.root, textvariable=self.V2Boff)
        self.EV2Boff.grid(row=7, column=1)
        self.EV3Boff = Entry(self.root, textvariable=self.V3Boff)
        self.EV3Boff.grid(row=7, column=2)
        self.V2Boff.set('0.1')
        self.V3Boff.set('0.2')
        self.V2Boff.trace('w', self.Ready)
        self.V3Boff.trace('w', self.Ready)

        # V2V3 aiming point
        Label(self.root, text='Aiming V2V3').grid(row=8)
        self.V2Aim_label = Label(self.root, text='V2')
        self.V2Aim_label.grid(row=8, column=1)
        self.V3Aim_label = Label(self.root, text='V3')
        self.V3Aim_label.grid(row=8, column=2)

        # widgets for guide star calculation
        Label(self.root, text='RA').grid(row=9, column=1)
        Label(self.root, text='Dec').grid(row=9, column=2)
        Label(self.root, text='PA').grid(row=9, column=3)
        Label(self.root, text='Guide Star').grid(row=10)
        self.RA = StringVar()
        self.Dec = StringVar()
        self.PA = StringVar()
        self.egs1 = Entry(self.root, textvariable=self.RA)
        self.egs1.grid(row=10, column=1)
        self.egs2 = Entry(self.root, textvariable=self.Dec)
        self.egs2.grid(row=10, column=2)
        self.egs3 = Entry(self.root, textvariable=self.PA)
        self.egs3.grid(row=10, column=3)

        self.RA.trace('w', self.Ready)
        self.Dec.trace('w', self.Ready)
        self.PA.trace('w', self.Ready)

        # Bottom row - Calculate, error message and Finish
        self.go = Button(self.root, text='Calculate', command=self.Calculate,
                         fg='green', state='disabled')
        self.go.grid(row=11)
        self.errmsg = Label(self.root, text='')
        self.errmsg.grid(row=11, column=1, columnspan=2)
        Button(self.root, text='Finish', command=self.Finish, fg='red').grid(row=11, column=3)

    def ChosenSeg(self, *args):
        '''
        1) Check that the user-provided segment ID number is valid (0 to 18)
        2) Calculate the central V2/V3 point of either the provided segment or
        the center of the segment array
        3) Update GUI with V2/V3 coordinates and Ideal angle coordinatess
        '''
        # Try to convert provided segment ID number to an integer
        try:
            segN = int(self.segNum.get())
        except AttributeError:
            segN = int(self.segNum)
        except ValueError:
            print('Unrecognized segment')
            if self.GUI:
                self.errmsg.configure(text='Unrecognized segment')
            return

        # If no segment ID number is provided, do nothing
        if segN == '':
            return

        # Refresh FGS data
        self.FGSsetup()

        # Ensure the provided segment ID is valid
        segMax = len(self.V2SegArray)
        if (segN < 0) or (segN > segMax):
            print('Segment number out of range')
            if self.GUI:
                self.errmsg.configure(text='Segment number out of range')
            return

        # Determine the central V2/V3 point from the given segment ID

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        if segN > 0:
            V2SegN = self.V2SegArray[segN - 1]
            V3SegN = self.V3SegArray[segN - 1]
        # Otherwise, if the input segment ID was 0, set the V2/V3 ref point to
        # be the mean of all segments' locations
        else:
            V2SegN = self.V2SegArray.mean()
            V3SegN = self.V3SegArray.mean()

        # Calculate aim
        self.V2Ref = float(self.fgsV2)   # obtained from FGS setup
        self.V3Ref = float(self.fgsV3)
        dV2Aim = V2SegN
        dV3Aim = V3SegN
        self.V2Aim = self.V2Ref + dV2Aim
        self.V3Aim = self.V3Ref + dV3Aim

        # Convert to Ideal coordinates
        xIdl, yIdl = self.fgs_siaf_aperture.tel_to_idl(self.V2Aim, self.V3Aim)

        if self.GUI:
            # Update GUI labels with V2/V3 location of chosen segment (or the
            # segment location mean, if 0 was input)
            self.chsegdv2.configure(text='%8.4f' %V2SegN)
            self.chsegdv3.configure(text='%8.4f' %V3SegN)
            self.V2Aim_label.configure(text='%8.4f' %self.V2Aim)
            self.V3Aim_label.configure(text='%8.4f' %self.V3Aim)
            self.segxidl.configure(text='%8.4f' %xIdl)
            self.segyidl.configure(text='%8.4f' %yIdl)

    def Ready(self, *args):
        if self.RA.get() == '' or self.Dec.get() == '' or self.PA.get() == '' or \
           self.V2Boff.get() == '' or self.V3Boff.get() == '':
            return

        else:
            try:
                segN = int(self.segNum.get())
            except AttributeError:
                segN = int(self.segNum)

            if segN in list(range(19)):
                if self.GUI:
                    self.go.config(state='normal')
                    self.errmsg.config(text='')
            else:
                if self.GUI:
                    self.go.config(state='disabled')
                    self.errmsg.config(text='Segment number out of range')
                print('Segment number out of range')

            return

    def FGSsetup(self, *args):
        '''Taking the current guider number (per the radio buttons on the GUI),
        extracts V2Ref, V3Ref, V3IdlYAngle, and VIdlParity from the respective
        aperture in the FGS SIAF.
        '''

        # Read chosen guider number
        if self.GUI:
            fgsN = self.fgsNum.get()
        else:
            fgsN = self.fgsNum
        det = 'FGS' + str(fgsN) + '_FULL_OSS'

        # Open SIAF aperture for appropriate guider with pysiaf
        self.fgs_siaf_aperture = FGS_SIAF[det]
        V2Ref = self.fgs_siaf_aperture.V2Ref
        V3Ref = self.fgs_siaf_aperture.V3Ref
        V3IdlYAngle = self.fgs_siaf_aperture.V3IdlYAngle
        VIdlParity = self.fgs_siaf_aperture.VIdlParity

        fgsV2 = '%10.4f' %V2Ref
        fgsV3 = '%10.4f' %V3Ref
        fgsAngle = '%10.4f' %V3IdlYAngle
        fgsParity = '%3d' %VIdlParity

        if self.GUI:
            # Update GUI and clear all calculations to force new calculation
            self.fgstitle.config(text='FGS' + fgsN)
            self.LfgsV2.config(text=fgsV2)  # Put results in Lfgs Label
            self.LfgsV3.config(text=fgsV3)
            self.LfgsA.config(text=fgsAngle)
            self.chsegdv2.configure(text='')
            self.chsegdv3.configure(text='')
            self.segxidl.configure(text='')
            self.segyidl.configure(text='')
            self.V2Aim_label.configure(text='')
            self.V3Aim_label.configure(text='')
            self.errmsg.configure(text='')  # Clear error message

        # Update class attributes
        self.fgsV2 = fgsV2
        self.fgsV3 = fgsV3
        self.fgsAngle = fgsAngle
        self.fgsParity = fgsParity

    def Calculate(self):
        # recall data read from V2V3offsets.txt
        SegIDs = self.SegIDArray
        V2Segs = self.V2SegArray
        V3Segs = self.V3SegArray
        nseg = len(SegIDs)

        # Convert V3/V3 coordinates to ideal coordinates
        xIdlSegs, yIdlSegs = self.fgs_siaf_aperture.tel_to_idl(V2Segs + self.V2Ref, V3Segs + self.V3Ref)

        # Get the guide star and boresight parameters
        V2B, V3B, gsRA, gsDec, gsPA, A = self.get_guidestar_params()

        # Get RA and Dec for each segment.
        self.SegRA = np.zeros(nseg)
        self.SegDec = np.zeros(nseg)
        for i in range(nseg):
            V2 = self.V2Ref + V2Segs[i]
            V3 = self.V3Ref + V3Segs[i]
            (self.SegRA[i], self.SegDec[i]) = rotations.pointing(A, V2, V3)
            #print('%2d %10.4f %10.4f %12.7f %12.7f' %(i+1, V2, V3, self.SegRA[i], self.SegDec[i]))

        if self.GUI:
            self.errmsg.configure(text='Calculation complete')

        # Convert segment coordinates to detector frame
        xDet, yDet = self.fgs_siaf_aperture.idl_to_det(xIdlSegs, yIdlSegs)

        # Check to make sure no segments are off the detector
        for x, y, i_seg in zip(xDet, yDet, SegIDs):
            if x < 0.5 or x > 2048.5:
                print('WARNING: %8s off detector in X direction' %i_seg)
            if y < 0.5 or y > 2048.5:
                print('WARNING: %8s off detector in Y direction' %i_seg)

        # Retrieve FGS parameters obtained in FGSsetup
        if self.GUI:
            fgsN = self.fgsNum.get()
        else:
            fgsN = self.fgsNum

        # Summary output
        print('\nSummary')
        print('Aperture FGS', fgsN)
        print('V2Ref %s V3Ref %s arc-sec IdlAngle %s degrees' %(self.fgsV2, self.fgsV3, self.fgsAngle))
        print('Used segment', self.segNum)
        print('Boresight offset', V2B, V3B, 'arc-sec')
        print('Guide star at RA %s  Dec %s degrees' %(gsRA, gsDec))
        print('Position angle %s degrees' %gsPA)
        print('\nSegment     dV2    dV3    xIdl   yIdl     RA         Dec         xDet     yDet')

        for p in range(nseg):
            print('%5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f' \
            %(SegIDs[p], V2Segs[p], V3Segs[p], xIdlSegs[p], yIdlSegs[p],
              self.SegRA[p], self.SegDec[p], xDet[p], yDet[p]))

        # Final output
        print('\nFinal Output')
        with open('SegmentGuiding.txt', 'w') as sg:
            rate = 0.0  # placeholder for count rate
            for p in range(nseg):
                part1 = '-star%02d = %12.6f %12.6f %8.2f  ' %(p + 1, self.SegRA[p],
                                                              self.SegDec[p], rate)
                print(part1, end='')
                sg.write(part1)
                onDet = []  # List of other segments on detector
                for q in range(nseg):
                    if (q != p) and (1 <= xDet[q] <= 2048) and (1 <= yDet[q] <= 2048):
                        onDet.append(q + 1)
                ns = len(onDet)
                for q in range(ns - 1):
                    part2 = '%2d, ' % onDet[q]
                    print(part2, end=' ')
                    sg.write(part2)
                part3 = '%2d' %onDet[ns - 1]
                print(part3)
                sg.write(part3 + '\n')

        self.plot_segments()

    def parse_infile(self, segment_infile):

        print('Segment data read from {}'.format(segment_infile))

        # If the input file is a .txt file, parse the file
        if segment_infile[-4:] == '.txt':
            read_table = asc.read(segment_infile)
            column_names = read_table.colnames

            if (any(['V2' in c for c in column_names])) and \
               (any(['V3' in c for c in column_names])):
                print('Input file = V2/V3 coordinates')
                segment_coords = read_table

            elif (any(['x' in c for c in column_names])) and \
                 (any(['y' in c for c in column_names])):
                print('Input file = pixel coordinates')
                # I need to convert from pixels to V2/V3 here??

                v2v3_table = np.copy(read_table)

                segment_coords = v2v3_table

        # if it is a FITS file...
        elif segment_infile[-5:] == '.fits':
            print('Input file = FITS image')

            fits_data = fits.getdata(segment_infile)
            fits_data[fits_data < 0] = 0.1

            # Run the tool!
            root = None
            guider = 1
            GA = True
            select_psfs.create_reg_file(fits_data, root, guider,
                                        '/Users/lchambers/TEL/FGS/Commissioning-tools/jwst_fgs_commissioning_tools/segment_guiding/',
                                        global_alignment=GA)
            all_psfs_file = '{0}_G{1}_ALLpsfs.txt'.format(segment_infile[:-5], guider)

            read_table = asc.read(all_psfs_file)

            # now translate the x/y cols into V2/V3
            v2v3_table = np.copy(read_table)

            segment_coords = v2v3_table

        # Determine the IDs and coordinates of all segments
        self.SegIDArray = segment_coords['SegID']
        self.V2SegArray = segment_coords['V2Seg']
        self.V3SegArray = segment_coords['V3Seg']
        print('{} Segments in input file.'.format(len(self.SegIDArray)))

    def plot_segments(self):
        # Plot segments in V2/V3 frame
        plt.figure(1)
        plt.clf()
        plt.plot(self.V2SegArray, self.V3SegArray, 'b*')
        plt.plot(self.V2SegArray.mean(), self.V2SegArray.mean(), 'ro')
        plt.grid(True)
        plt.axis([40.0, -40.0, -40.0, 40.0])  # V2 to the left
        plt.title('Segments')
        plt.xlabel('<-- Delta V2')
        plt.ylabel('Delta V3 -->')
        for s in range(len(self.V2SegArray)):
            plt.text(self.V2SegArray[s], self.V3SegArray[s], self.SegIDArray[s])
        plt.savefig('Segments.png')

        # Plot calculate segments' RA and Dec
        plt.figure(2)
        plt.clf()
        plt.plot(self.SegRA, self.SegDec, 'b*')
        plt.grid(True)
        plt.title('Segment RA and Dec')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        RAmean = self.SegRA.mean()
        Decmean = self.SegDec.mean()
        plt.plot(RAmean, Decmean, 'ro')
        if self.GUI:
            segN = int(self.segNum.get())
        else:
            segN = int(self.segNum)
        if segN > 0:
            plt.plot(self.SegRA[segN - 1], self.SegDec[segN - 1], 'mx', markersize=12)
        for i in range(len(self.V2SegArray)):
            plt.text(self.SegRA[i], self.SegDec[i], str(i + 1))
        plt.gca().invert_xaxis()
        plt.gca().ticklabel_format(useOffset=False)
        plt.savefig('SegmentSky.png')

    def get_guidestar_params(self):

        if self.GUI:
            V2B = self.V2Boff.get()
            V3B = self.V3Boff.get()
            gsRA = self.RA.get()
            gsDec = self.Dec.get()
            gsPA = self.PA.get()
        else:
            V2B = self.V2Boff
            V3B = self.V3Boff
            gsRA = self.RA
            gsDec = self.Dec
            gsPA = self.PA

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BORESIGHT OFFSET ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        # Guide star information
        msg = ["OK", "Boresight parameter conversion error",
               "Boresight parameter out of range"]

        errcode = checkout(V2B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.errmsg.configure(text=error)
            return error
        else:
            V2B = float(V2B)

        errcode = checkout(V3B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.errmsg.configure(text=error)
            return error
        else:
            V3B = float(V3B)

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ RA, DEC, AND POSITION ANGLE ~ ~ ~ ~ ~ ~

        msg = ["OK", "Guide Star parameter conversion error",
               "Guide Star parameter out of range"]

        errcode = checkout(gsRA, 0.0, 360.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.errmsg.configure(text=error)
            return error
        else:
            gsRA = float(gsRA)

        errcode = checkout(gsDec, -90.0, 90.0)
        error = msg[errcode]
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.errmsg.configure(text=error)
            return error
        else:
            gsDec = float(gsDec)

        errcode = checkout(gsPA, -180.0, 180.0)
        error = msg[errcode]
        if errcode != 0:
            print(error)
            if self.GUI:
                self.errmsg.configure(text=error)
            return
        else:
            if self.GUI:
                self.errmsg.configure(text='')  # Clear error message
            # Determine the attitude matrix
            A = rotations.attitude(self.V2Aim + V2B, self.V3Aim + V3B, gsRA,
                                   gsDec, float(gsPA))

        return V2B, V3B, gsRA, gsDec, gsPA, A

    def Show(self):
        self.root.mainloop()

    def Finish(self):
        print('Stopping')
        self.root.quit()  # frees up iPython window
        self.root.destroy()  # closes GUI

############################## End Class SegmentForm ###############################

def checkout(str, low, high):
    """Test conversion from string to float.
    if float conversion works, test range
    return errcode 0 for OK, 1 for conversion error, 2 for range error"""

    try:
        x = float(str)
        if low <= x <= high:
            errcode = 0
        else:
            errcode = 2
    except ValueError:
        errcode = 1
    return errcode

def run_tool(segment_infile, GUI=True, GS_params_dict=None):
    if not GS_params_dict and not GUI:
        GS_params_dict = {'V2Boff': 0.1,  # V2 boresight offset
                          'V3Boff': 0.2,  # V3 boresight offset
                          'fgsNum': 1,  # guider number
                          'RA': 30.,  # RA of guide star
                          'Dec': 50.,  # Dec of guide star
                          'PA': 2.,  # position angle of guide star
                          'segNum': 0}  # selected segment to guide on

    # Set up segment form class
    sgForm = SegmentForm(segment_infile, GUI=GUI, GS_params_dict=GS_params_dict)

    # Either run the GUI or run the calculation
    if GUI:
        sgForm.Show()
    else:
        sgForm.ChosenSeg()
        sgForm.Calculate()

############################## Main Program - Actions #####################################


if __name__ == '__main__':
    segment_infile = 'V2V3offsets.txt'
    run_tool(segment_infile)
