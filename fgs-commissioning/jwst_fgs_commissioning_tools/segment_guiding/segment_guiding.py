#!/usr/bin/env python

"""
segment_guiding.py

Segment Guiding Tool (SGT), optionally using tkinter GUI

Created by Colin Cox on 2017-05-01.
Modified by Lauren Chambers January 2018
"""

import os

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pysiaf
from pysiaf.utils import rotations
from tkinter import Tk, StringVar, Radiobutton, Label, Entry, Button
import numpy as np
from astropy.io import ascii as asc
from functools import partial

from jwst_fgs_commissioning_tools import coordinate_transforms

# Establish segment guiding files directory
LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PATH = os.path.split(LOCAL_PATH)[0]
SGT_FILES_PATH = os.path.join(os.path.split(PACKAGE_PATH)[0], 'segment_guiding_files')

# Open the SIAF with pysiaf
FGS_SIAF = pysiaf.Siaf('FGS')

class SegmentGuidingCalculator:
    def __init__(self, segment_infile, root=None, GUI=True, GS_params_dict=None):

        # If no root provided, name the output the same as the input
        if not root:
            self.root = os.path.basename(segment_infile).split('.')[0]
        else:
            self.root = root

        # Will the tool be run through the GUI?
        self.GUI = GUI
        if self.GUI:
            # If so, initialize the GUI object
            self.SegmentGuidingGUI = SegmentGuidingGUI(self)
        else:
            # If not running through the GUI, define necessary attributes
            if not GS_params_dict:
                raise ValueError('If running the tool outside of the GUI, must '
                                 'supply a dictionary of the required parameters '
                                 'to the GS_params_dict argument.')

            # Map GS_params_dict keys to attributes
            for attr_name in GS_params_dict.keys():
                setattr(self, attr_name, GS_params_dict[attr_name])

        # Get aperture parameters from FGS SIAF
        self.FGSsetup()

        # Parse the input file type
        self.parse_infile(segment_infile)

    def ChosenSeg(self, *args):
        '''
        1) Check that the user-provided segment ID number is valid (0 to 18)
        2) Calculate the central V2/V3 point of either the provided segment or
        the center of the segment array
        3) Update GUI with V2/V3 coordinates and Ideal angle coordinatess
        '''
        # Try to convert provided segment ID number to an integer
        try:
            segN = int(self.segNum)
        except ValueError:
            print('Unrecognized segment')
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text='Unrecognized segment')
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
                self.SegmentGuidingGUI.errmsg.configure(text='Segment number out of range')
            return

        # Determine the central V2/V3 point from the given segment ID

        # If a specific segment was provided, set the V2/V3 ref point to be
        # that segment's location
        if segN > 0:
            self.V2SegN = self.V2SegArray[segN - 1]
            self.V3SegN = self.V3SegArray[segN - 1]
        # Otherwise, if the input segment ID was 0, set the V2/V3 ref point to
        # be the mean of all segments' locations
        else:
            self.V2SegN = self.V2SegArray.mean()
            self.V3SegN = self.V3SegArray.mean()

        # Calculate aim
        self.V2Ref = float(self.fgsV2)   # obtained from FGS setup
        self.V3Ref = float(self.fgsV3)
        dV2Aim = self.V2SegN
        dV3Aim = self.V3SegN
        self.V2Aim = self.V2Ref + dV2Aim
        self.V3Aim = self.V3Ref + dV3Aim

        # Convert to Ideal coordinates
        self.xIdl, self.yIdl = self.fgs_siaf_aperture.tel_to_idl(self.V2Aim, self.V3Aim)

        if self.GUI:
            # Update GUI labels with V2/V3 location of chosen segment (or the
            # segment location mean, if 0 was input)
            self.SegmentGuidingGUI.redefine_vars('ChosenSeg')

    def FGSsetup(self, *args):
        '''Taking the current guider number (per the radio buttons on the GUI),
        extracts V2Ref, V3Ref, V3IdlYAngle, and VIdlParity from the respective
        aperture in the FGS SIAF.
        '''

        # Read chosen guider number
        # Ensure the guider number is valid
        if str(self.fgsNum) not in ['1', '2']:
            raise ValueError('Invalid guider number: "{}"'.format(self.fgsNum))

        det = 'FGS' + str(self.fgsNum) + '_FULL_OSS'

        # Open SIAF aperture for appropriate guider with pysiaf
        self.fgs_siaf_aperture = FGS_SIAF[det]
        V2Ref = self.fgs_siaf_aperture.V2Ref
        V3Ref = self.fgs_siaf_aperture.V3Ref
        V3IdlYAngle = self.fgs_siaf_aperture.V3IdlYAngle
        VIdlParity = self.fgs_siaf_aperture.VIdlParity

        self.fgsV2 = '%10.4f' % V2Ref
        self.fgsV3 = '%10.4f' % V3Ref
        self.fgsAngle = '%10.4f' % V3IdlYAngle
        self.fgsParity = '%3d' % VIdlParity

        if self.GUI:
            self.SegmentGuidingGUI.redefine_vars('FGSsetup')

    def Calculate(self):
        # recall data read from infile
        SegIDs = self.SegIDArray
        V2Segs = self.V2SegArray
        V3Segs = self.V3SegArray
        nseg = len(SegIDs)

        # Convert V3/V3 coordinates to ideal coordinates
        xIdlSegs, yIdlSegs = self.fgs_siaf_aperture.tel_to_idl(V2Segs + self.V2Ref,
                                                               V3Segs + self.V3Ref)

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
            self.SegmentGuidingGUI.errmsg.configure(text='Calculation complete')

        # Convert segment coordinates to detector frame
        xDet, yDet = self.fgs_siaf_aperture.idl_to_det(xIdlSegs, yIdlSegs)

        # Check to make sure no segments are off the detector
        for x, y, i_seg in zip(xDet, yDet, SegIDs):
            if x < 0.5 or x > 2048.5:
                print('WARNING: %8s off detector in X direction' % i_seg)
            if y < 0.5 or y > 2048.5:
                print('WARNING: %8s off detector in Y direction' % i_seg)

        # Summary output
        print('\nSummary')
        print('Aperture FGS', self.fgsNum)
        print('V2Ref %s V3Ref %s arc-sec IdlAngle %s degrees' % (self.fgsV2,
                                                                 self.fgsV3,
                                                                 self.fgsAngle))
        print('Used segment', self.segNum)
        print('Boresight offset', V2B, V3B, 'arc-sec')
        print('Guide star at RA %s  Dec %s degrees' % (gsRA, gsDec))
        print('Position angle %s degrees' % gsPA)
        print('\nSegment     dV2    dV3    xIdl   yIdl     RA         Dec         xDet     yDet')

        for p in range(nseg):
            print('%5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f' \
                  % (SegIDs[p], V2Segs[p], V3Segs[p], xIdlSegs[p], yIdlSegs[p],
                     self.SegRA[p], self.SegDec[p], xDet[p], yDet[p]))

        # Print and save final output
        out_file = os.path.join(SGT_FILES_PATH, self.root + '_segmentguiding_visit.txt')
        with open(out_file, 'w') as sg:
            print('\nFinal Output:')
            rate = 0.0  # placeholder for count rate
            for p in range(nseg):
                part1 = '-star%02d = %12.6f %12.6f %8.2f  ' % (p + 1, self.SegRA[p],
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
                part3 = '%2d' % onDet[ns - 1]
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

        else:
            raise TypeError('Incompatible file type: ', segment_infile)

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
        plt.savefig(os.path.join(SGT_FILES_PATH, self.root + '_V2V3segments.png'))

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

        segN = int(self.segNum)
        if segN > 0:
            plt.plot(self.SegRA[segN - 1], self.SegDec[segN - 1], 'mx', markersize=12)
        for i in range(len(self.V2SegArray)):
            plt.text(self.SegRA[i], self.SegDec[i], str(i + 1))
        plt.gca().invert_xaxis()
        plt.gca().ticklabel_format(useOffset=False)
        plt.savefig(os.path.join(SGT_FILES_PATH, self.root + '_RADecsegments.png'))

    def get_guidestar_params(self):
        V2B = self.V2Boff
        V3B = self.V3Boff
        gsRA = self.RA
        gsDec = self.Dec
        gsPA = self.PA

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BORESIGHT OFFSET ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        # Guide star information
        msg = ["OK", "Boresight parameter conversion error",
               "Boresight parameter out of range"]

        errcode = self.checkout(V2B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text=error)
            return error
        else:
            V2B = float(V2B)

        errcode = self.checkout(V3B, -10.0, 10.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text=error)
            return error
        else:
            V3B = float(V3B)

        # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ RA, DEC, AND POSITION ANGLE ~ ~ ~ ~ ~ ~

        msg = ["OK", "Guide Star parameter conversion error",
               "Guide Star parameter out of range"]

        errcode = self.checkout(gsRA, 0.0, 360.0)
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text=error)
            return error
        else:
            gsRA = float(gsRA)

        errcode = self.checkout(gsDec, -90.0, 90.0)
        error = msg[errcode]
        if errcode != 0:
            error = msg[errcode]
            print(error)
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text=error)
            return error
        else:
            gsDec = float(gsDec)

        errcode = self.checkout(gsPA, -180.0, 180.0)
        error = msg[errcode]
        if errcode != 0:
            print(error)
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text=error)
            return
        else:
            if self.GUI:
                self.SegmentGuidingGUI.errmsg.configure(text='')  # Clear error message
            # Determine the attitude matrix
            A = rotations.attitude(self.V2Aim + V2B, self.V3Aim + V3B, gsRA,
                                   gsDec, float(gsPA))

        return V2B, V3B, gsRA, gsDec, gsPA, A

    def checkout(self, str, low, high):
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

class SegmentGuidingGUI(SegmentGuidingCalculator):
    def __init__(self, calculator):
        # Check that you are using the right backend (or python will crash!)
        if matplotlib.get_backend() != 'TkAgg':
            errmsg = 'Cannot run GUI with current matplotlib backend (' + \
                     matplotlib.get_backend() + '). Please restart python ' \
                     'and load with TkAgg as matplotlib backend, or elect ' \
                     'to run the tool without the GUI.'
            raise EnvironmentError(errmsg)

        # Inherit attributes from calculator object
        self.calculator = calculator

        # Initialize the GUI
        self.initUI()

    def initUI(self):
        self.root_window = Tk()
        self.root_window.title('Segmented Guide Stars')
        self.root_window['bg'] = 'Snow2'

        # Choose FGS using Radio Buttons
        self.fgsNum = StringVar()
        self.fgsNum.trace('w', partial(self.callback, 'fgsNum'))
        r1 = Radiobutton(self.root_window, text='FGS1', variable=self.fgsNum,
                         value='1', command=self.calculator.FGSsetup)
        r1.grid(row=1)
        r2 = Radiobutton(self.root_window, text='FGS2', variable=self.fgsNum,
                         value='2', command=self.calculator.FGSsetup)
        r2.grid(row=1, column=1)
        r1.select()  # Default to FGS1

        # Boxes to receive FGS lookup
        Label(self.root_window, text='V2Ref').grid(row=2, column=1)
        Label(self.root_window, text='V3Ref').grid(row=2, column=2)
        Label(self.root_window, text='Angle').grid(row=2, column=3)
        self.fgsV2 = StringVar()
        self.fgsV3 = StringVar()
        self.fgsAngle = StringVar()
        self.fgsParity = StringVar()
        self.fgstitle = Label(self.root_window, text='FGS' + self.fgsNum.get())
        self.fgstitle.grid(row=3)
        self.LfgsV2 = Label(self.root_window, text=self.fgsV2.get())
        self.LfgsV2.grid(row=3, column=1)
        self.LfgsV3 = Label(self.root_window, text=self.fgsV3.get())
        self.LfgsV3.grid(row=3, column=2)
        self.LfgsA = Label(self.root_window, text=self.fgsAngle.get())
        self.LfgsA.grid(row=3, column=3)

        # Choose segment
        Label(self.root_window, text='Segment number').grid(row=4)
        self.segNum = StringVar()
        self.seg = Entry(self.root_window, textvariable=self.segNum)
        self.seg.grid(row=4, column=1)
        Label(self.root_window, text='  Segment 0 means use segment centroid').grid(row=4, column=2)
        self.segNum.trace('w', self.calculator.ChosenSeg)
        self.segNum.trace('w', partial(self.callback, 'segNum'))

        self.chSegdV2 = StringVar()
        self.chSegdV3 = StringVar()
        Label(self.root_window, text='Segment Offset').grid(row=5)
        self.chsegdv2 = Label(self.root_window, text='TBD')  # Chosen segment
        self.chsegdv2.grid(row=5, column=1)
        self.chsegdv3 = Label(self.root_window, text='TBD')
        self.chsegdv3.grid(row=5, column=2)

        # Ideal
        Label.idl = Label(self.root_window, text='Ideal Coordinates').grid(row=6)
        self.segxidl = Label(self.root_window, text='')
        self.segxidl.grid(row=6, column=1)
        self.segyidl = Label(self.root_window, text='')
        self.segyidl.grid(row=6, column=2)

        # Boresight Offset
        self.V2Boff = StringVar()
        self.V3Boff = StringVar()
        Label(self.root_window, text='Boresight Offset').grid(row=7)
        self.EV2Boff = Entry(self.root_window, textvariable=self.V2Boff)
        self.EV2Boff.grid(row=7, column=1)
        self.EV3Boff = Entry(self.root_window, textvariable=self.V3Boff)
        self.EV3Boff.grid(row=7, column=2)
        self.V2Boff.set('0.1')
        self.V3Boff.set('0.2')
        self.V2Boff.trace('w', self.Ready)
        self.V2Boff.trace('w', partial(self.callback, 'V2Boff'))
        self.V3Boff.trace('w', self.Ready)
        self.V3Boff.trace('w', partial(self.callback, 'V3Boff'))

        # V2V3 aiming point
        Label(self.root_window, text='Aiming V2V3').grid(row=8)
        self.V2Aim_label = Label(self.root_window, text='V2')
        self.V2Aim_label.grid(row=8, column=1)
        self.V3Aim_label = Label(self.root_window, text='V3')
        self.V3Aim_label.grid(row=8, column=2)

        # widgets for guide star calculation
        Label(self.root_window, text='RA').grid(row=9, column=1)
        Label(self.root_window, text='Dec').grid(row=9, column=2)
        Label(self.root_window, text='PA').grid(row=9, column=3)
        Label(self.root_window, text='Guide Star').grid(row=10)
        self.RA = StringVar()
        self.Dec = StringVar()
        self.PA = StringVar()
        self.egs1 = Entry(self.root_window, textvariable=self.RA)
        self.egs1.grid(row=10, column=1)
        self.egs2 = Entry(self.root_window, textvariable=self.Dec)
        self.egs2.grid(row=10, column=2)
        self.egs3 = Entry(self.root_window, textvariable=self.PA)
        self.egs3.grid(row=10, column=3)

        self.RA.trace('w', self.Ready)
        self.RA.trace('w', partial(self.callback, 'RA'))
        self.Dec.trace('w', self.Ready)
        self.Dec.trace('w', partial(self.callback, 'Dec'))
        self.PA.trace('w', self.Ready)
        self.PA.trace('w', partial(self.callback, 'PA'))

        # Bottom row - Calculate, error message and Finish
        self.go = Button(self.root_window, text='Calculate', command=self.calculator.Calculate,
                         fg='green', state='disabled')
        self.go.grid(row=11)
        self.errmsg = Label(self.root_window, text='')
        self.errmsg.grid(row=11, column=1, columnspan=2)
        Button(self.root_window, text='Finish', command=self.Finish, fg='red').\
            grid(row=11, column=3)

    def callback(self, var_name, *args):
        if var_name == 'fgsNum':
            self.calculator.fgsNum = self.fgsNum.get()

        elif var_name == 'V2Boff':
            self.calculator.V2Boff = self.V2Boff.get()

        elif var_name == 'V3Boff':
            self.calculator.V3Boff = self.V3Boff.get()

        elif var_name == 'RA':
            self.calculator.RA = self.RA.get()

        elif var_name == 'Dec':
            self.calculator.Dec = self.Dec.get()

        elif var_name == 'PA':
            self.calculator.PA = self.PA.get()

        elif var_name == 'segNum':
            self.calculator.segNum = self.segNum.get()

    def redefine_vars(self, function):
        # Update GUI and clear all calculations to force new calculation

        # In FGSsetup
        if function == 'FGSsetup':
            self.fgstitle.config(text='FGS' + self.calculator.fgsNum)
            self.LfgsV2.config(text=self.calculator.fgsV2)  # Put results in Lfgs Label
            self.LfgsV3.config(text=self.calculator.fgsV3)
            self.LfgsA.config(text=self.calculator.fgsAngle)

        # In chosenseg
        if function == 'ChosenSeg':
            self.chsegdv2.configure(text='%8.4f' % self.calculator.V2SegN)
            self.chsegdv3.configure(text='%8.4f' % self.calculator.V3SegN)
            self.V2Aim_label.configure(text='%8.4f' % self.calculator.V2Aim)
            self.V3Aim_label.configure(text='%8.4f' % self.calculator.V3Aim)
            self.segxidl.configure(text='%8.4f' % self.calculator.xIdl)
            self.segyidl.configure(text='%8.4f' % self.calculator.yIdl)

    def Ready(self, *args):
        if self.RA.get() == '' or self.Dec.get() == '' or self.PA.get() == '' or \
           self.V2Boff.get() == '' or self.V3Boff.get() == '':
            return

        else:
            segN = int(self.segNum.get())
            if segN in list(range(19)):
                self.go.config(state='normal')
                self.errmsg.config(text='')

    def Show(self):
        self.root_window.mainloop()

    def Finish(self):
        print('Stopping')
        self.root_window.quit()  # frees up iPython window
        self.root_window.destroy()  # closes GUI

############################## End Class SegmentForm ###############################


def run_tool(segment_infile, root=None, GUI=True, GS_params_dict=None):

    if not GS_params_dict and not GUI:
        GS_params_dict = {'V2Boff': 0.1,  # V2 boresight offset
                          'V3Boff': 0.2,  # V3 boresight offset
                          'fgsNum': 1,  # guider number
                          'RA': 30.,  # RA of guide star
                          'Dec': 50.,  # Dec of guide star
                          'PA': 2.,  # position angle of guide star
                          'segNum': 0}  # selected segment to guide on

    # Set up guiding calculator object
    sg = SegmentGuidingCalculator(segment_infile, root=root, GUI=GUI,
                                  GS_params_dict=GS_params_dict)

    # Either run the GUI or run the calculation
    if GUI:
        sg.SegmentGuidingGUI.Show()
    else:
        sg.ChosenSeg()
        sg.Calculate()


############################## Main Program - Actions #####################################


if __name__ == '__main__':
    segment_infile = os.path.join(SGT_FILES_PATH, 'SGTintegrationexample.txt')
    run_tool(segment_infile)
