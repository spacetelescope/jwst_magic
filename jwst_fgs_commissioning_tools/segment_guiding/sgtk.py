#!/usr/bin/env python

"""
sgtk.py

Segment Guiding tool using tkinter

Created by Colin Cox on 2017-05-01.
"""
from math import *
from tkinter import *
#from tkinter.ttk import *

import numpy as np
from astropy.io import ascii
import glob
import matplotlib.pyplot as pl
import xml.etree.cElementTree as etree

# Private packages 
import fulltransform.rotations
import fulltransform.polynomial

class SegmentForm :
    def __init__(self) :
        self.root = None
        
    def CreateForm(self) :
        self.root = Tk()
        self.root.title('Segmented guide stars')
        
        # Choose FGS using Radio Buttons
        self.fgsNum = StringVar()
        r1 = Radiobutton(self.root, text='FGS1', variable=self.fgsNum, value='1', command=FGSsetup)
        r1.grid(row=1)
        r2 = Radiobutton(self.root, text='FGS2', variable=self.fgsNum, value='2', command=FGSsetup)
        r2.grid(row=1,column=1)
        r1.select()  # Default to FGS1
        
        # Boxes to receive FGS lookup
        Label(self.root,text='V2Ref').grid(row=2,column=1)
        Label(self.root,text='V3Ref').grid(row=2,column=2)
        Label(self.root,text='Angle').grid(row=2,column=3)
        self.fgsV2=StringVar()
        self.fgsV3=StringVar()
        self.fgsAngle=StringVar()
        self.fgsParity=StringVar()
        self.fgstitle=Label(self.root, text='FGS'+ self.fgsNum.get())
        self.fgstitle.grid(row=3)
        self.LfgsV2=Label(self.root, text=self.fgsV2.get())
        self.LfgsV2.grid(row=3,column=1)
        self.LfgsV3 = Label(self.root, text=self.fgsV3.get())
        self.LfgsV3.grid(row=3,column=2)
        self.LfgsA =Label(self.root,text=self.fgsAngle.get())
        self.LfgsA.grid(row=3,column=3)

        fgsPars=FGSsetup()
        self.fgsV2 = fgsPars[0]
        self.fgsV3 = fgsPars[1]
        self.fgsAngle = fgsPars[2]
        self.fgsParity = fgsPars[3]

        # Choose segment
        self.segNum = StringVar()
        Label(self.root, text = 'Segment number').grid(row=4)
        Entry(self.root,textvariable=self.segNum).grid(row=4,column=1)
        Label(self.root, text='  Segment 0 means use segment centroid').grid(row=4,column=2)
        self.segNum.trace('w', self.Ready)
        
        # widgets for guide star calculation
        Label(self.root, text='RA').grid(row=5, column=1)
        Label(self.root, text='Dec').grid(row=5, column=2)
        Label(self.root, text='PA').grid(row=5, column=3)
        Label(self.root, text='Guide Star').grid(row=6)
        self.RA = StringVar()
        self.Dec = StringVar()
        self.pa = StringVar()
        self.egs1 = Entry(self.root, textvariable=self.RA)
        self.egs1.grid(row=6, column=1)
        self.egs2 = Entry(self.root, textvariable=self.Dec)
        self.egs2.grid(row=6, column=2)
        self.egs3 = Entry(self.root, textvariable=self.pa)
        self.egs3.grid(row=6, column=3)

        self.RA.trace('w', self.Ready)
        self.Dec.trace('w', self.Ready)
        self.pa.trace('w', self.Ready)
        
        # Bottom row - Calculate, error message and Finish        
        self.go = Button(self.root, text='Calculate', command=Calculate, fg='green', state='disabled')
        self.go.grid(row=7, column=0)
        self.errmsg = Label(self.root, text='')
        self.errmsg.grid(row=7, column=1, columnspan=2)
        Button(self.root, text='Finish', command=self.Finish, fg='red').grid(row=7,column=3)

        (self.SegIDArray, self.V2SegArray,self.V3SegArray) = showsegments()

        
    def Ready(self, *args):
        if self.segNum.get() == '' \
        or self.RA.get() == '' \
        or self.Dec.get() == '' \
        or self.pa.get() == '':  return
        
        else:
            self.seg =int(self.segNum.get())
            if self.seg in list(range(19)) :
                self.go.config(state='normal')
                self.errmsg.config(text = '')
                
            else:
                self.go.config(state='disabled')
                self.errmsg.config(text = 'Segment number out of range')
                #print('check values')
                return        

    def Show(self):
        self.root.mainloop()


    def Finish(self):  
        print ('Stopping')      
        self.root.destroy()
        

############################## End Class SegmentForm ###############################

def FGSsetup():
    global siafData, r
    # Chosen FGS and segment 
    fgsN = sgForm.fgsNum.get()
    print ('FGS', fgsN)
    det = 'FGS' + fgsN +'_FULL_OSS'
    print (det)

    xmlList = glob.glob('/itar/jwst/tel/share/SIAF_WG/Instruments/FGS/FGS_SIAF_20*.xml')
    if len(xmlList) < 1:
        print ('No xml file found')
        return
        
    xmlList.sort() # Will normally be in date order
    xmlFile = xmlList[-1]  
    print ('Using xml file ', xmlFile) 
    siaf = etree.parse(xmlFile)
    root = siaf.getroot()
    siafData = {}               # Dictionary to store all SIAF data

    # Read  FGS xml file for item in root[0]:
    for item in root[0]:
        columnName = item.tag
        columnList = []
        for e in root.iter(columnName): columnList.append(e.text)
        columnArray = np.array(columnList) # Convert to array for easier processing
        siafData[columnName] = columnArray    # Each dictionary column now an array indexed by column name
        
    row = np.where(siafData['AperName'] == det)
    r =  (row[0][0])
    #print ('Row', r)
    V2Ref = float(siafData['V2Ref'][r])
    V3Ref = float(siafData['V3Ref'][r])
    V3IdlYAngle = float(siafData['V3IdlYAngle'][r])
    VIdlParity = int(siafData['VIdlParity'][r])
    #print ('Reference point %10.4f %10.4f arcsec' %(V2Ref,V3Ref))
    #print ('Angle %10.4f degrees Parity %2d' %(V3IdlYAngle, VIdlParity))
    fgsV2 = '%10.4f' %V2Ref
    fgsV3 = '%10.4f' %V3Ref
    fgsAngle = '%10.4f' %V3IdlYAngle
    fgsParity = '%3d' %VIdlParity
    print (fgsV2,fgsV3,fgsAngle,fgsParity)
    sgForm.fgstitle.config(text='FGS'+ fgsN)
    sgForm.LfgsV2.config(text=fgsV2) # Put results in Lfgs Label
    sgForm.LfgsV3.config(text=fgsV3) 
    sgForm.LfgsA.config(text=fgsAngle) 
    return (fgsV2, fgsV3,fgsAngle, fgsParity)

def Calculate():
    global siafData, r
    # recall data read from V2V3offsets.txt
    SegIDs = sgForm.SegIDArray
    V2Segs = sgForm.V2SegArray
    V3Segs = sgForm.V3SegArray
    V2mean = V2Segs.mean()
    V3mean = V3Segs.mean()
    theta = radians(float(sgForm.fgsAngle))
    VIdlParity = int(sgForm.fgsParity)
    xIdlSegs = VIdlParity*(V2Segs*cos(theta) - V3Segs*sin(theta))
    yIdlSegs = V2Segs*sin(theta) + V3Segs*cos(theta)
    print('\n Segment Positions')
    print (' Segment     dV2        dV3        XIdl       YIdl')
    nseg = len(SegIDs)
    for p in range(nseg): print ('%8s %10.4f %10.4f %10.4f %10.4f' %(SegIDs[p],V2Segs[p], V3Segs[p], xIdlSegs[p], yIdlSegs[p]))
    print ('Mean position %8.3f %8.3f\n' %(V2mean,V3mean))
    

    # Guide star information
    msg = ["OK", "GS parameter conversion error", "GS parameter out of range"]
    gsRA = sgForm.RA.get() # This will be a text string
    errcode = checkout(gsRA, 0.0, 360.0)
    if errcode != 0:
        error = msg[errcode]
        print (error)
        sgForm.errmsg.configure(text=error)
        return error
    gsRA = float(gsRA)
               
    gsDec = sgForm.Dec.get() # This will be a text string
    errcode = checkout(gsDec, -90.0, 90.0)
    error = msg[errcode]
    if errcode != 0:
        print (error)
        sgForm.errmsg.configure(text=error)
        return error
    gsDec = float(gsDec)

    gspa = sgForm.pa.get() # This will be a text string
    errcode = checkout(gspa, -180.0, 180.0)
    error = msg[errcode]
    if errcode != 0:
        print (error)
        sgForm.errmsg.configure(text=error)
        return
    gspa = float(gspa)
    print ('Guide star', gsRA, gsDec, gspa)

    sgForm.errmsg.configure(text='') # Clear error message
    segN = int(sgForm.segNum.get())
    
    print ('segN', segN)    
    
    if segN > 0:
        dV2Aim = V2Segs[segN-1]
        dV3Aim = V3Segs[segN-1]
    else: # if SegID is now set to -1, use mean values
        dV2Aim = V2mean
        dV3Aim = V3mean
    print ('dV', dV2Aim, dV3Aim)

    V2Ref = float(sgForm.fgsV2)   # obtained from FGS setup
    V3Ref = float(sgForm.fgsV3)
    print ('Ref',V2Ref, V3Ref)
    
    V2Aim = V2Ref + dV2Aim
    V3Aim = V3Ref + dV3Aim    
        
    print ('V2V3 point %8.3f %8.3f' %(V2Aim, V3Aim))
    A = fulltransform.rotations.attitude(V2Aim, V3Aim, gsRA, gsDec, gspa)
        
    # Get RA and Dec for each segment.
    SegRA = np.zeros(nseg)
    SegDec = np.zeros(nseg) 
    for i in range(nseg):
        V2 = V2Ref + V2Segs[i]
        V3 = V3Ref + V3Segs[i]
        (SegRA[i],SegDec[i]) = fulltransform.rotations.pointing(A, V2, V3)
        print ('%2d %10.4f %10.4f %12.7f %12.7f' %(i+1,V2,V3,SegRA[i], SegDec[i]))
        
    # Plot results
    pl.figure(2)
    pl.clf()
    pl.plot(SegRA, SegDec, 'b*')
    pl.grid(True)
    pl.title('Segment RA and Dec')
    pl.xlabel('RA')
    pl.ylabel('Dec')
    RAmean = SegRA.mean()
    Decmean = SegDec.mean()
    pl.plot(RAmean, Decmean, 'ro')
    if segN > 0: pl.plot(SegRA[segN-1], SegDec[segN-1], 'bx', markersize=12)        
    for i in range(nseg):
        pl.text(SegRA[i],SegDec[i], str(i+1))
    pl.gca().invert_xaxis()
    pl.gca().ticklabel_format(useOffset=False)    
    pl.savefig('SegmentSky.png')
    
    # Other output
    # Pixel calculations
    #print ('Aperture parameters')
    XDetRef = float(siafData['XDetRef'][r])
    YDetRef = float(siafData['YDetRef'][r])
    XSciRef = float(siafData['XSciRef'][r])
    YSciRef = float(siafData['YSciRef'][r])    
    order = int(siafData['Sci2IdlDeg'][r])
    
    # From xml tree read Idl2Sci coefficients  
    terms = (order+1)*(order+2)//2  # Integer division
    C = np.zeros(terms)
    D = np.zeros(terms)
    k = 0
    for i in range(order+1):
        for j in range(i+1):
            suffix = str(i) + str(j)
            colName = 'Idl2SciX' + suffix # e,g 'Idl2SciX21'
            value = siafData[colName][r]
            C[k] = float(value)
            colName = 'Idl2SciY' + suffix
            value = siafData[colName][r]
            D[k] = float(value)
            k += 1

    #print ('C')
    #polynomial.triangle(C,order)
    #print ('D')
    #polynomial.triangle(D,order) 

    xDet = np.zeros(nseg)
    yDet = np.zeros(nseg)
    print('\nSegment x y positions')
    for p in range(nseg):
        # for OSS apertures Sci same as Det
        xDet[p] = XSciRef + fulltransform.polynomial.poly(C, xIdlSegs[p], yIdlSegs[p], order)
        yDet[p] = YSciRef + fulltransform.polynomial.poly(D, xIdlSegs[p], yIdlSegs[p], order)
        print ('%8s %8.2f %8.2f' %(SegIDs[p], xDet[p], yDet[p]))
        if xDet[p] < 0.5 or xDet[p] > 2048.5: print ('%8s off detector in X direction' %SegID[p])       
        if yDet[p] < 0.5 or yDet[p] > 2048.5: print ('%8s off detector in Y direction' %SegID[p])

    # Retrieve FGS parameters obtained in FGSsetup 
    fgsN = sgForm.fgsNum.get()
    V2Ref = float(sgForm.fgsV2)
    V3Ref = float(sgForm.fgsV3)
    V3IdlYAngle = float(sgForm.fgsAngle)

    
    # Summary output
    print ('\nSummary')
    print ('Aperture FGS', fgsN)
    print ('V2Ref %8.3f V3Ref %8.3f IdlAngle %8.4f' %(V2Ref, V3Ref, V3IdlYAngle))
    print ('Guide star at RA %10.6f  Dec %10.6f' %(gsRA, gsDec))
    print ('Position angle %8.4f' %gspa)
    print ('\nSegment    dV2    dV3    xIdl   yIdl     RA         Dec         xDet     yDet')


    for p in range(nseg):
        print ('%5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f' \
        %(SegIDs[p], V2Segs[p], V3Segs[p], xIdlSegs[p], yIdlSegs[p], SegRA[p], SegDec[p], xDet[p], yDet[p]))

    # Final output
    print ('\nFinal Output')
    sg = open('SegmentGuiding.txt', 'w')
    rate = 0.0  # placeholder for count rate
    for p in range(nseg):
        part1 = '-star%02d = %12.6f %12.6f %8.2f  ' %(p+1, SegRA[p], SegDec[p], rate)
        print (part1, end='')
        sg.write(part1)
        onDet = [] #List of other segments on detector
        for q in range(nseg):
            if (q != p) and (1 <= xDet[q] <= 2048) and (1 <= yDet[q] <= 2048) : onDet.append(q+1)
        ns = len(onDet)
        for q in range(ns-1): 
            part2 = '%2d,' % onDet[q]
            print (part2, end=' ')
            sg.write(part2)
        part3 = '%2d' %onDet[ns-1]    
        print (part3)
        sg.write(part3 + '\n')

    sg.close()    
    return

def showsegments():
    # Segments
    # Currently read from a local file V2V3offsets.txt.
    # Might later get from PPS database
    segments = ascii.read('V2V3offsets.txt')
    print ('Segment data read from V2V3offsets.txt')
    SegID = segments['SegID']
    V2Seg = segments['V2Seg']
    V3Seg = segments['V3Seg'] 
    V2mean = V2Seg.mean()
    V3mean = V3Seg.mean()
    nseg = len(V2Seg)
    print (nseg, ' Segments')

    pl.figure(1)
    pl.clf()
    pl.plot(V2Seg,V3Seg, 'b*')
    pl.plot(V2mean, V3mean, 'ro')
    pl.grid(True)
    pl.axis([40.0,-40.0,-40.0, 40.0]) # V2 to the left
    pl.title('Segments')
    pl.xlabel('<-- Delta V2')
    pl.ylabel('Delta V3 -->')
    for s in range(nseg):
        pl.text(V2Seg[s],V3Seg[s], SegID[s])
    pl.savefig('Segments.png')
    return (SegID,V2Seg, V3Seg)
    
    
def checkout(str, low, high):
    """Test conversion from string to float.
    if OK test range
    return errcode 0 for OK, 1 for conversion error, 2 for range error"""
    
    try:
        x = float(str)
        if low <= x <= high: errcode = 0
        else: errcode = 2
    except ValueError:
        errcode = 1
    return errcode        

############################## Main Program - Actions #####################################

sgForm = SegmentForm()
sgForm.CreateForm()
sgForm.Show()



