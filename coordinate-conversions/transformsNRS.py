from math import *
import numpy as NP
import xml.etree.ElementTree as ET
import os, glob
import polynomial as P
import rotations as R
from tkinter import *

global apPars

def findxml(apName):
    '''Locate latest SIAF xml file containing the aperture'''
    insName = {'FGS':'FGS', 'NRC':'NIRCam','NRS':'NIRSpec', 'NIS':'NIRISS', 'MIR':'MIRI'}
    insFolder = '/itar/jwst/tel/share/SIAF_WG/Instruments/'
    apName = apName.strip()
    try:
        ins = insName[apName[:3]]
        searchString = insFolder + ins +'/' + ins + '_SIAF_20*.xml'
        print (searchString)
        xmlList = glob.glob(searchString)
        files = len(xmlList)
        #print (files, 'xml files')
        if files > 0:
            xmlList.sort()
            xmlFile = xmlList[-1] # Latest version
            xmlSplit = xmlFile.split('/')
            print ('Use', xmlSplit[-1])
            message.configure(text='')
            
            return xmlFile
        else:
            print ('No xml files found')
            message.configure(text='No xml files found')
            
            return
    except KeyError:
        print ('No such instrument key')
        message.configure(text='No such instrument key')
        return
        
def xml(xmlfile):
    '''Reads in xml file version of spreadsheet'''
    siaf = ET.parse(xmlfile)  # This is where the work is done
    root = siaf.getroot()
    siafData = {}               # Dictionary to store all SIAF data

    for item in root[0]:
        columnName = item.tag
        columnList = []
        for e in root.iter(columnName):
            columnList.append(e.text)
            
        columnArray = NP.array(columnList) # Convert to array for easier indexing
        siafData[columnName] = columnArray    # Each dictionary column now an array indexed by column name

    return siafData



def apData(apName,xmlFile):
    siafData = xml(xmlFile)
    order = int(siafData['Sci2IdlDeg'][0])
    print ('Order',order)
    apertures = siafData['AperName']
    print (apName)
    row = NP.where(apertures == apName)[0]
    row = int(row)
    print ('Aperture row', row)
    apTypes = (siafData['AperType']) # The full array
    apType = apTypes[row] # Single instance
    print (apType)
    if apType in ['FULLSCA', 'SUBARRAY', 'ROI']:
        xDetRef = float(siafData['XDetRef'][row])
        yDetRef = float(siafData['YDetRef'][row])
        DetAngle = float(siafData['DetSciYAngle'][row])
        DetPar = int(siafData['DetSciParity'][row])
        xSciRef = float(siafData['XSciRef'][row])
        ySciRef = float(siafData['YSciRef'][row])
        print ('Detector', xDetRef, yDetRef)
        print ('Science ', xSciRef, ySciRef, DetAngle, DetPar)
        order = int(siafData['Sci2IdlDeg'][row])
        
    V2Ref = float(siafData['V2Ref'][row])
    V3Ref = float(siafData['V3Ref'][row])
    IdlAngle = float(siafData['V3IdlYAngle'][row])
    VParity = int(siafData['VIdlParity'][row])
    print ('Ideal   ', V2Ref, V3Ref, IdlAngle, VParity)


    if apType in ['FULLSCA', 'SUBARRAY', 'ROI'] :
        A = []
        B = []
        C = []
        D = []
        k = 0 
        for i in range(order+1):
            for j in range(i+1):
                suffix = str(i)+str(j)
                Aname = 'Sci2IdlX' + suffix
                Bname = 'Sci2IdlY' + suffix
                Cname = 'Idl2SciX' + suffix
                Dname = 'Idl2SciY' + suffix
                v = float(siafData[Aname][row])
                A.append(v)
                v = float(siafData[Bname][row])
                B.append(v)
                v = float(siafData[Cname][row])
                C.append(v)
                v = float(siafData[Dname][row])
                D.append(v)        
    
        if apName[:3] == 'NRS':
            print ('NIRSpec')
            A2=[]
            B2=[]
            C2=[]
            D2=[]
            # Get row CLEAR_GWA_OTE
            row = NP.where(apertures == 'CLEAR_GWA_OTE')[0]
            row = int(row)
            print ('Aperture row', row)

            for i in range(order+1):
                for j in range(i+1):
                    suffix = str(i)+str(j)
                    Aname = 'Sci2IdlX' + suffix
                    Bname = 'Sci2IdlY' + suffix
                    Cname = 'Idl2SciX' + suffix
                    Dname = 'Idl2SciY' + suffix
                    v = float(siafData[Aname][row])
                    A2.append(v)
                    v = float(siafData[Bname][row])
                    B2.append(v)
                    v = float(siafData[Cname][row])
                    C2.append(v)
                    v = float(siafData[Dname][row])
                    D2.append(v)
                
            print ('A2\n',A2) 
            return (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2)
        
        else: # Not NIRSpec
            return (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D)
                       
       
def Finish():
    root.destroy()



global xmlFile,siafData, apList

def apsetup(*args):
    global xmlFile,siafData, apList, apPars
    aperture = apName.get()
    #print (aperture)
    if len(aperture) == 3: # Get instrument and xml file
        xmlFile = findxml(aperture)
        siafData = xml(xmlFile)
        apList = siafData['AperName']
        print (apList)
        
    if len(aperture) > 3:    
        if aperture in apList: # Look for complete aperture name
            apPars = apData(aperture, xmlFile)
            print ('apPars',apPars)
            apType.set(apPars[0])
            xDetRef.set(apPars[1])
            yDetRef.set(apPars[2])
            DetYAngle.set(apPars[3])
            DetParity.set(apPars[4])
            xSciRef.set(apPars[5])
            ySciRef.set(apPars[6])
            IdlAngle.set(apPars[7])
            V2Ref.set(apPars[8])
            V3Ref.set(apPars[9])
            IdlParity.set(apPars[10]) 
            
            # Only when aperture is selected can calculations be done.
            message.configure(text=aperture)
            print ('Type', apPars[0])
            if apPars[0] in ['FULLSCA', 'SUBARRAY', 'ROI']:
                dButton.configure(state='normal')
                sButton.configure(state='normal')
            iButton.configure(state='normal')
            vButton.configure(state='normal')
            
            return apPars

# Transform modules
def Det2Sci():
    """Detector to Science"""
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars
    print ('Det2Sci apType', apType)
    if apType not in ['FULLSCA', 'SUBARRAY', 'ROI']:
        return
    try:
        x = xDet.get()
    except:
        x = 0.0    
    try:
        y = yDet.get()
    except:
        y = 0.0
            
    yFlip = cos(radians(DetAngle))
    xFlip = DetPar*yFlip
    xout = xSciRef + xFlip*(x-xDetRef)
    yout = ySciRef + yFlip*(y-yDetRef)
    xSci.set(round(xout,2))
    ySci.set(round(yout,2))
    return

def Sci2Idl():
    """Science to Ideal"""
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars

    if apType not in ['FULLSCA', 'SUBARRAY', 'ROI']:
        return
    try:
        x = xSci.get() - xSciRef
    except:
        x = 0.0    
    try:
        y = ySci.get() - ySciRef
    except:
        y = 0.0
    xout = P.poly(A,x,y)
    yout = P.poly(B,x,y)
    xIdl.set(round(xout,4)) # Only keep 4 decimal places
    yIdl.set(round(yout,4))   
    return    


def Idl2Sci():
    """Science to Ideal"""
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars

    if apType not in ['FULLSCA', 'SUBARRAY', 'ROI']:
        return
    try:
        x = xIdl.get()
    except:
        x = 0.0    
    try:
        y = yIdl.get()
    except:
        y = 0.0
    xout = P.poly(C,x,y) + xSciRef
    yout = P.poly(D,x,y) + ySciRef
    xSci.set(round(xout,2))
    ySci.set(round(yout,2))    
    return    



def Idl2V():
    '''Ideal to V2V3'''
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars

    try:
        x = xIdl.get()
    except:
        x = 0.0    
    try:
        y = yIdl.get()
    except:
        y = 0.0
    
    a = radians(IdlAngle)    
    v2out = V2Ref + VParity*x*cos(a) + y*sin(a)
    v3out = V3Ref - VParity*x*sin(a) + y*cos(a)
    V2.set(round(v2out,4))
    V3.set(round(v3out,4))
    return    

def V2Idl():
    '''V2V3 to Ideal'''
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars

    try:
        x = V2.get()
    except:
        x = 0.0    
    try:
        y = V3.get()
    except:
        y = 0.0

    a = radians(IdlAngle)
    dx = x - V2Ref
    dy = y - V3Ref
    xout = VParity*(dx*cos(a) - dy*sin(a))
    yout = dx*sin(a) + dy*cos(a)
    xIdl.set(round(xout,4))
    yIdl.set(round(yout,4))
    return    

    
def idealrotate(V2,V3,pa):
    """Rotation matrix for FGS Ideal"""
    r1 = R.rotate(3,-V2/3600)
    r2 = R.rotate(2, V3/3600)
    r3 = R.rotate(1, pa)
    M = NP.dot(r2,r1)
    M = NP.dot(r3,M)
    return M
    

def V2FGS():
    """From V2V3 to FGS1 and 2 Ideal"""
    # Make up rotation matrices
    v21 = f1V2Ref.get()
    v31 = f1V3Ref.get()
    pa = f1V3Angle.get()
    M1 = idealrotate(v21,v31,pa)
    
    v22 = f2V2Ref.get()
    v32 = f2V3Ref.get()
    pa = f2V3Angle.get()
    M2 = idealrotate(v22,v32,pa)
    
    v2 = V2.get()
    v3 = V3.get()
    u = R.unit(v2/3600, v3/3600)
    
    # FGS1
    v = NP.dot(M1,u)
    (w1,w2) = R.v2v3(v)
    P = f1Parity.get()
    w1 = P*w1
    f1xIdl.set(round(w1,4))
    f1yIdl.set(round(w2,4))
    
    # FGS2
    v = NP.dot(M2,u)
    (w1,w2) = R.v2v3(v)
    P = f2Parity.get()
    w1 = P*w1
    f2xIdl.set(round(w1,4))
    f2yIdl.set(round(w2,4))
    
    return
    
def Sci2Det():    
    """Science to Detector"""
    global apPars
    if len(apPars) == 16:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D) = apPars
    else:
        (apType,xDetRef, yDetRef, DetAngle, DetPar, xSciRef, ySciRef, IdlAngle, V2Ref, V3Ref, VParity, order, A,B,C,D, A2,B2,C2,D2) = apPars
                
    try:
        x = xSci.get()
    except:
        x = 0.0    
    try:
        y = ySci.get()
    except:
        y = 0.0
        
    yFlip = cos(radians(DetAngle))
    xFlip = DetPar*yFlip
    xout = xDetRef + xFlip*(x-xSciRef)
    yout = yDetRef + yFlip*(y-ySciRef)
    xDet.set(round(xout,2))
    yDet.set(round(yout,2))
    return

def textout():
    """Output results"""
    convName = 'conversions.txt'
    outString = 12*'%10.4f' %(xDet.get(), yDet.get(),xSci.get(), ySci.get(),xIdl.get(), yIdl.get(),
                            V2.get(), V3.get(), f1xIdl.get(), f1yIdl.get(), f2xIdl.get(), f2yIdl.get())
    print (outString)                        
    if os.path.exists(convName): # Make new file and write header row
        print ('Append to existing file')
        conv = open(convName,'a')
        conv.write(outString + '\n')

    else:    
        print ('New text file')
        conv = open(convName, 'w')
        conv.write('   Xdet      YDet      Xsci      Ysci      Xidl      YIdl        V2        V3      F1XIdl    F1YIdl    F2XIdl    F2YIdl \n')
        conv.write(outString + '\n')
    
    conv.close()
    return       
    
# Transform combinations
def Detector():
    Det2Sci()
    Sci2Idl()
    Idl2V()
    V2FGS()
    textout()
    return
    
def Science():
    Sci2Det()
    Sci2Idl()
    Idl2V()
    V2FGS()
    textout()
    return

def Ideal():
    Idl2Sci()
    Sci2Det()
    Idl2V()
    V2FGS()
    textout()
    return
    
def V2V3():
    V2FGS()
    V2Idl()
    Idl2Sci()
    Sci2Det()
    textout()
    return            

# Interface layout

root = Tk()
root.title("Colin Cox's Coordinate Conversions")
shade = 'moccasin'
root['bg'] = shade

# Aperture Entry
Label(root, text='Aperture name', bg=shade).grid(row=0)
apName = StringVar()
Entry(root, textvariable = apName).grid(row=0,column=1)
apName.trace('w', apsetup)
apType = StringVar() 
Entry(root, textvariable=apType,bg=shade).grid(row=0,column=2)
message = Label(root, text='message', bg=shade)
message.grid(row=0, column=3, columnspan=3)

# Titles
Label(root,text='X', bg=shade).grid(row=1, column=1)
Label(root,text='Y', bg=shade).grid(row=1, column=2)
Label(root,text='XRef', bg=shade).grid(row=1, column=4)
Label(root,text='YRef', bg=shade).grid(row=1, column=5)
Label(root, text='Angle', bg=shade).grid(row=1,column=6)
Label(root, text='Parity', bg=shade).grid(row=1,column=7)

#Detector frame
Label(root, text='Detector', bg=shade).grid(row=2)
xDet = DoubleVar()
yDet = DoubleVar()
xDetRef = DoubleVar()
yDetRef = DoubleVar()
DetYAngle = DoubleVar()
DetParity = IntVar()
Entry(root, textvariable=xDet).grid(row=2,column=1)
Entry(root, textvariable=yDet).grid(row=2,column=2)
dButton = Button(root, text='GO', fg='red', command=Detector, state='disabled')
dButton.grid(row=2,column=3)
Entry(root, textvariable=xDetRef, bg='lightgrey').grid(row=2,column=4)
Entry(root, textvariable=yDetRef, bg='lightgrey').grid(row=2,column=5)
Entry(root, textvariable=DetYAngle, bg='lightgrey').grid(row=2,column=6)
Entry(root, textvariable=DetParity, bg='lightgrey').grid(row=2,column=7)


#Science frame
Label(root, text='Science', bg=shade).grid(row=3)
xSci = DoubleVar()
ySci = DoubleVar()
xSciRef = DoubleVar()
ySciRef = DoubleVar()
Entry(root, textvariable=xSci).grid(row=3,column=1)
Entry(root, textvariable=ySci).grid(row=3,column=2)
sButton = Button(root, text='GO', command=Science, state='disabled')
sButton.grid(row=3,column=3)
Entry(root, textvariable=xSciRef, bg='lightgrey').grid(row=3,column=4)
Entry(root, textvariable=ySciRef, bg='lightgrey').grid(row=3,column=5)

#Ideal Frame
Label(root, text='Ideal', bg=shade).grid(row=4)
xIdl = DoubleVar()
yIdl = DoubleVar()
V2Ref = DoubleVar()
V3Ref = DoubleVar()
IdlAngle = DoubleVar()
IdlParity = DoubleVar()
Entry(root, textvariable=xIdl).grid(row=4,column=1)
Entry(root, textvariable=yIdl).grid(row=4,column=2)
iButton = Button(root, text='GO', command=Ideal, state='disabled')
iButton.grid(row=4,column=3)
Entry(root, textvariable=V2Ref, bg='lightgrey').grid(row=4,column=4)
Entry(root, textvariable=V3Ref, bg='lightgrey').grid(row=4,column=5)
Entry(root, textvariable=IdlAngle, bg='lightgrey').grid(row=4,column=6)
Entry(root, textvariable=IdlParity, bg='lightgrey').grid(row=4,column=7)

#V2V3 frame
Label(root, text='V2 V3', bg=shade).grid(row=5)
V2 = DoubleVar()
V3 = DoubleVar()
Entry(root, textvariable=V2).grid(row=5,column=1)
Entry(root, textvariable=V3).grid(row=5,column=2)
vButton = Button(root, text='GO', command=V2V3, state='disabled')
vButton.grid(row=5,column=3)

# FGS Ideal
f1xIdl = DoubleVar()
f1yIdl = DoubleVar()
f1V2Ref = DoubleVar()
f1V3Ref = DoubleVar()
f1V3Angle = DoubleVar()
f1Parity = IntVar()
Label(root, text='FGS1 Ideal', bg=shade).grid(row=6)
Entry(root, textvariable=f1xIdl, bg='lightgrey').grid(row=6,column=1)
Entry(root, textvariable=f1yIdl, bg='lightgrey').grid(row=6,column=2)
Entry(root, textvariable=f1V2Ref, bg='lightgrey').grid(row=6,column=4)
Entry(root, textvariable=f1V3Ref, bg='lightgrey').grid(row=6,column=5)
Entry(root, textvariable=f1V3Angle, bg='lightgrey').grid(row=6,column=6)
Entry(root, textvariable=f1Parity, bg='lightgrey').grid(row=6,column=7)

f2xIdl = DoubleVar()
f2yIdl = DoubleVar()
f2V2Ref = DoubleVar()
f2V3Ref = DoubleVar()
f2V3Angle = DoubleVar()
f2Parity = IntVar()
Label(root, text='FGS2 Ideal', bg=shade).grid(row=7)
Entry(root, textvariable=f2xIdl, bg='lightgrey').grid(row=7,column=1)
Entry(root, textvariable=f2yIdl, bg='lightgrey').grid(row=7,column=2)
Entry(root, textvariable=f2V2Ref, bg='lightgrey').grid(row=7,column=4)
Entry(root, textvariable=f2V3Ref, bg='lightgrey').grid(row=7,column=5)
Entry(root, textvariable=f2V3Angle, bg='lightgrey').grid(row=7,column=6)
Entry(root, textvariable=f2Parity, bg='lightgrey').grid(row=7,column=7)


Button(root, text='Finish', command=Finish, bg='red').grid(row=8,column=3)

# End of layout definition

# Action starts here

# Install FGS data
xmlFile = findxml('FGS1_FULL')
FGS1par = apData('FGS1_FULL', xmlFile)
print ('FGS1par', FGS1par)
f1V2Ref.set(FGS1par[8])
f1V3Ref.set(FGS1par[9])
f1V3Angle.set(FGS1par[7])
f1Parity.set(FGS1par[10])

FGS2par = apData('FGS2_FULL', xmlFile) # Use same xml file
f2V2Ref.set(FGS2par[8])
f2V3Ref.set(FGS2par[9])
f2V3Angle.set(FGS2par[7])
f2Parity.set(FGS2par[10])

message.configure(text='Enter aperture name')
root.mainloop()

    
    