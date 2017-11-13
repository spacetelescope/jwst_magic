#!/usr/bin/env python
# encoding: utf-8

'''SegmentGuiding.py
Used with guide stars before mirror segments are aligned'''

# Setup
from math import *
import numpy as np
from astropy.io import ascii
import glob
import matplotlib.pyplot as pl
import xml.etree.cElementTree as etree

# Private packages 
import fulltransform.rotations
import fulltransform.polynomial


print ('\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n')
print ('FGS Geometry Tool for unaligned mirror segments')

# Segments
# Currently read from a local file V2V3offsets.txt.
# Might later get from PPS database
segments = ascii.read('V2V3offsets.txt')
print ('Segment data read from V2V3offsets.txt')
SegID = segments['SegID']
V2seg = segments['V2seg']
V3seg = segments['V3seg']
V2mean = V2seg.mean()
V3mean = V3seg.mean()
nseg = len(V2seg)
print (nseg, ' Segments')

pl.figure(1)
pl.clf()
pl.plot(V2seg,V3seg, 'b*')
pl.plot(V2mean, V3mean, 'ro')
pl.grid(True)
pl.axis([80.0,-80.0,-80.0, 80.0]) # V2 to the left
pl.title('Segments')
pl.xlabel('<-- Delta V2')
pl.ylabel('Delta V3 -->')
for s in range(nseg):
    pl.text(V2seg[s],V3seg[s], SegID[s])
pl.savefig('Segments.png')

# Convert to FGS Ideal positions
# Read local FGS xml file

xmlList = glob.glob('/itar/jwst/tel/share/SIAF_WG/Instruments/FGS/FGS_SIAF_20*.xml')
if len(xmlList) < 1:
    print ('No xml file found')
    
xmlList.sort() # Will normally be in date order
print ('Using xml file\n', xmlList[-1]) 
xmlFile = xmlList[-1]  
siaf = etree.parse(xmlFile)
root = siaf.getroot()
siafData = {}               # Dictionary to store all SIAF data

for item in root[0]:
    columnName = item.tag
    columnList = []
    for e in root.iter(columnName): columnList.append(e.text)
    columnArray = np.array(columnList) # Convert to array for easier processing
    siafData[columnName] = columnArray    # Each dictionary column now an array indexed by column name

OK = False # Until proved True
while not OK:
    try:
        fgsN = input('FGS number (1 or 2):  ')
        fgsN = fgsN.strip()
        if fgsN == '1' or fgsN == '2': OK = True
        else: print ('Incorrect FGS number') 
    except: print ('Data input error')
    
det = 'FGS' + fgsN +'_FULL_OSS'
print (det)

row = np.where(siafData['AperName'] == det)
r =  (row[0][0])
print ('Row', r)

V2Ref = float(siafData['V2Ref'][r])
V3Ref = float(siafData['V3Ref'][r])
V3IdlYAngle = float(siafData['V3IdlYAngle'][r])
VIdlParity = int(siafData['VIdlParity'][r])
print ('Reference point %10.4f %10.4f arcsec' %(V2Ref,V3Ref))
print ('Angle %10.4f degrees Parity %2d' %(V3IdlYAngle, VIdlParity))

# FGS calculations
theta = radians(V3IdlYAngle)

xIdlSeg = VIdlParity*(V2seg*cos(theta) - V3seg*sin(theta))
yIdlSeg = V2seg*sin(theta) + V3seg*cos(theta)
print('\n Segment Positions')
print (' Segment     dV2        dV3        XIdl       YIdl')
for p in range(nseg): print ('%8s %10.4f %10.4f %10.4f %10.4f' %(SegID[p],V2seg[p], V3seg[p], xIdlSeg[p], yIdlSeg[p]))
print ('Mean position %8.3f %8.3f\n' %(V2mean,V3mean))

# RA and Dec calculations
# Nominal position for guide star
OK = False # Data input not verified
while not OK:
    try:
        starPos = input('Guide Star RA and Dec (decimal degrees)  ')
        [RA,Dec] = starPos.split()
        a0 = float(RA); d0 = float(Dec)
        if (0.0 <= a0 <= 360.0) and (-90.0 <= d0 <= 90.0): OK = True
        else: print ('Data out of range') 
    except: print ('Data input error')
            
OK = False # Data input not verified
while not OK:
    try:
        pa = float(input('Position Angle  '))
        if (-360.0 <= pa <= 360.0): OK = True
        else: print ('Data out of range')
    except: print ('Data input error')     

OK = False
while not OK:
    try:
        UseSeg = int(input('Number of segment to use (zero to use centroid)  '))
        if (0 <= UseSeg <= 18): OK = True
        else: print ('Data out of range')
    except: print ('Data input error')     

if UseSeg == 0:
    V2Mean = V2seg.mean()
    V3Mean = V3seg.mean()
    V2Aim = V2Ref + V2Mean
    V3Aim = V3Ref + V3Mean
else:
    SegIndex = UseSeg -1   #Allow for zero indexing    
    V2Aim = V2Ref + V2seg[SegIndex]
    V3Aim = V3Ref + V3seg[SegIndex]
    print ('Ref %8.4f %8.4f' %(V2Ref, V3Ref))
    print ('Aim %8.4f %8.4f' %(V2Aim, V3Aim))

# Calculate attitude to place this star at mean segment position
A = fulltransform.rotations.attitude(V2Aim, V3Aim, a0, d0, pa)
#print ('Attitude Matrix')
#print (A)
(a1,d1) = fulltransform.rotations.pointing(A, V2Aim, V3Aim)
# Check we have the right pointing
print ('Point to %8.4f %8.4f'%(a1, d1))

pl.figure(2)
pl.clf()
pl.axis('equal')
pl.plot(a1, d1, 'ro', markersize=10.0)

# Get RA and Dec for each segment

print ('Segment      V2        V3        RA           Dec')
apoint = np.zeros(nseg)
dpoint = np.zeros(nseg)
for s in range(nseg):
    if UseSeg > 0:
        V2center = V2seg[SegIndex]
        V3center = V3seg[SegIndex]
    else:
        V2center = V2Mean
        V3center = V3Mean    
    V2 = V2Aim + V2seg[s] - V2center
    V3 = V3Aim + V3seg[s] - V3center
    (a,d) = fulltransform.rotations.pointing(A, V2, V3)
    apoint[s] = a
    dpoint[s] = d
    print ('%8s %8.4f %8.4f %12.6f %12.6f' %(SegID[s], V2, V3, a, d))
    da = (a-a0)*cos(radians(d))
    pl.plot(a0+da,d, 'b*')
    pl.text(a0+da,d, SegID[s])
# Reverse x axis - RA increasing to the left
u = pl.axis()
v = (u[1],u[0],u[2],u[3])
pl.axis(v)
pl.ticklabel_format(useOffset=False)
pl.axis('equal')
titleString = 'FGS %1d Segment %2d\nTarget (%10.6f %10.6f)   PA %9.4f' %(int(fgsN),UseSeg,a0,d0,pa)
pl.title(titleString)
pl.xlabel('<-- RA')
pl.ylabel('Dec -->')
pl.savefig('RA_Dec.png')


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
    xDet[p] = XSciRef + fulltransform.polynomial.poly(C, xIdlSeg[p], yIdlSeg[p], order)
    yDet[p] = YSciRef + fulltransform.polynomial.poly(D, xIdlSeg[p], yIdlSeg[p], order)
    print ('%8s %8.2f %8.2f' %(SegID[p], xDet[p], yDet[p]))
    if xDet[p] < 0.5 or xDet[p] > 2048.5: print ('%8s off detector in X direction' %SegID[p])       
    if yDet[p] < 0.5 or yDet[p] > 2048.5: print ('%8s off detector in Y direction' %SegID[p])
    
# Summary output
print ('\nSummary')
print ('Aperture FGS', fgsN)
print ('V2Ref %8.3f V3Ref %8.3f IdlAngle %8.4f' %(V2Ref, V3Ref, V3IdlYAngle))
print ('Guide star at RA %10.6f  Dec %10.6f' %(a0, d0))
print ('Position angle %8.4f' %pa)
print ('\nSegment    dV2    dV3    xIdl   yIdl     RA         Dec         xDet     yDet')


for p in range(nseg):
    print ('%5s    %6.2f %6.2f  %6.2f %6.2f  %10.6f %10.6f  %8.2f %8.2f' \
    %(SegID[p], V2seg[p], V3seg[p], xIdlSeg[p], yIdlSeg[p], apoint[p], dpoint[p], xDet[p], yDet[p]))
    
# Final output
print ('\nFinal Output')
sg = open('SegmentGuiding.txt', 'w')
rate = 0.0  # placeholder for count rate
for p in range(nseg):
    part1 = '-star%02d = %12.6f %12.6f %8.2f  ' %(p+1, apoint[p], dpoint[p], rate)
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
        

# Quick test that separations remain constant for any pointing
# To activate uncomment testsep on last line

def sep(a1,d1,a2,d2):
    cossep = cos(radians(a1-a2))*cos(radians(d1))*cos(radians(d2)) + sin(radians(d1))*sin(radians(d2))
    s = 3600*degrees(acos(cossep))
    return s

def testsep():
    print ('\nSeparations calculated from RAs and Decs')
    print ('          ', end='')
    for i in range(nseg-1): print('%11s' %SegID[i], end='')
    print ()
    for i in range(1, nseg):
        print ('%11s  '  %SegID[i], end='')
        for j in range(i):
            s = sep(apoint[i], dpoint[i], apoint[j], dpoint[j])
            print ('%10.4f ' %s, end='')    
    return
    
print ('\nPlots are saved in Segments.png and RA_Dec.png')
print ('Summary data saved in SegmentGuiding.txt\n')
#testsep()
input('Hit return to clear screen images and exit')    

    