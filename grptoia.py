import numpy as np

#real pixel to ideal angle

#xOffset = 1023.5
#yOffset = 1023.5
#plateScale = 0.06738281367
#gain=1.5


## read in X Y in pixels

def g1RPtoIA(yRealPixel,xRealPixel,xOffset=1023.5, yOffset=1023.5,
             plateScale=0.06738281367, gain=1.5):

    # the following coeffs are from julia
    # dated 17-Feb-2012
    realPToIdealPXCoeff0 = -2.3132211E+01
    realPToIdealPXCoeff1 = 9.9858571E-01
    realPToIdealPXCoeff2 = 1.0458177E-02
    realPToIdealPXCoeff3 = 2.6914738E-06
    realPToIdealPXCoeff4 = 6.7167416E-06
    realPToIdealPXCoeff5 = 9.9063452E-07
    realPToIdealPXCoeff6 = 1.2100144E-09
    realPToIdealPXCoeff7 = -3.4359146E-11
    realPToIdealPXCoeff8 = 1.2718376E-09
    realPToIdealPXCoeff9 = -2.0353075E-11
    realPToIdealPXCoeff10 = 5.3149090E-14
    realPToIdealPXCoeff11 = 9.3076979E-14
    realPToIdealPXCoeff12 = 8.3907183E-14
    realPToIdealPXCoeff13 = 9.7074620E-14
    realPToIdealPXCoeff14 = 1.9366816E-14

    realPToIdealPYCoeff0 = -2.6649388E+01
    realPToIdealPYCoeff1 = -2.5585494E-03
    realPToIdealPYCoeff2 = 1.0114078E+00
    realPToIdealPYCoeff3 = 2.4056198E-06
    realPToIdealPYCoeff4 = 2.0788818E-06
    realPToIdealPYCoeff5 = 9.3110110E-06
    realPToIdealPYCoeff6 = -1.6052215E-11
    realPToIdealPYCoeff7 = 1.2730773E-09
    realPToIdealPYCoeff8 = -9.1059163E-11
    realPToIdealPYCoeff9 = 1.3291468E-09
    realPToIdealPYCoeff10 = 2.1479653E-14
    realPToIdealPYCoeff11 = 5.2933520E-14
    realPToIdealPYCoeff12 = 1.4401378E-13
    realPToIdealPYCoeff13 = 7.4524608E-14
    realPToIdealPYCoeff14 = 1.2816500E-13


    # Reverse calculation of the coordinate conversion
    # use shifted coords for these calculations
    xIdealPixel = realPToIdealPXCoeff0
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff1*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff2*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff3*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff4*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff5*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff6*xRealPixel*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff7*xRealPixel*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff8*xRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff9*yRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff10*xRealPixel*xRealPixel*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff11*xRealPixel*xRealPixel*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff12*xRealPixel*xRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff13*xRealPixel*yRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff14*yRealPixel*yRealPixel*yRealPixel*yRealPixel

    yIdealPixel = realPToIdealPYCoeff0
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff1*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff2*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff3*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff4*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff5*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff6*xRealPixel*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff7*xRealPixel*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff8*xRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff9*yRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff10*xRealPixel*xRealPixel*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff11*xRealPixel*xRealPixel*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff12*xRealPixel*xRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff13*xRealPixel*yRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff14*yRealPixel*yRealPixel*yRealPixel*yRealPixel

    xAngle = (xIdealPixel - xOffset)*plateScale
    yAngle = (yOffset - yIdealPixel)*plateScale

    return xAngle, yAngle


def g2RPtoIA(yRealPixel,xRealPixel,xOffset=1023.5, yOffset=1023.5,
             plateScale=0.06738281367, gain=1.5):

    tint=0.338

    # the following coeffs are from julia
    # dated 17-Feb-2012
    realPToIdealPXCoeff0 = -3.2653125E+01
    realPToIdealPXCoeff1 = 1.0343933E+00
    realPToIdealPXCoeff2 = 2.1128135E-02
    realPToIdealPXCoeff3 = -9.0969494E-06
    realPToIdealPXCoeff4 = -1.4320516E-05
    realPToIdealPXCoeff5 = -3.9403361E-06
    realPToIdealPXCoeff6 = 1.6122843E-09
    realPToIdealPXCoeff7 = 5.9346314E-10
    realPToIdealPXCoeff8 = 2.0380335E-09
    realPToIdealPXCoeff9 = 2.0974142E-10
    realPToIdealPXCoeff10 = -2.8502771E-14
    realPToIdealPXCoeff11 = -8.3569450E-14
    realPToIdealPXCoeff12 = -5.0976659E-14
    realPToIdealPXCoeff13 = -9.1800561E-14
    realPToIdealPXCoeff14 = -8.9563271E-15

    realPToIdealPYCoeff0 = -7.7165127E+01
    realPToIdealPYCoeff1 = 2.9152735E-02
    realPToIdealPYCoeff2 = 1.0773690E+00
    realPToIdealPYCoeff3 = -6.2761496E-06
    realPToIdealPYCoeff4 = -7.8529777E-06
    realPToIdealPYCoeff5 = -2.1370284E-05
    realPToIdealPYCoeff6 = 2.0146586E-10
    realPToIdealPYCoeff7 = 2.0163576E-09
    realPToIdealPYCoeff8 = 6.6916423E-10
    realPToIdealPYCoeff9 = 2.4400761E-09
    realPToIdealPYCoeff10 = -1.9936426E-14
    realPToIdealPYCoeff11 = -3.2487010E-14
    realPToIdealPYCoeff12 = -1.3263445E-13
    realPToIdealPYCoeff13 = -4.2058076E-14
    realPToIdealPYCoeff14 = -1.2323597E-13

    #  Reverse calculation of the coordinate conversion
    xIdealPixel = realPToIdealPXCoeff0
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff1*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff2*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff3*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff4*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff5*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff6*xRealPixel*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff7*xRealPixel*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff8*xRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff9*yRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff10*xRealPixel*xRealPixel*xRealPixel*xRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff11*xRealPixel*xRealPixel*xRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff12*xRealPixel*xRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff13*xRealPixel*yRealPixel*yRealPixel*yRealPixel
    xIdealPixel = xIdealPixel + realPToIdealPXCoeff14*yRealPixel*yRealPixel*yRealPixel*yRealPixel

    yIdealPixel = realPToIdealPYCoeff0
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff1*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff2*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff3*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff4*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff5*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff6*xRealPixel*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff7*xRealPixel*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff8*xRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff9*yRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff10*xRealPixel*xRealPixel*xRealPixel*xRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff11*xRealPixel*xRealPixel*xRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff12*xRealPixel*xRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff13*xRealPixel*yRealPixel*yRealPixel*yRealPixel
    yIdealPixel = yIdealPixel + realPToIdealPYCoeff14*yRealPixel*yRealPixel*yRealPixel*yRealPixel

    xAngle = (xIdealPixel - xOffset)*plateScale
    yAngle = (yIdealPixel - yOffset)*plateScale

    return xAngle, yAngle


def write_to_file(xAngle, yAngle):

    with open('ideal.tmp','w') as f:
        for i in range(len(xAngle)):
           f.writeto(i-1,xAngle[i],yAngle[i])#'(i,2x,f10.4,2x,f10.4,f,f)')
