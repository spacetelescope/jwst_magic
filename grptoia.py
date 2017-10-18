import numpy as np

## pass in X Y in pixels
def g1RPtoIA(x,y,xOffset=1023.5, yOffset=1023.5,
             plateScale=0.06738281367, gain=1.5):

    yrealpixel, xrealpixel = x, y # THIS IS NECESSARY, DO NOT QUESTION

    # the following coeffs are from julia
    # dated 17-Feb-2012
    real_p_to_ideal_p_xcoeff0 = -2.3132211E+01
    real_p_to_ideal_p_xcoeff1 = 9.9858571E-01
    real_p_to_ideal_p_xcoeff2 = 1.0458177E-02
    real_p_to_ideal_p_xcoeff3 = 2.6914738E-06
    real_p_to_ideal_p_xcoeff4 = 6.7167416E-06
    real_p_to_ideal_p_xcoeff5 = 9.9063452E-07
    real_p_to_ideal_p_xcoeff6 = 1.2100144E-09
    real_p_to_ideal_p_xcoeff7 = -3.4359146E-11
    real_p_to_ideal_p_xcoeff8 = 1.2718376E-09
    real_p_to_ideal_p_xcoeff9 = -2.0353075E-11
    real_p_to_ideal_p_xcoeff10 = 5.3149090E-14
    real_p_to_ideal_p_xcoeff11 = 9.3076979E-14
    real_p_to_ideal_p_xcoeff12 = 8.3907183E-14
    real_p_to_ideal_p_xcoeff13 = 9.7074620E-14
    real_p_to_ideal_p_xcoeff14 = 1.9366816E-14

    real_p_to_ideal_p_ycoeff0 = -2.6649388E+01
    real_p_to_ideal_p_ycoeff1 = -2.5585494E-03
    real_p_to_ideal_p_ycoeff2 = 1.0114078E+00
    real_p_to_ideal_p_ycoeff3 = 2.4056198E-06
    real_p_to_ideal_p_ycoeff4 = 2.0788818E-06
    real_p_to_ideal_p_ycoeff5 = 9.3110110E-06
    real_p_to_ideal_p_ycoeff6 = -1.6052215E-11
    real_p_to_ideal_p_ycoeff7 = 1.2730773E-09
    real_p_to_ideal_p_ycoeff8 = -9.1059163E-11
    real_p_to_ideal_p_ycoeff9 = 1.3291468E-09
    real_p_to_ideal_p_ycoeff10 = 2.1479653E-14
    real_p_to_ideal_p_ycoeff11 = 5.2933520E-14
    real_p_to_ideal_p_ycoeff12 = 1.4401378E-13
    real_p_to_ideal_p_ycoeff13 = 7.4524608E-14
    real_p_to_ideal_p_ycoeff14 = 1.2816500E-13


    # Reverse calculation of the coordinate conversion
    # use shifted coords for these calculations
    xidealpixel = real_p_to_ideal_p_xcoeff0
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff1*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff2*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff3*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff4*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff5*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff6*xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff7*xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff8*xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff9*yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff10*xrealpixel*xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff11*xrealpixel*xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff12*xrealpixel*xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff13*xrealpixel*yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff14*yrealpixel*yrealpixel*yrealpixel*yrealpixel

    yidealpixel = real_p_to_ideal_p_ycoeff0
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff1*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff2*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff3*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff4*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff5*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff6*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff7*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff8*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff9*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff10*xrealpixel*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff11*xrealpixel*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff12*xrealpixel*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff13*xrealpixel*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff14*yrealpixel*yrealpixel*yrealpixel*yrealpixel

    xangle = (xidealpixel - xOffset)*plateScale
    yangle = (yOffset - yidealpixel)*plateScale

    return xangle, yangle


def g2RPtoIA(x, y, xOffset=1023.5, yOffset=1023.5,
             plateScale=0.06738281367, gain=1.5):

    yrealpixel, xrealpixel = x, y # THIS IS NECESSARY, DO NOT QUESTION
    # the following coeffs are from julia
    # dated 17-Feb-2012
    real_p_to_ideal_p_xcoeff0 = -3.2653125E+01
    real_p_to_ideal_p_xcoeff1 = 1.0343933E+00
    real_p_to_ideal_p_xcoeff2 = 2.1128135E-02
    real_p_to_ideal_p_xcoeff3 = -9.0969494E-06
    real_p_to_ideal_p_xcoeff4 = -1.4320516E-05
    real_p_to_ideal_p_xcoeff5 = -3.9403361E-06
    real_p_to_ideal_p_xcoeff6 = 1.6122843E-09
    real_p_to_ideal_p_xcoeff7 = 5.9346314E-10
    real_p_to_ideal_p_xcoeff8 = 2.0380335E-09
    real_p_to_ideal_p_xcoeff9 = 2.0974142E-10
    real_p_to_ideal_p_xcoeff10 = -2.8502771E-14
    real_p_to_ideal_p_xcoeff11 = -8.3569450E-14
    real_p_to_ideal_p_xcoeff12 = -5.0976659E-14
    real_p_to_ideal_p_xcoeff13 = -9.1800561E-14
    real_p_to_ideal_p_xcoeff14 = -8.9563271E-15

    real_p_to_ideal_p_ycoeff0 = -7.7165127E+01
    real_p_to_ideal_p_ycoeff1 = 2.9152735E-02
    real_p_to_ideal_p_ycoeff2 = 1.0773690E+00
    real_p_to_ideal_p_ycoeff3 = -6.2761496E-06
    real_p_to_ideal_p_ycoeff4 = -7.8529777E-06
    real_p_to_ideal_p_ycoeff5 = -2.1370284E-05
    real_p_to_ideal_p_ycoeff6 = 2.0146586E-10
    real_p_to_ideal_p_ycoeff7 = 2.0163576E-09
    real_p_to_ideal_p_ycoeff8 = 6.6916423E-10
    real_p_to_ideal_p_ycoeff9 = 2.4400761E-09
    real_p_to_ideal_p_ycoeff10 = -1.9936426E-14
    real_p_to_ideal_p_ycoeff11 = -3.2487010E-14
    real_p_to_ideal_p_ycoeff12 = -1.3263445E-13
    real_p_to_ideal_p_ycoeff13 = -4.2058076E-14
    real_p_to_ideal_p_ycoeff14 = -1.2323597E-13

    #  Reverse calculation of the coordinate conversion
    xidealpixel = real_p_to_ideal_p_xcoeff0
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff1*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff2*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff3*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff4*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff5*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff6*xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff7*xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff8*xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff9*yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff10*xrealpixel*xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff11*xrealpixel*xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff12*xrealpixel*xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff13*xrealpixel*yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + real_p_to_ideal_p_xcoeff14*yrealpixel*yrealpixel*yrealpixel*yrealpixel

    yidealpixel = real_p_to_ideal_p_ycoeff0
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff1*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff2*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff3*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff4*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff5*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff6*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff7*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff8*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff9*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff10*xrealpixel*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff11*xrealpixel*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff12*xrealpixel*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff13*xrealpixel*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + real_p_to_ideal_p_ycoeff14*yrealpixel*yrealpixel*yrealpixel*yrealpixel

    xangle = (xidealpixel - xOffset)*plateScale
    yangle = (yidealpixel - yOffset)*plateScale

    return xangle, yangle

def RPtoIA(guider, yrealpixel, xrealpixel, xOffset=1023.5, yOffset=1023.5,
           plateScale=0.06738281367, gain=1.5):

    if guider == 1:
        g = get_g1_conversion_factors()
        guider_scale = 1
    elif guider == 2:
        g = get_g2_conversion_factors()
        guider_scale = -1

    #  Reverse calculation of the coordinate conversion
    xidealpixel = g['X']['Coeff0']
    xidealpixel = xidealpixel + g['X']['Coeff1'] * xrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff2'] * yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff3'] * xrealpixel**2 #xrealpixel*xrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff4'] * xrealpixel * yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff5'] * yrealpixel**2 #yrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff6'] * xrealpixel**3 #xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff7'] * xrealpixel**2 * yrealpixel#xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff8'] * xrealpixel * yrealpixel**2 #xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff9'] * yrealpixel**3 #yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff10'] * xrealpixel**4 #xrealpixel*xrealpixel*xrealpixel*xrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff11'] * xrealpixel**3 * yrealpixel #xrealpixel*xrealpixel*xrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff12'] * xrealpixel**2 * yrealpixel**2 #xrealpixel*xrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff13'] * xrealpixel * yrealpixel**3 #xrealpixel*yrealpixel*yrealpixel*yrealpixel
    xidealpixel = xidealpixel + g['X']['Coeff14'] * yrealpixel**4 #yrealpixel*yrealpixel*yrealpixel*yrealpixel

    yidealpixel = g['Y']['Coeff0']
    yidealpixel = yidealpixel + g['Y']['Coeff1']*xrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff2']*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff3']*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff4']*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff5']*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff6']*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff7']*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff8']*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff9']*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff10']*xrealpixel*xrealpixel*xrealpixel*xrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff11']*xrealpixel*xrealpixel*xrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff12']*xrealpixel*xrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff13']*xrealpixel*yrealpixel*yrealpixel*yrealpixel
    yidealpixel = yidealpixel + g['Y']['Coeff14']*yrealpixel*yrealpixel*yrealpixel*yrealpixel

    xangle = (xidealpixel - xOffset)*plateScale
    yangle = (yidealpixel - yOffset)*plateScale*guider_scale

    return xangle, yangle

def get_g1_conversion_factors():
    g1 = {}
    g1['X'] = {}
    g1['X']['Coeff0'] = -2.3132211E+01
    g1['X']['Coeff1'] = 9.9858571E-01
    g1['X']['Coeff2'] = 1.0458177E-02
    g1['X']['Coeff3'] = 2.6914738E-06
    g1['X']['Coeff4'] = 6.7167416E-06
    g1['X']['Coeff5'] = 9.9063452E-07
    g1['X']['Coeff6'] = 1.2100144E-09
    g1['X']['Coeff7'] = -3.4359146E-11
    g1['X']['Coeff8'] = 1.2718376E-09
    g1['X']['Coeff9'] = -2.0353075E-11
    g1['X']['Coeff10'] = 5.3149090E-14
    g1['X']['Coeff11'] = 9.3076979E-14
    g1['X']['Coeff12'] = 8.3907183E-14
    g1['X']['Coeff13'] = 9.7074620E-14
    g1['X']['Coeff14'] = 1.9366816E-14

    g1['Y'] = {}
    g1['Y']['Coeff0'] = -2.6649388E+01
    g1['Y']['Coeff1'] = -2.5585494E-03
    g1['Y']['Coeff2'] = 1.0114078E+00
    g1['Y']['Coeff3'] = 2.4056198E-06
    g1['Y']['Coeff4'] = 2.0788818E-06
    g1['Y']['Coeff5'] = 9.3110110E-06
    g1['Y']['Coeff6'] = -1.6052215E-11
    g1['Y']['Coeff7'] = 1.2730773E-09
    g1['Y']['Coeff8'] = -9.1059163E-11
    g1['Y']['Coeff9'] = 1.3291468E-09
    g1['Y']['Coeff10'] = 2.1479653E-14
    g1['Y']['Coeff11'] = 5.2933520E-14
    g1['Y']['Coeff12'] = 1.4401378E-13
    g1['Y']['Coeff13'] = 7.4524608E-14
    g1['Y']['Coeff14'] = 1.2816500E-13

    return g1


def get_g2_conversion_factors():
    g2 = {}
    g2['X'] = {}
    g2['X']['Coeff0'] = -3.2653125E+01
    g2['X']['Coeff1']  = 1.0343933E+00
    g2['X']['Coeff2']  = 2.1128135E-02
    g2['X']['Coeff3']  = -9.0969494E-06
    g2['X']['Coeff4']  = -1.4320516E-05
    g2['X']['Coeff5']  = -3.9403361E-06
    g2['X']['Coeff6']  = 1.6122843E-09
    g2['X']['Coeff7']  = 5.9346314E-10
    g2['X']['Coeff8']  = 2.0380335E-09
    g2['X']['Coeff9'] = 2.0974142E-10
    g2['X']['Coeff10'] = -2.8502771E-14
    g2['X']['Coeff11'] = -8.3569450E-14
    g2['X']['Coeff12'] = -5.0976659E-14
    g2['X']['Coeff13'] = -9.1800561E-14
    g2['X']['Coeff14'] = -8.9563271E-15

    g2['Y'] = {}
    g2['Y']['Coeff0'] = -7.7165127E+01
    g2['Y']['Coeff1'] = 2.9152735E-02
    g2['Y']['Coeff2'] = 1.0773690E+00
    g2['Y']['Coeff3'] = -6.2761496E-06
    g2['Y']['Coeff4'] = -7.8529777E-06
    g2['Y']['Coeff5'] = -2.1370284E-05
    g2['Y']['Coeff6'] = 2.0146586E-10
    g2['Y']['Coeff7'] = 2.0163576E-09
    g2['Y']['Coeff8'] = 6.6916423E-10
    g2['Y']['Coeff9'] = 2.4400761E-09
    g2['Y']['Coeff10'] = -1.9936426E-14
    g2['Y']['Coeff11'] = -3.2487010E-14
    g2['Y']['Coeff12'] = -1.3263445E-13
    g2['Y']['Coeff13'] = -4.2058076E-14
    g2['Y']['Coeff14'] = -1.2323597E-13

    return g2


def write_to_file(xangle, yangle):

    with open('ideal.tmp','w') as f:
        for i in range(len(xangle)):
           f.writeto(i-1,xangle[i],yangle[i])#'(i,2x,f10.4,2x,f10.4,f,f)')
