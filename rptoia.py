'''
Convert real pixel coordinates to ideal angle
'''
#THIRD PARTY
import numpy as np

#LOCAL
import log

XOFFSET = 1023.5
YOFFSET = 1023.5
PLATESCALE = 0.06738281367
GAIN = 1.5


def rptoia(xarr, yarr, guider):
    '''
    Pass in X Y in real pixels and get out X Y in ideal angle
    '''
    yrealpixel, xrealpixel = np.asarray(xarr), np.asarray(yarr) # THIS IS NECESSARY, DO NOT QUESTION
    coeff_dict = get_conversion_factors(guider)

    #  Reverse calculation of the coordinate conversion
    xidealpixel = convert_to_ideal_pixel(coeff_dict, 'X', xrealpixel, yrealpixel)
    yidealpixel = convert_to_ideal_pixel(coeff_dict, 'Y', xrealpixel, yrealpixel)

    xangle = (xidealpixel - XOFFSET) * PLATESCALE
    if guider == 1:
        yangle = (YOFFSET - yidealpixel) * PLATESCALE
    elif guider == 2:
        yangle = (yidealpixel - YOFFSET)*PLATESCALE

    return xangle, yangle

def convert_to_ideal_pixel(coeff_dict, xy_str, xrealpixel, yrealpixel):
    '''
    Convert x- or y- real pixel to ideal pixel

    Parameters
    ----------
    coeff_dict: dictionary
        Dictionary of all the X and Y coefficients for both guider1 and guider2
    xy_str: str
        "X" or "Y" based on if you are converting the x or y coordinate
    xrealpixel: float
        The x value passed in by the user
    yrealpixel: float
        The y value passed in by the user

    Returns
    -------
    idealpixel: float
        The converted coordinate in ideal pixels
    '''

    idealpixel = coeff_dict[xy_str]['Coeff0']
    idealpixel += coeff_dict[xy_str]['Coeff1'] * xrealpixel
    idealpixel += coeff_dict[xy_str]['Coeff2'] * yrealpixel
    idealpixel += coeff_dict[xy_str]['Coeff3'] * xrealpixel**2
    idealpixel += coeff_dict[xy_str]['Coeff4'] * xrealpixel * yrealpixel
    idealpixel += coeff_dict[xy_str]['Coeff5'] * yrealpixel**2
    idealpixel += coeff_dict[xy_str]['Coeff6'] * xrealpixel**3
    idealpixel += coeff_dict[xy_str]['Coeff7'] * xrealpixel**2 * yrealpixel
    idealpixel += coeff_dict[xy_str]['Coeff8'] * xrealpixel * yrealpixel**2
    idealpixel += coeff_dict[xy_str]['Coeff9'] * yrealpixel**3
    idealpixel += coeff_dict[xy_str]['Coeff10'] * xrealpixel**4
    idealpixel += coeff_dict[xy_str]['Coeff11'] * xrealpixel**3 * yrealpixel
    idealpixel += coeff_dict[xy_str]['Coeff12'] * xrealpixel**2 * yrealpixel**2
    idealpixel += coeff_dict[xy_str]['Coeff13'] * xrealpixel * yrealpixel**3
    idealpixel += coeff_dict[xy_str]['Coeff14'] * yrealpixel**4

    return idealpixel

def get_conversion_factors(guider):
    '''Get the coefficents for the conversion to ideal pixel for guiders 1 and 2
    '''
    coeff_dict = {}
    coeff_dict['X'] = {}
    coeff_dict['Y'] = {}

    if guider == 1:
        coeff_dict['X']['Coeff0'] = -2.3132211E+01
        coeff_dict['X']['Coeff1'] = 9.9858571E-01
        coeff_dict['X']['Coeff2'] = 1.0458177E-02
        coeff_dict['X']['Coeff3'] = 2.6914738E-06
        coeff_dict['X']['Coeff4'] = 6.7167416E-06
        coeff_dict['X']['Coeff5'] = 9.9063452E-07
        coeff_dict['X']['Coeff6'] = 1.2100144E-09
        coeff_dict['X']['Coeff7'] = -3.4359146E-11
        coeff_dict['X']['Coeff8'] = 1.2718376E-09
        coeff_dict['X']['Coeff9'] = -2.0353075E-11
        coeff_dict['X']['Coeff10'] = 5.3149090E-14
        coeff_dict['X']['Coeff11'] = 9.3076979E-14
        coeff_dict['X']['Coeff12'] = 8.3907183E-14
        coeff_dict['X']['Coeff13'] = 9.7074620E-14
        coeff_dict['X']['Coeff14'] = 1.9366816E-14

        coeff_dict['Y']['Coeff0'] = -2.6649388E+01
        coeff_dict['Y']['Coeff1'] = -2.5585494E-03
        coeff_dict['Y']['Coeff2'] = 1.0114078E+00
        coeff_dict['Y']['Coeff3'] = 2.4056198E-06
        coeff_dict['Y']['Coeff4'] = 2.0788818E-06
        coeff_dict['Y']['Coeff5'] = 9.3110110E-06
        coeff_dict['Y']['Coeff6'] = -1.6052215E-11
        coeff_dict['Y']['Coeff7'] = 1.2730773E-09
        coeff_dict['Y']['Coeff8'] = -9.1059163E-11
        coeff_dict['Y']['Coeff9'] = 1.3291468E-09
        coeff_dict['Y']['Coeff10'] = 2.1479653E-14
        coeff_dict['Y']['Coeff11'] = 5.2933520E-14
        coeff_dict['Y']['Coeff12'] = 1.4401378E-13
        coeff_dict['Y']['Coeff13'] = 7.4524608E-14
        coeff_dict['Y']['Coeff14'] = 1.2816500E-13

    elif guider == 2:
        coeff_dict['X']['Coeff0'] = -3.2653125E+01
        coeff_dict['X']['Coeff1'] = 1.0343933E+00
        coeff_dict['X']['Coeff2'] = 2.1128135E-02
        coeff_dict['X']['Coeff3'] = -9.0969494E-06
        coeff_dict['X']['Coeff4'] = -1.4320516E-05
        coeff_dict['X']['Coeff5'] = -3.9403361E-06
        coeff_dict['X']['Coeff6'] = 1.6122843E-09
        coeff_dict['X']['Coeff7'] = 5.9346314E-10
        coeff_dict['X']['Coeff8'] = 2.0380335E-09
        coeff_dict['X']['Coeff9'] = 2.0974142E-10
        coeff_dict['X']['Coeff10'] = -2.8502771E-14
        coeff_dict['X']['Coeff11'] = -8.3569450E-14
        coeff_dict['X']['Coeff12'] = -5.0976659E-14
        coeff_dict['X']['Coeff13'] = -9.1800561E-14
        coeff_dict['X']['Coeff14'] = -8.9563271E-15

        coeff_dict['Y'] = {}
        coeff_dict['Y']['Coeff0'] = -7.7165127E+01
        coeff_dict['Y']['Coeff1'] = 2.9152735E-02
        coeff_dict['Y']['Coeff2'] = 1.0773690E+00
        coeff_dict['Y']['Coeff3'] = -6.2761496E-06
        coeff_dict['Y']['Coeff4'] = -7.8529777E-06
        coeff_dict['Y']['Coeff5'] = -2.1370284E-05
        coeff_dict['Y']['Coeff6'] = 2.0146586E-10
        coeff_dict['Y']['Coeff7'] = 2.0163576E-09
        coeff_dict['Y']['Coeff8'] = 6.6916423E-10
        coeff_dict['Y']['Coeff9'] = 2.4400761E-09
        coeff_dict['Y']['Coeff10'] = -1.9936426E-14
        coeff_dict['Y']['Coeff11'] = -3.2487010E-14
        coeff_dict['Y']['Coeff12'] = -1.3263445E-13
        coeff_dict['Y']['Coeff13'] = -4.2058076E-14
        coeff_dict['Y']['Coeff14'] = -1.2323597E-13

    else:
        log.error("Do not recognize guider {}. Exiting.".format(guider))
        raise StandardError("Guider not recognized.")

    return coeff_dict


def write_to_file(xangle, yangle):
    ''' Write ideal angles to file'''
    with open('ideal.tmp', 'w') as f:
        for i, (xang, yang) in enumerate(zip(xangle, yangle)):
            f.writeto(i-1, xang, yang)
