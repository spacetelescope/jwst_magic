"""
SIAF related tools to support OTE commissioning programs.

Authors
-------
    Johannes Sahlmann

Content
-------
    ote_2_nircam_to_fgs : support OTE-02
        Example output to stdout:
        ==================================================
        Science pixel coordinates X=300 Y=400 in NRCA1_FULL correspond to
        Ideal frame coordinates X_idl=60.50 arcsec Y_idl=151.95 arcsec in FGS1_FULL
        ==================================================

"""
import numpy as np
import pylab as pl

import pysiaf

def ote_2_nircam_to_fgs(nircam_science_pixel_coordinates):
    """Use operational SIAF in pysiaf to transform between science and ideal frames of two apertures.

    Parameters
    ----------
    nircam_science_pixel_coordinates : numpy array with two elements

    """
    # get SIAF
    nircam_siaf = pysiaf.Siaf('nircam')
    fgs_siaf = pysiaf.Siaf('fgs')

    nircam_aperture_name = 'NRCA1_FULL'
    fgs_aperture_name = 'FGS1_FULL'

    # get aperture objects
    nircam_aperture = nircam_siaf[nircam_aperture_name]
    fgs_aperture =fgs_siaf[fgs_aperture_name]

    # transform to V2,V3
    nircam_tel = nircam_aperture.sci_to_tel(nircam_science_pixel_coordinates[0], nircam_science_pixel_coordinates[1])

    # transform to FGS ideal frame
    fgs_idl = fgs_aperture.tel_to_idl(*nircam_tel)

    # print results
    print('='*50)
    print('Science pixel coordinates X={0[0]} Y={0[1]} in {1} correspond to'.format(nircam_science_pixel_coordinates, nircam_aperture_name))
    print('Ideal frame coordinates X_idl={0[0]:2.2f} arcsec Y_idl={0[1]:2.2f} arcsec in {1}'.format(fgs_idl, fgs_aperture_name))
    print('='*50)

    # make figure
    pl.close('all')
    pl.figure()
    nircam_aperture.plot(label=True, mark_ref=True)
    fgs_aperture.plot(label=True, mark_ref=True)
    pl.plot(nircam_tel[0], nircam_tel[1], 'bo')
    pl.plot([fgs_aperture.V2Ref, nircam_tel[0]], [fgs_aperture.V3Ref, nircam_tel[1]], 'k--')
    pl.show()

# run code with NIRCam input coordinates determined as part of OTE-02
nircam_science_pixel_coordinates = np.array([300, 400])
ote_2_nircam_to_fgs(nircam_science_pixel_coordinates)
