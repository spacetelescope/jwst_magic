"""
Script to analyse OTE-01 mosaic and support OTE-02 procedure.

The script's output consists of ideal frame offsets that can be
entered in the observations of the OTE-02 APT program for visual
confirmation.
The returned Ra, Dec will be used to determine the equivalent
pixel coordinates in the individual OTE-01 frame.

Authors
-------
    Johannes Sahlmann


"""
from astropy.table import Table
# from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import pysiaf

print('='*50)
print('OTE-01 mosaic analysis')


########## EXTERNAL INPUTS ##########

# (average) position angle of OTE-01 observations
pa_v3_deg = 0.

# HD 84406
ote01_target = SkyCoord('09h47m30.5542 +63d14m52.11')

# mosaic FITS file
ote01_complete_mosaic = '/Users/jsahlmann/jwst/tel/jwst/rehearsals/team_exercise_january_2019/' \
                        'ote01_analysis/Msc_img.fits'

# pixel positions in ote01_complete_mosaic of the 18 spots
# These have to by X,Y FITS coordinates as determined using ds9, Aladin, or similar
ote01_spot_positions = Table()
ote01_spot_positions['number'] = np.array([1, 8, 18])
ote01_spot_positions['x_pixel'] = [5892, 6440, 7276]  # APT XY FITS coordinates
ote01_spot_positions['y_pixel'] = [10070, 10440, 10637]  # APT XY FITS coordinates

# ote01_spot_positions = Table.read('/Users/jsahlmann/jwst/tel/jwst/rehearsals/team_exercise_january_2019/ote01_analysis/spot_positions.txt', format='ascii.basic', delimiter='\t', names=('number', 'x_pixel', 'y_pixel'), data_start=0, guess=False)

# SIAF aperture defining the ideal frame to work in
aperture_name = 'NRCALL_FULL'
# aperture_name = 'NRCAS_FULL'

########################################

# get SIAF aperture
nircam_siaf = pysiaf.Siaf('NIRCam')
aperture = nircam_siaf[aperture_name]

# get WCS information
mosaic_wcs = WCS(ote01_complete_mosaic)
# mosaic_header = fits.getheader(ote01_complete_mosaic)

# convert spot pixel positions to RA, Dec
sky_coordinates = mosaic_wcs.wcs_pix2world(ote01_spot_positions['x_pixel'], ote01_spot_positions['y_pixel'], 1, ra_dec_order=True)
ote01_spot_positions['ra_deg'] = sky_coordinates[0]
ote01_spot_positions['dec_deg'] = sky_coordinates[1]

# add target to table
ote01_target_pixels = mosaic_wcs.wcs_world2pix(ote01_target.ra.value, ote01_target.dec.value, 1, ra_dec_order=True)
ote01_spot_positions.add_row([0, ote01_target_pixels[0], ote01_target_pixels[1], ote01_target.ra.value, ote01_target.dec.value])

# compute delta_RA, delta_Dec (unused)
spot_catalog = SkyCoord(ra=ote01_spot_positions['ra_deg']*u.deg, dec=ote01_spot_positions['dec_deg']*u.deg)
separations = ote01_target.spherical_offsets_to(spot_catalog)
ote01_spot_positions['delta_rastar_arcsec'], ote01_spot_positions['delta_dec_arcsec'] = separations[0].to(u.arcsec), separations[1].to(u.arcsec)

# approximate attitude computation
attitude = pysiaf.rotations.attitude(aperture.V2Ref, aperture.V3Ref, ote01_target.ra.value, ote01_target.dec.value, pa_v3_deg)

# compute ideal coordinates of spot locations
ote01_spot_positions['delta_x_idl'], ote01_spot_positions['delta_y_idl'] = aperture.tel_to_idl(*pysiaf.rotations.getv2v3(attitude, ote01_spot_positions['ra_deg'], ote01_spot_positions['dec_deg']))

# print to screen
ote01_spot_positions.pprint()

# compute averages
# TODO: implement weighting scheme
spot_index = np.where(ote01_spot_positions['number']!=0)[0]
average_offset_ideal = np.array([np.mean(ote01_spot_positions['delta_x_idl'][spot_index]), np.mean(ote01_spot_positions['delta_y_idl'][spot_index])])

print('='*50)

# for offset special requirement
print('average_correction_ideal (for input as APT offsets) X={0[0]:3.2f} Y={0[1]:3.2f}'.format(-1*average_offset_ideal))

# for next analysis step
average_offset_sky = pysiaf.rotations.pointing(attitude, *aperture.idl_to_tel(*average_offset_ideal))
print('Corresponding RA={0[0]:3.2f} Dec={0[1]:3.2f}'.format(average_offset_sky))
new_pointing = SkyCoord(ra=average_offset_sky[0]*u.deg, dec=average_offset_sky[1]*u.deg)
print('Corresponding {}'.format(new_pointing.to_string(style='hmsdms')))

print('='*50)
