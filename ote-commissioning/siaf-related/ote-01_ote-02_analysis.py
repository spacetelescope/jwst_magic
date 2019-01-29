from astropy.io import fits
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import pysiaf

# ote01 target
# TYC-5309-277-1 ?
# ote01_target = SkyCoord('03h36m43.5650 -14d31m39.45')

print('='*50)
print('OTE-01 mosaic analysis')

########## EXTERNAL INPUTS ##########

pa_v3_deg = 0.

# HD 84406
ote01_target = SkyCoord('09h47m30.5542 +63d14m52.11')

ote01_complete_mosaic = '/Users/jsahlmann/jwst/tel/jwst/rehearsals/team_exercise_january_2019/' \
                        'ote01_analysis/Msc_img.fits'
# ote01_complete_mosaic = '/Users/jsahlmann/jwst/tel/jwst/rehearsals/team_exercise_january_2019/' \
#                         'ote01_analysis/Full_mosaic_Jan2019.fits'

# pixel positions in ote01_complete_mosaic
# ote01_spot_positions = Table.read('/Users/jsahlmann/jwst/tel/jwst/rehearsals/team_exercise_january_2019/ote01_analysis/spot_positions.txt', format='ascii.basic', delimiter='\t', names=('number', 'x_pixel', 'y_pixel'), data_start=0, guess=False)
ote01_spot_positions = Table()
ote01_spot_positions['number'] = np.array([1, 8, 18])

ote01_spot_positions['x_pixel'] = [5892, 6440, 7276]  # APT XY FITS coordinates
ote01_spot_positions['y_pixel'] = [10070, 10440, 10637]  # APT XY FITS coordinates

# aperture_name = 'NRCALL_FULL'
aperture_name = 'NRCAS_FULL'
nircam_siaf = pysiaf.Siaf('NIRCam')
aperture = nircam_siaf[aperture_name]

########################################



mosaic_header = fits.getheader(ote01_complete_mosaic)

mosaic_wcs = WCS(ote01_complete_mosaic)

# convert to RA, Dec
# pixel_positions = np.array(ote01_spot_positions['x_pixel', 'y_pixel'].to_pandas())
# print(mosaic_wcs.wcs_pix2world(pixel_positions), 1)
sky_coordinates = mosaic_wcs.wcs_pix2world(ote01_spot_positions['x_pixel'], ote01_spot_positions['y_pixel'], 1, ra_dec_order=True)
ote01_spot_positions['ra_deg'] = sky_coordinates[0]
ote01_spot_positions['dec_deg'] = sky_coordinates[1]


# add target to table
ote01_target_pixels = mosaic_wcs.wcs_world2pix(ote01_target.ra.value, ote01_target.dec.value, 1, ra_dec_order=True)
ote01_spot_positions.add_row([0, ote01_target_pixels[0], ote01_target_pixels[1], ote01_target.ra.value, ote01_target.dec.value])
# print('Target location in ote01_complete_mosaic pixels: X={} Y={}'.format())

spot_catalog = SkyCoord(ra=ote01_spot_positions['ra_deg']*u.deg, dec=ote01_spot_positions['dec_deg']*u.deg)
# print(`separation_3d(spot_catalog))
separations = ote01_target.spherical_offsets_to(spot_catalog)


ote01_spot_positions['delta_rastar_arcsec'], ote01_spot_positions['delta_dec_arcsec'] = separations[0].to(u.arcsec), separations[1].to(u.arcsec)



# attitude computation
attitude = pysiaf.rotations.attitude(aperture.V2Ref, aperture.V3Ref, ote01_target.ra.value, ote01_target.dec.value, pa_v3_deg)

ote01_spot_positions['delta_x_idl'], ote01_spot_positions['delta_y_idl'] = aperture.tel_to_idl(*pysiaf.rotations.getv2v3(attitude, ote01_spot_positions['ra_deg'], ote01_spot_positions['dec_deg']))
# ote01_spot_positions['V2_arcsec'], ote01_spot_positions['V3_arcsec'] = aperture.tel_to_idl(pysiaf.rotations.getv2v3(attitude, ote01_spot_positions['ra_deg'], ote01_spot_positions['dec_deg']))


ote01_spot_positions.pprint()

spot_index = np.where(ote01_spot_positions['number']!=0)[0]
average_offset_ideal = np.array([np.mean(ote01_spot_positions['delta_x_idl'][spot_index]), np.mean(ote01_spot_positions['delta_y_idl'][spot_index])])
print('='*50)
print('average_correction_ideal (for APT) X={0[0]:3.2f} Y={0[1]:3.2f}'.format(-1*average_offset_ideal))
average_offset_sky = pysiaf.rotations.pointing(attitude, *aperture.idl_to_tel(*average_offset_ideal))
print('Corresponding RA={0[0]:3.2f} Dec={0[1]:3.2f}'.format(average_offset_sky))
new_pointing = SkyCoord(ra=average_offset_sky[0]*u.deg, dec=average_offset_sky[1]*u.deg)
print('Corresponding {}'.format(new_pointing.to_string(style='hmsdms')))


print('='*50)


