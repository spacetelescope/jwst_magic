import os

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import from_levels_and_colors
from matplotlib import cm
import numpy as np
import pandas as pd

from jwst_magic import convert_image


# DO NOT CHANGE THESE LOCATIONS
SEGMENTS = ["A1-1", "A2-2", "A3-3", "A4-4", "A5-5", "A6-6",
            "B1-7","B2-9", "B3-11", "B4-13", "B5-15", "B6-17",
            "C1-8", "C2-10", "C3-12", "C4-14", "C5-16", "C6-18"]

GA_PSF_LOCATIONS = {'A1-1':[1015, 1210], 'A2-2':[856, 1116],   'A3-3':[864, 931],
                    'A4-4':[1023, 837],  'A5-5':[1180, 936],   'A6-6':[1185 , 1126],
                    'B1-7':[1011, 1402], 'B2-9':[691, 1208],   'B3-11':[686, 822],
                    'B4-13':[1027, 676], 'B5-15':[1346, 843],  'B6-17':[1335, 1213],
                    'C1-8':[854, 1300],  'C2-10':[695, 1019],  'C3-12':[866, 742],
                    'C4-14':[1188, 744], 'C5-16':[1342, 1033], 'C6-18':[1172, 1304]
                    }


def distance_to_boresight(x, y, truth_x, truth_y):
    '''
    calculate the distance from the PSF in OTE-01 to it's location in OTE-06/LOS-02
    truth_x and truth_y are the location of the boresight
    '''

    return np.sqrt( (truth_x - x)**2 + (truth_y - y)**2 )


def get_position_from_magic(image, smoothing='default'):
    """ Get x and y postitions of each identifiable PSF in the MOSAIC. The identified PSFs might not all be
    from the target star and some from the target star might be missing"""
    xs, ys, _, _ = convert_image.convert_image_to_raw_fgs.create_all_found_psfs_file(image,
                                                                                     guider=1,
                                                                                     root=None,
                                                                                     out_dir=None,
                                                                                     smoothing=smoothing,
                                                                                     save=False)
    return xs, ys


def read_shadow_log(filename, params_to_include=['x', 'y', 'fwhm', 'fwhm_x', 'fwhm_y', 'ellipse']):
    """
    From the txt file provided from Shadow, pull out helpful information.

    The default is to pull out the x & y locations of each identified PSF, the average FWHM, the FWHM in x
    and y,
    and the ellipticity.
    """
    labels = ['ra_txt', 'dec_txt', 'equinox', 'x', 'y', 'fwhm', 'fwhm_x', 'fwhm_y', 'starsize',
              'ellipse', 'background', 'skylevel', 'brightness', 'time_local', 'time_ut', 'ra_deg',
              'dec_deg']

    info = pd.read_csv(filename, delimiter=' ', skiprows=[0, 1], names=labels)

    info_squeezed = info[params_to_include].copy()
    return info_squeezed


def match_psf_params_to_segment(info_df, matching_dictionary, target_location):
    """
    Match the location and other PSF parameters to the segment based on the matching dictionary. Also
    calculate the distance from each segment to the boresight
    """
    truth_x, truth_y = target_location
    seg_location_dictionary = {matching_dictionary[ind] : {'location': (int(info_df['x'].values[ind]),
                                                                        int(info_df['y'].values[ind])),
                                                           'fwhm': info_df['fwhm'].values[ind],
                                                           'fwhm_x': info_df['fwhm_x'].values[ind],
                                                           'fwhm_y': info_df['fwhm_y'].values[ind],
                                                           'ellipticity': info_df['ellipse'].values[ind],
                                                           'distance_to_target': distance_to_boresight(info_df['x'].values[ind],
                                                                                 info_df['y'].values[ind],
                                                                                 truth_x,
                                                                                 truth_y)
                                                       } for ind in matching_dictionary.keys()}
    return seg_location_dictionary


def get_values(seg_location_dictionary, parameter, decimals=2):
    """
    Get the PSF parameters from the dictionary
    """
    values = []
    for key in seg_location_dictionary:
        values.append(np.round(seg_location_dictionary[key][parameter], decimals))
    return values


def grab_vmin_vmax(param):
    """ Give a vmin & vmax depending on what parameter we are plotting
    """
    if 'fwhm' in param:
        vmin = 10 # TBD
        vmax = 20 # TBD
        cmap = cm.RdYlGn_r
    elif param == 'ellipse':
        vmin = 0
        vmax = 1
        cmap = cm.RdYlGn
    else:
        vmin = None
        vmax = None
        cmap = cm.RdYlGn_r

    return vmin, vmax, cmap


def plot_mosaic_with_psfs(mosaic, car, x_list=None, y_list=None, est_target_location=None, xlim=None, ylim=None):
    """
    Plot out the mosaic image with the identified PSFs and estimated target location plotted
    """
    plt.figure(figsize=(14, 14))

    plt.imshow(mosaic, norm=LogNorm(), origin='lower')
    if x_list is not None and y_list is not None and est_target_location is not None:
        plt.scatter(est_target_location[0], est_target_location[1], s=500, marker='*', c='white',
                    edgecolor='C1', label='Estimated Target Location')
        for j, (x, y) in enumerate(zip(x_list, y_list)):
            plt.annotate(j, (x, y), (x-15, y-20), c='black')
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        plt.legend()

    plt.title(f"{car} Mosaic")
    plt.show()


def plot_each_identified_segment_psf(image, car, x_list, y_list, window_size=40, num_psf_per_row=6):
    """
    Cut each segment out of the mosaic/image provided based on the x and y locations given with
    x_list, and y_list
    """
    radius = window_size//2

    rows = len(x_list)//num_psf_per_row+1
    columns = len(x_list) if len(x_list) < num_psf_per_row else num_psf_per_row
    fig, ax = plt.subplots(rows, columns, figsize=(6*columns, 6*rows))

    i = 0
    for j in range(rows):
        for k in range(columns):
            try:
                cutout = image[int(y_list[i])-radius:int(y_list[i])+radius, int(x_list[i])-radius:int(x_list[i])+radius]
                ax[j, k].imshow(cutout, origin='lower')
                ax[j, k].set_title(f'PSF {i}')
                i+=1
            except IndexError:
                ax[j, k].axis('off')

    plt.show()


def plot_parameters(xs, ys, values, title, cmap=cm.RdYlGn_r):
    """
    Plot the PSF parameters
    """
    center_x = np.median(xs)
    center_y = np.median(ys)

    # plot out new image with distance color coded by how far from boresight
    plt.figure(figsize=(10, 10))
    plt.title(title)
    plt.scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap)
    plt.colorbar()
    plt.ylim(center_y+1024, center_y-1024)
    plt.xlim(center_x-1024, center_x+1024)
    plt.show()


def plot_multiple_parameters(parameters, info_dictionary, xs=None, ys=None):

    # Check that all parameters are in info_dictionary
    parameters = set(parameters) & set(info_dictionary.keys())
    if not parameters:
        print('Nothing to plot')
    else:
        fig, ax = plt.subplots(1, len(parameters), figsize=(6*len(parameters), 8))
        segs = info_dictionary.keys()
        if xs is None and ys is None:
            xs = get_values(info_dictionary, 'x')
            ys = get_values(info_dictionary, 'y')
        for i, (seg, param) in enumerate(zip(segs, parameters)):
            values = get_values(info_dictionary, param)
            for j, (x, y) in enumerate(zip(xs, ys)):
                ax[i].annotate(j, (x, y), (x+15, y+15))
            vmin, vmax, cmap = grab_vmin_vmax(param)
            pl = ax[i].scatter(xs, ys, marker='o', s=100, c=values, cmap=cmap, vmin=vmin, vmax=vmax)
            fig.colorbar(pl, ax=ax[i], orientation='horizontal')
            ax[i].set_title(param)
