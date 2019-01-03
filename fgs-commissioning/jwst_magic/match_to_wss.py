"""For a given image with an image array configuration, match WSS segment numbers to PSFs

Taking in an image, match the JWST pupil to the PSF configuration and use this
match to determine which segment each PSF is, according to WSS numbering.

Authors
-------
    - Keira Brooks

Use
---
    This module can be imported in a Python shell as such:
    ::
        from jwst_magic import match_to_wss

Notes
-----
    This currently is only set up to work for GA, small image array, large image
    array, and small/large image array with one segment kicked out. In order to
    match WSS segment numbering to more complicated PSF configurations we will
    need to use the OPDs output by WSS of the current mirror configuration and
    use those to extrapolate where we expect the PSF to be in the detector (code
    that will do this is currently being written by L. Chambers for MIRaGe - we
    hope to use this module). One the locations in the detector are extrapolated
    we can create regions where we expect the PSF to lie. This will have to deal
    with half stacked (or fully stacked) PSFs.

    TODO: We also need to be cautious of rotations based on which detector you
    start with. The JWST pupil image should be oriented to the FGS raw frame
    since this is the frame we expect the input image to be in.
"""
# Standard Library Imports
import logging
import os

# Third Party Imports
from astropy.io import fits
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

# Local Imports
from . import utils

# Constants
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
JWST_PUPIL = os.path.join(__location__, 'data', 'JWST_pupil_no_struts_fgs_frame.fits')

# Start logger
LOGGER = logging.getLogger(__name__)



class MatchToWss(object):
    def __init__(self, coords, npix_im=2048, match_top=True, match_left=True,
                 plot=False):
        '''
        Given an image (GA, Image Array, some CMIMF images as of 12/06/2018),
        match each PSF with it's WSS segment number and #.

        This code is meant to be matched with an image in the FGS raw frame.

        parameters
        ==========
        data: str or array-like
        '''

        # Define variables
        self.match_top = match_top
        self.match_left = match_left

        self.coords = coords
        self.npix_im = npix_im
        # Using the coordinates, find the extent of the PSF array
        coords_array = self.check_for_outliers()
        self.define_edges_of_image_array(coords_array)

        # Get pupil information
        self.pupil = fits.getdata(JWST_PUPIL)
        npix_mask = np.shape(self.pupil)[0]
        self.top_seg, self.bottom_seg, self.left_seg, self.right_seg = MatchToWss.determine_edge_segments(self.pupil)
        ratio = self.define_scaling_factor_for_pupil(self.pupil)

        # Resize the pupil based on scaling factor
        self.pupil_scaled = utils.resize_array(self.pupil, int(np.round(npix_mask*ratio)),
                                               int(np.round(npix_mask*ratio)))

        self.full_pupil = self.resize_pupil_to_match_im(self.pupil_scaled)
        # Shift it!
        self.matched_pupil = self.shift_pupil_to_match_im(self.full_pupil)
        self.center_of_array = ndimage.measurements.center_of_mass(self.matched_pupil != 0)

        # Grab dictionary and then update it
        wss_segs_dict = MatchToWss.create_wss_seg_dict()
        self.dictionary = self.match_seg_to_psf(wss_segs_dict, self.matched_pupil)

        if plot:
            self.plot_found_segs(self.matched_pupil)


    @staticmethod
    def create_wss_seg_dict():
        '''
        Create a dictionary of WSS segment numbering that also matches pupil image.
        The coords of the PSF will be added to this dictionary.
        '''
        seg_nums = [1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 17, 8, 10, 12, 14, 16, 18]
        seg_ids = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5',
                   'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']
        wss_segs_dict = {}
        for (seg_num, seg_id) in zip(seg_nums, seg_ids):
            wss_segs_dict[seg_num] = {}
            wss_segs_dict[seg_num]['segid'] = seg_id
        return wss_segs_dict


    def check_for_outliers(self):
        # Check for outliers
        avg = np.mean(self.coords, axis=0)
        std = np.std(self.coords, axis=0)
        self.outlier = [c for c in self.coords if (np.array(c) > (avg+std)).all() or (np.array(c) < (avg-std)).all()]
        coords_array = [c for c in self.coords if (np.array(c) <= (avg+std)).all() or (np.array(c) >= (avg-std)).all()]

        if len(self.outlier) > 1:
            LOGGER.warning("WSS Matching: More than two PSFs lie outside the pupil mask. " +
                           "These will be considered missing.")
        return coords_array

    def define_edges_of_image_array(self, coords_array):
        """
        Use the image and found coordinates in order to determine the edges (and
        therefore scale and extent) of the PSF array in the image
        """
        # Check to see if coords_array is different than the list of all coords
        if coords_array == self.coords:
            coords_array = self.coords
        #Unzip coordinates
        x, y = zip(*coords_array)

        self.top_y_im = np.asarray(y).max()
        self.bottom_y_im = np.asarray(y).min()
        self.right_x_im = np.asarray(x).max()
        self.left_x_im = np.asarray(x).min()

    @staticmethod
    def determine_edge_segments(pupil):
        '''Find the top, bottom, left, right segments'''
        pupil_shape = np.shape(pupil)[0]
        vertical = pupil[:, pupil_shape//2] # left, right
        bottom_seg = int(next((value for value in vertical if value != 0), None))
        top_seg = int(next((value for value in vertical[::-1] if value != 0), None))

        horizontal = pupil[pupil_shape//2, :] # bottom, top
        left_seg = int(next((value for value in horizontal if value != 0), None))
        right_seg = int(next((value for value in horizontal[::-1] if value != 0), None))

        return top_seg, bottom_seg, left_seg, right_seg

    def define_edges_of_pupil_mask(self, pupil):
        """
        Use scipy's center of mass function to find the centers of the top-, bottom-,
        left-, and right-most segments. This will give boundaries on the shape of the
        pupil. This makes the assumption that the PSF will approximately line up with
        the center of the segment based on the rotation of the FGS raw frame.

        THIS IS ONLY VALID FOR THE FOLLOWING PUPIL CONFIGURATION:
                    15 16 17              B5  C5  B6
                  14  5   6  18          C4  A5  A6  C6
        +Y      13  14      1  7       B4  A4      A1  B1
        ^         12  3   2   8          C3  A3  A2  C1
        |           11  10  9              B3  C2  B2
         -> +X

        """
        # Center_of_mass returns values in y, x
        top_y, top_x = ndimage.measurements.center_of_mass(pupil == self.top_seg)
        bottom_y, bottom_x = ndimage.measurements.center_of_mass(pupil == self.bottom_seg)
        if top_x == bottom_x:
            top_pupil, bottom_pupil = top_y, bottom_y
        elif top_y == bottom_y:
            top_pupil, bottom_pupil = top_x, bottom_x

        left_y, left_x = ndimage.measurements.center_of_mass(pupil == self.left_seg)
        right_y, right_x = ndimage.measurements.center_of_mass(pupil == self.right_seg)
        if left_x == right_x:
            left_pupil, right_pupil = left_y, right_y
        elif left_y == right_y:
            left_pupil, right_pupil = left_x, right_x


        return top_pupil, bottom_pupil, left_pupil, right_pupil

    def define_scaling_factor_for_pupil(self, pupil):
        '''
        Find a scaling factor to go between the pupil and the dimentions of the image
        array being fit.
        '''
        top_ma, bottom_ma, left_ma, right_ma = self.define_edges_of_pupil_mask(pupil)
        mask_vert_dist = int(np.round(top_ma - bottom_ma))
        mask_hor_dist = int(np.round(right_ma - left_ma))

        im_vert_dist = self.top_y_im - self.bottom_y_im
        im_hor_dist = self.right_x_im - self.left_x_im

        # Find the ratios in x and y, they may not match, that's okay, I hope
        ratio_vert = im_vert_dist / mask_vert_dist
        ratio_hor = im_hor_dist / mask_hor_dist

        # Find maximum ratio - pupil should be bigger rather than smaller
        ratio = np.max([ratio_vert, ratio_hor])

        return ratio


    def resize_pupil_to_match_im(self, pupil):
        """
        Take the (scaled) pupil and place in array of size of the image
        """
        full_pupil = np.zeros((self.npix_im, self.npix_im))
        npix_mask_new = np.shape(pupil)[0]
        diff = self.npix_im - npix_mask_new

        full_pupil[diff//2:diff//2 + npix_mask_new,
                   diff//2:diff//2 + npix_mask_new] = pupil

        return full_pupil


    def shift_pupil_to_match_im(self, full_pupil):
        """
        Roll the full_pupil to matched top and left edges of image. Make sure to check
        to see if you should be matching the top and left sides of this image with
        the pupil. *This is particularly important when segments are missing or have
        been kicked out of the array.*

        Parameters
        ==========
        full_pupil: array-like
            Array the same size as the image the pupil is to be matched with with the
            JWST pupil in the middle, as defined by resize_pupil_to_match_im.
        y_im: int/float
            The y coordinate of the top- or bottom-most PSF in the PSF array in the
            image, depending on how you are matching your images. (see match_top)
        x_im: int/float
            The x coordinate of the left- or right-most PSF of the PSF array in the
            image, depending on how you are matching your images. (see match_left)
        match_top: bool, optional
            Match the top of the pupil with the top of the PSF array or match the
            bottom (match_top=False). Default: True
        match_left: bool, optional
            Match the left of the pupil with the left of the PSF array or match the
            right side (match_left=False). Default: True
        """
        top_pup, bottom_pup, left_pup, right_pup = self.define_edges_of_pupil_mask(full_pupil)
        if self.match_top:
            y_im = self.top_y_im
            y_pupil = top_pup
        else:
            y_im = self.bottom_y_im
            y_pupil = bottom_pup

        if self.match_left:
            x_im = self.left_x_im
            x_pupil = left_pup
        else:
            x_im = self.right_x_im
            x_pupil = right_pup

        # Shift
        shiftx = np.roll(full_pupil, int((x_im - x_pupil)), axis=0)
        shifted_pupil = np.roll(shiftx, int((y_im - y_pupil)), axis=1)

        return shifted_pupil


    def match_seg_to_psf(self, wss_segs_dict, pupil):
        '''
        Use the pupil to pull out the segment number for each PSF.
        '''
        missing = 0
        for (x, y) in self.coords:
            try:
                seg = int(pupil[int(y), int(x)]) # Yes, this should be y, x
                wss_segs_dict[seg]['coords'] = (int(x), int(y))
            except KeyError:
                missing += 1


        # If there was an error assigning coords, then either a segment has been
        # kicked out, segments are missing (off the detector or worse), or this
        # image is not going to work for WSS seg matching
        if missing == 1:
            for k in wss_segs_dict:
                try:
                    wss_segs_dict[k]['coords']
                except KeyError:
                    if self.outlier:
                        wss_segs_dict[k]['coords'] = (self.outlier[0][0],
                                                      self.outlier[0][1])
                    else:
                        wss_segs_dict[k]['coords'] = None
                        LOGGER.warning("WSS Matching: Segment %(k)d (%(seg)s) is missing!",
                                       extra={'k':k, 'seg':wss_segs_dict['segid']})
        elif missing > 1:
            LOGGER.error("There is more than one missing segment. Cannot"+
                         " match the WSS segment number to the PSFs accurately.")

        return wss_segs_dict

    def plot_found_segs(self, pupil):
        top_final, bottom_final, left_final, right_final = self.define_edges_of_pupil_mask(pupil)

        plt.figure(figsize=(10, 8))
        plt.imshow(pupil, origin='lower')
        plt.axhline(top_final, color='C4', label='pupil edges')
        plt.axhline(bottom_final, color='C4')
        plt.axvline(left_final, color='C4')
        plt.axvline(right_final, color='C4')

        #plt.imshow(self.data, norm=LogNorm(), alpha=0.5, origin='lower')
        plt.axhline(self.top_y_im, color='C1', label='image edges')
        plt.axhline(self.bottom_y_im, color='C1')
        plt.axvline(self.right_x_im, color='C1')
        plt.axvline(self.left_x_im, color='C1')

        for c in self.coords:
            plt.scatter(c[0], c[1], color='C1')

        # plt.ylim(self.bottom_y_im-30, self.top_y_im+30)
        # plt.xlim(self.left_x_im-30, self.right_x_im+30)
        plt.legend()
        plt.show()
