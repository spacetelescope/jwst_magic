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

# Local Imports
from . import utils
from .star_selector import select_psfs

# Constants
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
JWST_PUPIL = os.path.join(__location__, 'data', 'JWST_pupil_no_struts.fits')

# Start logger
LOGGER = logging.getLogger(__name__)



class MatchToWss(object):
    def __init__(self, data, global_alignment, match_top=True, match_left=True):
        '''
        Given an image (GA, Image Array, some CMIMF images as of 12/06/2018),
        match each PSF with it's WSS segment number and # IDEA:

        parameters
        ==========
        data: str or array-like
        '''
        # Get data information
        if isinstance(data, str):
            self.data = fits.getdata(data)
            header = fits.getheader(data)
            # Check that data is in raw FGS image frame
            if not (header['INSTRUME'] == 'FGS' or header['INSTRUME'] == 'GUIDER') and header['FILETYPE'] == 'raw':
                raise TypeError("This image is not in the FGS raw frame. Cannot continue.")
        else:
            self.data = data
            LOGGER.warning("Match WSS: If data is not in the FGS raw frame, the " +
                           "matching will NOT be correct.")

        # Define variables
        self.npix_im = np.shape(self.data)[0]
        self.match_top = match_top
        self.match_left = match_left

        self.coords = self.get_coords(global_alignment=global_alignment)
        self.npix_im = np.shape(data)[0]
        self.top_y_im, self.bottom_y_im, self.left_x_im, self.right_x_im, self.outliers = self.define_edges_of_psfs_in_image()


        # Get pupil information
        pupil = fits.getdata(JWST_PUPIL)
        npix_mask = np.shape(pupil)[0]
        ratio = self.define_scaling_factor_for_pupil(pupil)
        # Resize the pupil based on scaling factor
        pupil_scaled = utils.resize_array(pupil, int(np.round(npix_mask*ratio)),
                                          int(np.round(npix_mask*ratio)))

        full_pupil = self.resize_pupil_to_match_im(pupil_scaled)

        # Shift it!
        self.matched_pupil = self.shift_pupil_to_match_im(full_pupil)

        # Grab dictionary and then update it
        wss_segs_dict = MatchToWss.create_wss_seg_dict()
        self.dictionary = self.match_seg_to_psf(wss_segs_dict, self.matched_pupil)


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


    def get_coords(self, global_alignment=False):
        '''
        From start selector get the coordinates of the PSFs associated with the PM
        '''
        _, self.coords, _, _ = select_psfs.manual_star_selection(self.data,
                                                                 global_alignment=global_alignment,
                                                                 testing=True)



    def define_edges_of_psfs_in_image(self):
        """
        Use the image and found coordinates in order to determine the edges (and
        therefore scale and extent) of the PSF array in the image
        """
        # Check for outliers
        avg = np.mean(self.coords, axis=0)
        std = np.std(self.coords, axis=0)
        self.outlier = [c for c in self.coords if (np.array(c) > (avg+std)).all() or (np.array(c) < (avg-std)).all()]
        coords_array = list(set(self.coords)^set(self.outlier))

        if len(self.outlier) > 1:
            LOGGER.warning("WSS Matching: More than two PSFs lie outside the pupil mask. " +
                             "These will be considered missing.")

        #Unzip coordinates
        x, y = zip(*coords_array)

        self.top_y_im = np.asarray(y).max()
        self.bottom_y_im = np.asarray(y).min()
        self.right_x_im = np.asarray(x).max()
        self.left_x_im = np.asarray(x).min()


    def define_edges_of_pupil_mask(self, pupil, top=7, bottom=13, left=16, right=10):
        """
        Use scipy's center of mass function to find the centers of the top-, bottom-,
        left-, and right-most segments. This will give boundaries on the shape of the
        pupil. This makes the assumption that the PSF will approximately line up with
        the center of the segment.

        WARNING: This may not acutally be the case and we need to make sure we test
        many different GA images

        A1 on top of image:
        Top segment = 7
        Bottom segment = 13
        Left segment (middle) = 16
        Right segment (middle) = 10
        """
        # Center_of_mass returns values in y, x
        top_y, _ = ndimage.measurements.center_of_mass(pupil == top)
        bottom_y, _ = ndimage.measurements.center_of_mass(pupil == bottom)

        _, left_x = ndimage.measurements.center_of_mass(pupil == left)
        _, right_x = ndimage.measurements.center_of_mass(pupil == right)

        return top_y, bottom_y, left_x, right_x

    def define_scaling_factor_for_pupil(self, pupil):
        '''
        Find a scaling factor to go between the pupil and the dimentions of the image
        array being fit.
        '''
        top_y_ma, bottom_y_ma, left_x_ma, right_x_ma = self.define_edges_of_pupil_mask(pupil)
        mask_y_dist = int(np.round(top_y_ma - bottom_y_ma))
        mask_x_dist = int(np.round(right_x_ma - left_x_ma))

        im_y_dist = self.top_y_im - self.bottom_y_im
        im_x_dist = self.right_x_im - self.left_x_im

        # Find the ratios in x and y, they may not match, that's okay, I hope
        ratio_y = im_y_dist / mask_y_dist
        ratio_x = im_x_dist / mask_x_dist
        # Find maximum ratio - pupil should be bigger rather than smaller
        ratio = np.max([ratio_y, ratio_x])

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
        shifted_pupil = ndimage.shift(full_pupil, (y_im - y_pupil, x_im - x_pupil),
                                      mode='constant', cval=0.0)

        return shifted_pupil


    def match_seg_to_psf(self, wss_segs_dict, pupil):
        '''
        Use the pupil to pull out the segment number for each PSF.
        '''
        for (x, y) in self.coords:
            try:
                seg = int(pupil[y, x]) # Yes, this should be y, x
                wss_segs_dict[seg]['coords'] = (x, y)
                error = False
            except KeyError:
                error = True

        # If there was an error assigning coords, then either a segment has been
        # kicked out, segments are missing (off the detector or worse), or this
        # image is not going to work for WSS seg matching
        if error:
            for k in wss_segs_dict:
                try:
                    wss_segs_dict[k]['coords']
                except KeyError:
                    if self.outlier:
                        wss_segs_dict[k]['coords'] = (self.outlier[0][0], self.outlier[0][1])
                    else:
                        wss_segs_dict[k]['coords'] = None
                        LOGGER.warning("WSS Matching: Segment %(k)d (%(seg)s) is missing!",
                                       extra={'k':k, 'seg':wss_segs_dict['segid']})

        return wss_segs_dict
