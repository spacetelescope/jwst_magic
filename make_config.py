''' Make JSON config files
Master file contains the following:

ID (identificaiton):
nx, ny = 2048
nramps = 2
height = 64
overlap = 8
(nstrips = [32,33,34,35,36,37,39,40,42] for overlap = [0,2,4,6,8,10,12,14,16,18])
tcds = 0.338
bias = True
cds image = True
strips image = True

ACQ1 (acquisition 1):
nx, ny = 128
nramps = 6
tcds = 0.3612
bias = True
cds image = True
strips image = False

ACQ2 (acquisition 2):
nx, ny = 32
nramps = 5
tcds = 0.0516
bias = True
cds image = True
strips image = False

TRK (track):
nx, ny = 32
nramps = 5000
tcds = 0.0256
bias = True
cds image = False
strips image = False

LOSTRK (line of sight track):
nx, ny = 43
nramps = 1
tcds = 0.0256
bias = False
cds image = False
strips image = False
'''
import os
import datetime

import pickle
import numpy as np


TCDSID = 0.338    # neil 18 may 12
TCDSACQ1 = 0.3612 # 20 apr 17 keira # 1cds time = 2*(128*128*10.e-6) = 0.32768s
TCDSACQ2 = 0.0516 # 20 apr 17 keira brooks # 2cds time = 4*(32*32*10e-6) = 0.04096s
TCDSTRK = 0.0256  # neil 18 may 12
TCDSFG = 0.0512   # neil 18 may 12

LOCAL_PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(LOCAL_PATH, 'data')

def makeconfig():
    now = datetime.datetime.now()
    filename = "{}{:02d}{:02d}_config.p".format(now.year, now.month, now.day)
    data = build_default_config()
    pickle.dump(data, open(os.path.join(DATA_PATH, filename), 'wb'))

def build_default_config():
    """
    Build the dictionary for the config file
    """
    data = {}
    data['id_dict'] = build_dict(step='ID', imgsize=2048, nramps=2, tcds=TCDSID,
                                 bias=True, cdsimg=True, stripsimg=True)
    data['id_dict'] = update_id_dict(data['id_dict'], overlap=8, height=64)

    data['acq1_dict'] = build_dict(step='ACQ1', imgsize=128, nramps=6,
                                   tcds=TCDSACQ1, bias=True, cdsimg=True,
                                   stripsimg=False)

    data['acq2_dict'] = build_dict(step='ACQ2', imgsize=32, nramps=5,
                                   tcds=TCDSACQ2, bias=True, cdsimg=True,
                                   stripsimg=False)

    data['trk_dict'] = build_dict(step='TRK', imgsize=32, nramps=5000,
                                  tcds=TCDSTRK, bias=True, cdsimg=False,
                                  stripsimg=False)

    data['lostrk_dict'] = build_dict(step='LOSTRK', imgsize=43, nramps=1,
                                     tcds=TCDSTRK, bias=False, cdsimg=False,
                                     stripsimg=False)

    return data


def build_dict(step, imgsize, nramps, tcds, bias=True, cdsimg=True, stripsimg=False):
    '''
    Build the dictionary for any step
    '''
    dictionary = {}

    dictionary['step'] = step
    dictionary['tcds'] = tcds

    dictionary['imgsize'] = imgsize
    dictionary['nramps'] = nramps

    dictionary['bias'] = bias
    dictionary['cdsimg'] = cdsimg
    dictionary['stripsimg'] = stripsimg

    return dictionary

def update_id_dict(id_dict, overlap, height):
    '''
    Add the necessary components to the ID dictionary
    '''
    nstrips_arr = np.array([32, 33, 34, 35, 36, 37, 39, 40, 42])

    id_dict['overlap'] = overlap
    id_dict['height'] = height

    id_dict['nstrips'] = nstrips_arr[int(overlap/2)]

    return id_dict
