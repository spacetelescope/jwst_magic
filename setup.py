#!/usr/bin/env python
"""
Script to build the jwst_magic package

Authors
-------
    Lauren Chambers
    Shannon Osborne

Use
---
    To install in developer mode:
    ::
        pip install -e .

    while inside the root directory (jwst_magic)
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from setuptools import setup, find_packages
import socket

VERSION = '2.0.0'

INSTALL_REQUIRES = [
    'numpy',
    'astropy',
    'ipython',
    'fgscountrate',
    'matplotlib',
    'notebook',
    'pysiaf>=0.7.1',
    'pyyaml',
    'requests',
    'pytest',
    'photutils',
    'scipy',
    'pytest-qt',
    'pandas',
]

# Determine if PyQt5 needs to be included in the install_requires
try:
    import PyQt5
except ImportError:
    INSTALL_REQUIRES.append('PyQt5')

# Only install pytest-qt if not on SOGS (not available for install)
if "sogs" not in socket.gethostname():
    INSTALL_REQUIRES.append('pytest-qt')

setup(name='jwst_magic',
      version=VERSION,
      description='Multi-Application Guiding Interface for Commissioning (MAGIC)',
      long_description='Interactive tools to simulate fine guidance sensor data '
                       'and facilitate guiding operations during wavefront '
                       'commissioning of JWST.',
      classifiers=[
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy'
      ],
      keywords='jwst fgs',
      url='https://github.com/spacetelescope/jwst_magic',
      author='Keira Brooks, Lauren Chambers, Shannon Osborne, Sherie Holfeltz, Kathryn St. Laurent',
      license='BSD',
      packages=find_packages(),  # How will this work with subpackages?
      install_requires=INSTALL_REQUIRES,
      include_package_data=True,
      zip_safe=False)
