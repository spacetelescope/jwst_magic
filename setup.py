#!/usr/bin/env python
"""
Script to build the jwst_magic package

Authors
-------
    Lauren Chambers

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

INSTALL_REQUIRES = [
    'numpy',
    'astropy',
    'ipython',
    'fgscountrate',
    'matplotlib',
    'notebook',
    'pysiaf',
    'pyyaml',
    'requests',
    'pytest',
    'photutils',
    'scipy'
]

# Determine if PyQt5 needs to be included in the install_requires
try:
    import PyQt5
except ImportError:
    INSTALL_REQUIRES.append('PyQt5')

setup(name='jwst_magic',
      version='0.0',
      description='Multi-Application Guiding Interface for Commissioning (MAGIC)',
      long_description='Interactive tools to simulate fine guidance sensor data '
                       'and facilitate guiding operations during wavefront '
                       'commissioning of JWST.',
      classifiers=[
        # 'Development Status :: 3 - Alpha',
        # 'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy'
      ],
      keywords='jwst fgs',
      url='https://github.com/spacetelescope/jwst_magic',
      author='Keira Brooks, Lauren Chambers, Kathryn St. Laurent',
      # author_email='flyingcircus@example.com',
      # license='MIT',
      packages=find_packages(),  # How will this work with subpackages?
      install_requires=INSTALL_REQUIRES,
      include_package_data=True,
      zip_safe=False)
