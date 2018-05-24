#!/usr/bin/env python
'''
Script to build the jwst_magic package

Authors
-------
    Lauren Chambers

Use
---
    To install in developer mode:
    ::
        pip install -e .

    while inside the root directory (fgs-commissioning)
'''

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from setuptools import setup, find_packages

setup(name='jwst_magic',
      version='0.0',
      description='Multifunctional wavefront Guiding Interface for Commissioning (MaGIC)',
      long_description='Interactive tools to simulate fine guidance sensor data '
                       'and facilitate guiding operations during wavefront '
                       'commissioning of JWST.',
      classifiers=[
        # 'Development Status :: 3 - Alpha',
        # 'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Astronomy'
      ],
      keywords='jwst fgs',
      url='https://grit.stsci.edu/wfsc/tools/fgs-commissioning/',
      author='Keira Brooks, Lauren Chambers, Kathryn St. Laurent',
      # author_email='flyingcircus@example.com',
      # license='MIT',
      packages=find_packages(),  # How will this work with subpackages?
      install_requires=[
          'numpy',
          'astropy',
          'ipython',
          'matplotlib',
          'PyQt5',
          'pysiaf',
          'pyyaml',
          'requests',
          'pytest',
          'photutils'
      ],
      include_package_data=True,
      zip_safe=False)