#!/usr/bin/env python

# To install:
#
# >>>  pip install -e .
#
# while inside the root directory (Commissioning-tools)

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from setuptools import setup, find_packages

setup(name='jwst_fgs_commissioning_tools',
      version='0.0',
      description='Tools for the fine guidance sensor to use during commissioning of JWST.',
      # long_description='Really, the funniest around.',
      classifiers=[
        # 'Development Status :: 3 - Alpha',
        # 'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Astronomy'
      ],
      keywords='jwst fgs',
      url='https://grit.stsci.edu/wfsc/tools',
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