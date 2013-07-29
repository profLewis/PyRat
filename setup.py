#!/usr/bin/env python

from setuptools import setup, find_packages
import sys, os

DISTNAME = "PyRat"
DESCRIPTION = "Python ray tracing"
LONG_DESCRIPTION = open('README.md').read()
MAINTAINER = 'P Lewis'
MAINTAINER_EMAIL = "p.lewis@ucl.ac.uk"
URL = 'https://github.com/profLewis/PyRat'
LICENSE = 'Undecided'
VERSION = "1.1"
DOWNLOAD_URL="https://github.com/profLewis/PyRat/zipball/master"

setup(name='PyRat',
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      keywords='',
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license='',
      packages=['PyRat'],
      zip_safe=False,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Intended Audience :: Developers',
              'License :: OSI Approved',
              'Programming Language :: Fortran',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'
              ]
)
