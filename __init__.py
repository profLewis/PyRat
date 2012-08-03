#!/usr/bin/env python
from PyRatBox import *
from PyRatEllipsoid import PyRatEllipsoid
from PyRatRay import PyRatRay
from PyRatCylinder import PyRatCylinder
from PyRatFacet import PyRatFacet
from PyRatSpheroid import PyRatSpheroid
from PyRatDisk import PyRatDisk
from PyRatPlane import PyRatPlane
from PyRatObjParser import PyRatObjParser
from PyRatClone import PyRatClone

try:
  import pp
  PyRatIsPP = True
except:
  PyRatIsPP = False

import numpy as np

def main():
  '''
  Test code
  '''
  filename = 'PyRat/spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
  p = PyRatObjParser(filename,reportingFrequency=100,verbose=True)

if __name__ == "__main__":
    main()

