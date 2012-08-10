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
  from PyRatObjParser import PyRatObjParser
  filename = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
  try:
    world = PyRatObjParser.load(filename +'.npz')
  except:
    world = PyRatObjParser(filename,verbose=True)
    world.dump(filename +'.npz')

  print 'ok'


def dump(self,filename):
  '''
  Dump a numpy representation
  doesnt work right now
  '''
  np.savez(filename,self=self)

def load(self,filename):
  '''
  unDump a numpy representation
  '''
  return np.load(filename)




if __name__ == "__main__":
    main()

