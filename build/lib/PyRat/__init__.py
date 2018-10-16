#!/usr/bin/env python
from PyRat.PyRatBox import *
from PyRat.PyRatEllipsoid import *
from PyRat.PyRatPlane import *
from PyRat.PyRatClone import *
from PyRat.PyRatFacet import *
from PyRat.PyRatRay import *
from PyRat.PyRatCylinder import *
from PyRat.PyRatSpheroid import *
from PyRat.PyRatDisk import *
from PyRat.PyRatObjParser import *
from PyRat.PyRatVisualise import *
from PyRat.PyRatGL import *

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

  print ('ok')


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

