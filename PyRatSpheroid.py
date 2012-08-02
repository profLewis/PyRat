from  PyRatEllipsoid import *

class PyRatSpheroid(PyRatEllipsoid):
  '''
  PyRatSpheroid: A PyRat object

  P. Lewis 30/6/2012

  For a sphere
  '''
  def __init__(self,base,radius,contents=None,material=None,info=None):
    '''
    Load the object:

      base : vector that define the spheroid centre
      radius : vector (for ellipsoid) or scalar for sphere

      OPTIONAL:
      info['coords']:
        3 x 2 vectors that define a coordinate for each vertex

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

      Note the the centre of the object is stored as self.base

    '''
    base[2] -= radius
    PyRatEllipsoid.__init__(self,base,radius,contents=contents,material=material,info=info)

def main():
  '''
  A simple test of the Spheroid algorithm
 
  A scan over a Spheroid is made and an image produced
  tests/PyRatSpheroid-near.png with the distances.

  It should be 1 in the centre (since the camera is located at z=4)
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  from PyRatBox import test

  base = np.array([0,0,2.])
  radius = 1.0
  info = {'verbose':True}

  name = str(globals()['__file__'].split('.')[0])
  test(base,radius,info=info,type=name)

if __name__ == "__main__":
    main()


