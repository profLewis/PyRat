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

      base : vector that define the ellipsoid base (N.B. not centre)
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
    PyRatEllipsoid.__init__(self,base,radius,contents=contents,material=material,info=info)

def main():
  '''
  A simple test of the Spheroid algorithm

  A scan over an ellipsoid is made and an image produced
  tests/PyRatSpheroid-near.png with the distances.

  It should be 1 in the centre (since the camera is located at z=2) for the near
  intersection and 2.0 for the far.
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a facet
  base = np.array([0,0,0.])
  radius = 0.5

  sph = PyRatSpheroid(base,radius)

  # ray direction
  direction = np.array([0,0,-1])

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,2]).astype(float)
  # dimensions of the image in physical units
  dimensions = [2,2]

  result1 = np.zeros(size)
  result2 = np.zeros(size)

  o = origin.copy()
  ray = PyRatRay(o,direction)
  for ix in xrange(size[0]):
    o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]
    for iy in xrange(size[1]):
      o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1]
      ray.length = PyRatBig
      if sph.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))
        distance = ray.tfar * ray.direction
        result2[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatSpheroid-near.png')
  plt.clf()
  plt.imshow(result2,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRatSpheroid-far.png')


if __name__ == "__main__":
    main()



