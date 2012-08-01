from PyRatPlane import *

class PyRatEllipsoid(PyRatPlane):
  '''
  PyRatEllipsoid: A PyRat object

  P. Lewis 30/6/2012

  For an ellipsoid or sphere
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
    self.radius = np.abs(np.atleast_1d(radius))
    if len(self.radius) == 1:
      self.radius = np.abs(np.array([radius,radius,radius]))
 
    if len(self.radius) != 3 or (self.radius <=0).any():
      self.error('non positive radius defined for ellipsoid object')
      self.empty = True
      return self

    self.centre = base + np.array([0,0,self.radius[2]])

    # ellipsoid surface area 
    # approximation http://en.wikipedia.org/wiki/Ellipsoid 
    p = 1.6075
    rp = self.radius ** p
    self.size = np.pi * 4 * ((rp[0]*rp[1]+rp[0]*rp[2]+rp[1]*rp[2])/3.) ** (1./p)

    PyRatPlane.__init__(self,self.centre,None,\
                             contents=contents,material=material,info=info)
    self.empty = False
    self.base = self.centre

    V = np.array([self.centre+self.radius,self.centre-self.radius])

    self.min = np.min(V,axis=0)
    self.max = np.max(V,axis=0)
    self.extent = self.max - self.min

  def intersect(self,ray,closest=True):
    '''
    Intersect ray with ellipsoid

    glassner p.36

    ray  : ray descriptor

    Set ray.tnear to ray length

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    if self.empty:
      return False
    A = ray.origin - self.base
    B = A/self.radius
    C = ray.direction/self.radius
    p = np.dot(B,B)-1.
    q = np.dot(B,C)*2.
    r = np.dot(C,C)
    b = q*q-4*p*r
    if b < 0:
      return False
    b = np.sqrt(b)
    p1 = (-q+b)/(2*r)
    p2 = (-q-b)/(2*r)
    ray.tnear = ray.tfar = 0.
    if p1<p2 and p1 > 0:
      ray.tnear = p1
      ray.tfar = p2
    else:
      ray.tnear = p2
      ray.tfar = p1
    if ray.tnear < 0 or (closest and ray.tnear >= ray.length):
      return False
    ray.rayLengthThroughObject = np.abs(np.max([0,p1])-np.max([0,p2]))
    ray.tfar = ray.tnear + ray.rayLengthThroughObject
    return True

def main():
  '''
  A simple test of the Ellisoid algorithm
 
  A scan over an ellipsoid is made and an image produced
  tests/PyRatEllipsoid-near.png with the distances.

  It should be 1 in the centre (since the camera is located at z=2) for the near
  intersection and 2.0 for the far.
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a facet
  base = np.array([0,0,0.])
  radius = np.array([0.5,1,0.5])

  ell = PyRatEllipsoid(base,radius)

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
      if ell.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))
        distance = ray.tfar * ray.direction
        result2[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatEllipsoid-near.png')
  plt.clf()
  plt.imshow(result2,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRatEllipsoid-far.png')


if __name__ == "__main__":
    main()
 
