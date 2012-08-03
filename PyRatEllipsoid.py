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

    PyRatPlane.__init__(self,self.centre,None,\
                             contents=contents,material=material,info=info)
    self.size = np.pi * 4 * ((rp[0]*rp[1]+rp[0]*rp[2]+rp[1]*rp[2])/3.) ** (1./p)
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
  A simple test of the Ellipsoid algorithm
 
  A scan over a disk is made and an image produced
  tests/PyRatEllipsoid-near.png with the distances.

  It should be 1 in the centre (since the camera is located at z=4)
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a disk
  from PyRatBox import test
  base = np.array([0,0,1])
  radius = np.array([0.5,1,1])
  info = {'verbose':True}

  name = str(globals()['__file__'].split('.')[0])
  test(base,radius,info=info,type=name)

if __name__ == "__main__":
    main()
 
