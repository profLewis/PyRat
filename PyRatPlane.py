#!/usr/bin/env python
from PyRatBox import *

class PyRatPlane(PyRatBox):
  '''
  PyRatPlane: A PyRat object

  P. Lewis 30/6/2012

  For an infinite plane
  '''
  def __init__(self,base,normal,contents=None,material=None,info=None):
    '''
    Load the object:

      base  : vector of a point on the plane
      normal: normal to the plane (gets normalised)

      min,max etc dont have much meaning in this object

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

    '''
    PyRatBox.__init__(self,base,None,contents=contents,material=material,info=info)
    self.empty = False
    if normal != None:
      self.normal = np.array(normal).astype(float)
      d = dot(self.normal,self.normal)
    else:
      self.normal = np.array([0,0,1.])
      d = 1.
    if d == 0:
      self.empty = True
    else:
      self.normal /= sqrt(d)
    self.dw = -dot(base,self.normal)

  def copy(self):
    '''
    Copy the object and combine the
    contents into a flat list.

    It returns a new instance of this class.
    '''
    return self.__class__(self.min.copy(),self.normal,\
               contents=self.__tryCopy__(self.contents),\
               material=self.__tryCopy__(self.material),info=self.info)

  def surfaceNormal(self,ray,true=True):
    '''
    Return the local surface normal
    where ray intersects object

    all we use in ray is:
      ray.origin
      ray.direction

    Options:
      true  : set True if you want the treu (rather
              than interpolated) surface normal

    '''
    d = dot(ray.direction,self.normal)
    if d > 0:
      ray.localNormal = -self.normal
      return -self.normal
    ray.localNormal = self.normal
    return self.normal

  def intersect(self,ray,closest=True):
    '''
    Intersect ray with infinite plane

    ray  : ray descriptor

    Set ray.tnear to ray length

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    if self.empty:
      return False
    a = np.dot(ray.direction,self.normal)
    if a == 0:
       # ray parallel to plane
      return False
    b = np.dot(ray.origin,self.normal) + self.dw
    if b == 0:
      # ray is on plane
      return False
    D = -b/a
    if D <= 0 or (closest and D >= ray.length):
      return False
    ray.tnear = D
    if ray.tnear < PyRatRayTol or (closest and ray.tnear >= ray.length):
      ray.tnear = ray.tfar = ray.length = PyRatBig
      return False
    return True 

def main():
  '''
  A simple test of the Plane algorithm
 
  A scan over a Plane is made and an image produced
  tests/PyRatPlane-near.png with the distances.

  It should be 1 in the centre (since the camera is located at z=4)
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  from PyRatBox import test
  min = [0.,0.,3]
  normal = [1,1,1]
  info = {'verbose':True}

  name = str(globals()['__file__'].split('/')[-1].split('.')[0])
  test(min,normal,info=info,type=name)

if __name__ == "__main__":
    main()

