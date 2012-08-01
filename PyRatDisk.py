from PyRatPlane import *

class PyRatDisk(PyRatPlane):
  '''
  PyRatDisk: A PyRat object

  P. Lewis 30/6/2012

  For a flat disk
  '''
  def __init__(self,base,normal,contents=None,material=None,info=None):
    '''
    Load the object:

      centre : vectors that define the disk centre

      normal: vector that defines the normal

      REQUIRED:
      info['radius']: disk radius
 
      OPTIONAL:
      info['coords']:
        3 x 2 vectors that define a coordinate for each vertex

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

    '''
    try:
      self.radius = info['radius']
    except:
      self.error('no radius defined for disk object')
      self.empty = True
      return self

    if self.radius <= 0:
      self.error('non positive radius defined fro disk object')
      self.empty = True
      return self

    self.r2 = self.radius *self.radius 
    self.size = np.pi * self.r2

    PyRatPlane.__init__(self,base,normal,\
                             contents=contents,material=material,info=info)
    R = np.ones(3)*self.radius
    
    self.min = np.min([R+base,base-R],axis=0)
    self.max = np.max([R+base,base-R],axis=0)
    self.extent = self.max - self.min
    self.base = base

  def rayToPlane(self,ray):
    '''
    Use the PyRatPlane intersect() method
    '''
    return PyRatPlane.intersect(self,ray)

  def intersect(self,ray,closest=True):
    '''
    Intersect ray with disk

    ray  : ray descriptor

    Set ray.tnear to ray length

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    import numpy as np
    if self.empty:
      return False
    A = self.base - ray.origin
    a = np.dot(A,self.normal)
    b = np.dot(ray.direction,self.normal)
    if b == 0:
      return False
    ray.tnear = a/b
    if ray.tnear < 0 or (closest and ray.tnear >= ray.length):
      return False
    # radial vector
    R = ray.direction*ray.tnear - A
    r2 = np.dot(R,R)
    if r2>self.r2:
      return False
    self.r = np.sqrt(r2)
    return True

  def triangulate(self,N=8):
    '''
    Triangulate the object.

    This method forms a new (traingulated) version of the object
    with N segments.

    It returns a PyRatBox object which contains the set of N PyRatFacet
    objects.

    Options:
      N   : number of segments
    '''
    from PyRatBox import PyRatBox
    from PyRatFacet import PyRatFacet

    base = self.base
    radius = self.radius
    normal = self.normal
    # we need 2 vectors orthogonal to normal
    basis1 = np.array([0.,0.,1.]) 
    basis2 = np.cross(normal,basis1)
    size2 = np.sqrt(np.dot(basis2,basis2))
    if size2 == 0:
      basis1 = np.array([1.,1.,1.])
      basis2 = np.cross(normal,basis1)
      size2 = np.sqrt(np.dot(basis2,basis2))
    basis2 /= size2
    basis1 = np.cross(basis2,normal)
    dTheta = 2. * np.pi/float(N)
    vec = []
    for i in xrange(N):
      Rc=radius*np.cos(i*dTheta)
      Rs=radius*np.cos(i*dTheta)
      vec.append(base + Rc*basis1 + Rs*basis2)
    vec = np.array(vec)
    out = PyRatBox(base,None)
    out.contents = []
    for i in xrange(N):
      v0 = base
      v1 = vec[i-1]
      v2 = vec[i]
      out.contents.append(PyRatFacet(np.array([v0,v1,v2]),normal))
    out.sortContent()
    return out    
     
def main():
  '''
  A simple test of the Facet algorithm
 
  A scan over a disk is made and an image produced
  tests/PyRatDisk-near.png with the distances.

  It should be 2 in the centre (since the camera is located at z=2)
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a facet
  base = np.array([0,0,0.])
  normal = np.array([0,1,1.])
  info = {'radius':1.0}

  disk = PyRatDisk(base,normal,info=info)

  # ray direction
  direction = np.array([0,0,-1])

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,2]).astype(float)
  # dimensions of the image in physical units
  dimensions = [2,2]

  result1 = np.zeros(size)

  o = origin.copy()
  ray = PyRatRay(o,direction)
  for ix in xrange(size[0]):
    o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]
    for iy in xrange(size[1]):
      o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1]
      ray.length = PyRatBig
      if disk.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatDisk-near.png')
  plt.clf()


if __name__ == "__main__":
    main()
 
