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
    self.normal = np.array(normal).astype(float)
    d = np.dot(self.normal,self.normal)
    if d == 0:
      self.empty = True
    self.normal /= np.sqrt(d)
    self.dw = -np.dot(base,self.normal)

  def copy(self):
    '''
    Copy the object and combine the
    contents into a flat list.

    It returns a new instance of this class.
    '''
    return self.__class__(self.min.copy(),self.normal,\
               contents=self.__tryCopy__(self.contents),\
               material=self.__tryCopy__(self.material),info=self.info)

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
    return True 

def main():
  '''
  A simple test of the Plane algorithm
 
  A scan over a plane with normal [1,1,1] is made and an image produced
  tests/PyRatPlane-near.png with the distances.

  It should be 10 in the centre (since the camera is located at z=10)
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a cube
  min = [0.,0.,0]
  normal = [1,1,1]

  box = PyRatPlane(min,normal)

  # ray direction
  direction = np.array([0,0,-1])

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,10]).astype(float)
  # dimensions of the image in physical units
  dimensions = [4,4]

  result1 = np.zeros(size)

  o = origin.copy()
  ray = PyRatRay(o,direction)
  for ix in xrange(size[0]):
    o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]
    for iy in xrange(size[1]):
      o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1]
      ray.length = PyRatBig
      if box.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatPlane-near.png')
  plt.clf()


if __name__ == "__main__":
    main()
 
