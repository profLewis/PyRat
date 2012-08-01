from PyRatPlane import *

class PyRatFacet(PyRatPlane):
  '''
  PyRatPlane: A PyRat object

  P. Lewis 30/6/2012

  For a triangular facet
  '''
  def __init__(self,vertices,normal,contents=None,material=None,info=None):
    '''
    Load the object:

      vertices : 3 vectors that define the facet

      normal: 3 x 3 vectors that define the normals at the vertices
              or None

      info['coords']:
        3 x 2 vectors that define a coordinate for each vertex

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

    '''
    self.vertices = vertices

    du = vertices[1] - vertices[0]
    dv = vertices[2] - vertices[0]
    n = np.cross(du,dv)
    mod_normal = np.dot(n,n)
    # area of the triangle
    self.size = np.sqrt(mod_normal*0.5)

    if mod_normal == 0:
      self.error('zero-sized facet') 
      self.empty = True
      return self
    # normalise
    n /= np.sqrt(mod_normal)

    # find max dimension of normal -> facet orientation
    nn = np.abs(n)
    self.orientation = np.where(nn == np.max(nn))[0][0] - 1
    r = np.arange(3)
    # fbase is of dimension (2,3)
    tmp = 0.0
    count = 0
    if tmp == 0:
      self.orientation = (self.orientation+1)%3
      count+=1
      if count > 4:
        self.error('funny-sized facet')
        self.empty = True
        return self
      self.fbase = vertices[0][r!=self.orientation]
      self.Du = du[r!=self.orientation]
      self.Dv = dv[r!=self.orientation]
      tmp = self.Du[0]*self.Dv[1] - self.Dv[0]*self.Du[1]
    self.store = np.array([self.Dv[1],-self.Dv[0],-self.Du[1],self.Du[0]])/tmp
    self.scale = 1./tmp
      
    PyRatPlane.__init__(self,vertices[0],n,\
                             contents=contents,material=material,info=info)
    self.min = np.min(vertices,axis=0)
    self.max = np.max(vertices,axis=0)
    self.extent = self.max - self.min
    self.base = self.min

  def copy(self):
    '''
    Copy the object and combine the
    contents into a flat list.

    It returns a new instance of this class.
    '''
    return  self.__class__(self.vertices,self.normal,\
               contents=self.__tryCopy__(self.contents),\
               material=self.__tryCopy__(self.material),info=self.info)

  def rayToPlane(self,ray,closest=True):
    '''
    Use the PyRatPlane intersect() method
    '''
    if self.empty:
      return False
    return PyRatPlane.intersect(self,ray,closest=closest)

  def intersect(self,ray,closest=True):
    '''
    Intersect ray with facet

    ray  : ray descriptor

    Set ray.tnear to ray length

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length

    '''
    if not self.rayToPlane(ray):
      return False
    if closest and ray.tnear >= ray.length:
      return False
    ray.hitPoint = self.hit(ray,ok=True)
    r = np.arange(3)
    v = ray.hitPoint[r!=self.orientation] - self.fbase
    # local coords
    uvY = (self.store[0]*v[0] + self.store[1]*v[1])*self.scale
    if uvY <0 or uvY > 1:
      return False
    uvX = (self.store[2]*v[0] + self.store[3]*v[1])*self.scale
    if uvX<0 or uvX>1 or uvX+uvY>1+PyRatUnityTol :
      return False
    ray.uv = np.array([uvX,uvY])
    ray.length = ray.tnear

    return True

def main():
  '''
  A simple test of the Facet algorithm
 
  A scan over a trangle is made and an image produced
  tests/PyRatFacet-near.png with the distances.

  It should be 2 in the centre (since the camera is located at z=2)
  and 1 at the lower edge
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a facet
  vertices = np.array([[0,0,0.],[0,1,0],[1,0,1]])

  facet = PyRatFacet(vertices,None)

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
      if facet.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatFacet-near.png')
  plt.clf()


if __name__ == "__main__":
    main()
 
