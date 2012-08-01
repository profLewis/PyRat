from PyRatPlane import *
from PyRatDisk import PyRatDisk
from PyRatBox import PyRatTol
try:
  import pp
except:
  pass
import numpy as np

class PyRatCylinder(PyRatPlane):
  '''
  PyRatCylinder: A PyRat object

  P. Lewis 30/6/2012

  For a cylinder
  '''
  def __init__(self,base,tip,contents=None,material=None,info=None):
    '''
    Load the object:

      base : vector of base of cylinder
      tip : vector of top of cylinder

      REQUIRED:
      info['radius']: radius

      OPTIONAL:
      info['caps']: True to put disks on the ends of the cylinder

      info['coords']:
        3 x 2 vectors that define a coordinate for each vertex

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

      Note the the centre of the object is stored as self.base

    '''
    try:
      self.radius = np.abs(info['radius'])
    except:
      self.error('non defined radius defined for cylinder object')
      self.empty = True
      return self

    try:
      self.doCaps = info['caps']
    except:
      self.doCaps = False

    self.base = base
    self.tip = tip

    self.normal = tip - base
    self.length = np.sqrt(np.dot(self.normal,self.normal))
    if self.radius == 0 or self.length == 0:
      self.error('non defined radius or length for cylinder object')
      self.empty = True
      return self

    self.normal /= self.length 
    self.size = 2. * np.pi * self.radius * self.length
    if self.doCaps:
      self.size += 2. * np.pi * self.radius * self.radius
    if self.size <= 0:
      self.empty = True
      return self

    PyRatPlane.__init__(self,self.base,self.normal,\
                             contents=contents,material=material,info=info)
    self.empty = False
    self.base = base

    self.R2 = self.radius*self.radius
    
    V = np.array([self.base+self.radius,self.base-self.radius,\
                  self.tip+self.radius,self.tip-self.radius])

    self.min = np.min(V,axis=0)
    self.max = np.max(V,axis=0)
    self.extent = self.max - self.min

    if self.doCaps:
      info['radius'] = self.radius
      self.caps = [PyRatDisk(self.base,self.normal,material=material,info=info),\
                   PyRatDisk(self.tip,self.normal,material=material,info=info)]

  def intersect(self,ray,closest=True):
    '''
    Intersect ray with cylinder

    ray  : ray descriptor

    Set ray.tnear to ray length

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    import numpy as np
    if self.empty:
      return False
    if not self.ray_on_infinite_cylinder(ray,closest=closest):
      return False
    hitLengths =  ray.tnear
    isHit      =  ray.tfar
    projection = np.array([1.0,0.0,-1,-1])
    for i,ok in enumerate(isHit[2:]): 
      if ok:
        HP = ray.origin+ray.direction*hitLengths[i+2]
        projection[i+2] = np.dot(HP-self.base,self.normal)/self.length 
        if projection[i+2] < 0 or projection[i+2] > 1:
          isHit[i+2] = False
    valid = np.sort(hitLengths[isHit])
    lv = len(valid)
    if lv == 0:
      return False
    if lv >= 1:
      ray.tnear = valid[0]
    if lv >= 2:
      #import pdb;pdb.set_trace()
      ray.tfar = valid[1]
    else:
      ray.tfar = ray.tnear

    ray.rayLengthThroughObject=ray.tfar-ray.tnear
    return True


  def ray_on_infinite_cylinder(self,ray,closest=True):
    '''
    Ray on infinite cylinder test
    '''
    import numpy as np
    len = np.array([0.]*4)
    ok  = np.array([False]*4)

    # deal with caps
    if self.doCaps:
      for i in xrange(2):
        ok[i] = self.caps[i].intersect(ray)
        if ok[i]:
          len[i] = ray.tnear
    v = ray.direction
    p = ray.origin
    pa = self.base
    va = self.normal
    
    deltaP = p - pa
    vva = np.dot(v,va)
    dpva = np.dot(deltaP,va)
    a = v - vva * va
    b = deltaP - dpva * va
    A = np.dot(a,a)
    B = 2.*np.dot(a,b)
    C = np.dot(b,b) - self.R2

    q = B*B - 4.*A*C

    if q<0.0 and (ok[0] or ok[1]):
      ray.tnear = len
      ray.tfar = ok
      return True
    if q<0.:
      return False

    if A == 0:
      A = PyRatTol
    if q == 0:
      len[2] = -B /(2*A)
      ok[2] = True
    else:
      t = np.sqrt(q)
      len[2] = (-B +t) /(2*A)
      len[3] = (-B -t) /(2*A)
      ok[2] = ok[3] = True

    ww = np.where(len<=0)[0]
    ok[ww] = False
    ray.tnear = len
    ray.tfar = ok
    return ok.any()

def main():
  '''
  A simple test of the cylinder algorithm
 
  A scan over a cylinder is made and an image produced
  tests/PyRatCylinder-near.png and tests/PyRatCylinder-far.png with the distances.

  It should be close to 1 at the nearest point (since the camera is located at z=2) for the near
  intersection and 2.0 for the far.
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  from PyRatCylinder import PyRatCylinder
  import pylab as plt

  # set up a test object: a facet
  base = np.array([-0.5,0.,0.])
  tip = np.array([0.5,0.,1.0])
  info = {'radius':0.5,'caps':True}

  cyl = PyRatCylinder(base,tip,info=info)

  # ray direction
  direction = np.array([0,0,-1])

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,2]).astype(float)
  o = origin.copy()
  ray = PyRatRay(o,direction)

  # dimensions of the image in physical units
  dimensions = [2,2]

  result0 = np.zeros(size)
  result1 = np.zeros(size)
  result2 = np.zeros(size)

  # tuple of all parallel python servers to connect with
  isPP = globals().get('pp',False)
  if isPP:
    ppservers = ()
    if len(sys.argv) > 1:
      ncpus = int(sys.argv[1])
      # Creates jobserver with ncpus workers
      job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
      # Creates jobserver with automatically detected number of workers
      job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"
    proc = cyl

    results = [] 
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1]
        ray.length = PyRatBig
        f = job_server.submit(cyl.intersects,(ray,))
        results.append((ix,iy,f))

    for ix,iy,f in results: 
      val = f()
      if val[0]:
        ray = val[1]
        result1[ix,iy] = ray.tnear
        result2[ix,iy] = ray.tfar
  else:
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1]
        ray.length = PyRatBig
        if cyl.intersects(ray):
          n = 0.
          result1[ix,iy] = ray.tnear
          result2[ix,iy] = ray.tfar

  plt.imshow(results0,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRatCylinder-near.png')
  plt.clf()
  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRatCylinder-near.png')
  plt.clf()
  plt.imshow(result2,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRatCylinder-far.png')


if __name__ == "__main__":
    main()
 
