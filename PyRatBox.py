"""
A PyRat object has methods that allow for:
    - loading (e.g. from a txt string)
    - intersection tests

  This prototype object is defined for an axis-aligned box.

  The core properties are:
    origin : vector of the box origin
    extent : vector of extent

  From which we derive:
    min,max : bound extents

        box aligned to axes

"""

import numpy as np
import sys

# a set of globals 
PyRatTol = 1e-20
PyRatMinExtent = 1e-20
PyRatBig = 1e20
PyRatRayTol = 1e-20
PyRatUnityTol = 0.001



class PyRatBox(object):
  '''
  PyRatBox: A PyRat object 

  P. Lewis 30/6/2012

  A PyRat object has methods that allow for:
    - loading (e.g. from a txt string)
    - intersection tests

  This prototype object is defined for an axis-aligned box.

  The core properties are:
    origin : vector of the box origin
    extent : vector of extent

  From which we derive:
    min,max : bound extents
        
        box aligned to axes
  '''

  def __init__(self,base,extent,contents=None,material=None,info=None):
    '''
    Load the object:

      base  : vector of the minimum of the box
      extent: vector of the box extent

      set min and max, accounting for the fact that
      extent might be set negative

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

    '''
    # load the core descriptors
    self.info = info or {}

    self.base = np.array(base).astype(float)
    if np.array(extent).size == 3:
      self.extent = extent
    else:
      self.extent = np.array(np.zeros(3)+PyRatMinExtent).astype(float)
    max = self.base + self.extent
    self.min = np.min([self.base,max],axis=0)
    self.max = np.max([self.base,max],axis=0)
    self.extent = self.getExtent()
    self.base = self.min

    if 'transparent' in self.info:
      self.invisible = True
    else:
      self.invisible = False

    # check if object is volumetric
    if 'lad' in self.info:
      self.lad = self.info['lad']
      if 'g' in self.info:
        self.g = g
      else:
        self.g = 0.5 #self.spherical

    # some information about what we contain
    self.contents = contents or []
    self.material = material or {}

    # extent can be None
    self.empty = False
    self.size = 0.
    try:
      if not self.extent:
        self.empty = True
        self.extent = np.zeros(3)+PyRatMinExtent
    except:
      if (self.extent == PyRatMinExtent).all():
        self.empty = True
        self.extent = np.zeros(3)+PyRatMinExtent
    if not self.empty:
      self.size = np.prod(self.extent)

  def error(self,msg):
    '''
    Error reporting
    '''
    sys.stderr.write(msg+'\n')

  def copy(self):
    '''
    Copy the object and combine the 
    contents into a flat list.

    It returns a new instance of this class.
    '''
    return self.__class__(self.min.copy(),self.getExtent(),\
               contents=self.__tryCopy__(self.contents),\
               material=self.__tryCopy__(self.material),info=self.info)

  def __tryCopy__(self,x):
    '''
    Try a copy or return the original
    '''
    try:
      return x.copy()
    except:
      return x

  def __iadd__(self,other):
    '''
    The += operator for this class

    Sets this object to one that encompases the limits of self
    and other and contains the combined contents.

    The material dictionary is updated by other.
    '''
    self.min = np.min([self.min,other.min],axis=0)
    self.max = np.max([self.max,other.max],axis=0)
    self.extent = self.getExtent()
    self.base = self.min
    self.empty = not( not self.empty or not other.empty)

    test = [self.contents.append(i) for i in other.contents]
    # update the material
    self.material.update(other.material)

    return self

  def __add__(self,other):
    '''
    The + operator for this class

    Return a *new* instance of this class that encompases the limits of self
    and other and contains the combined contents.

    The material dictionary is updated by other.
    '''
    min = np.min([self.min,other.min],axis=0)
    max = np.max([self.max,other.max],axis=0) 
    contents = self.__tryCopy__(self.contents)
    try:
      test = [contents.append(i) for i in other.contents]
    except:
      test = contents or other.contents

    # update the material
    material = self.__tryCopy__(self.material)
    try:
      material.update(other.material)
    except:
      material = material or other.material

    return self.__class__(min,self.getExtent(max=max,min=min),\
                          contents=contents,material=material,info=self.info)

  def overlap(self,other,all=False):
    '''
    Test to see if other is within bounds of self

    In the case of a Box object, the overlap box is returned

    Arguments:
      other   : the other object to produce the overlap between self and other

    Options:
      all     : when testing contents, we need to decide if an object that is partially
                within scope is allowed. If all is set to False, then if any part of
                the object bounding box is within scope, we include it.
                Since all=False allows ragged edges, we finally must update
                the box extent to admit this.
    '''
    newMin = np.max([self.min,other.min],axis=0)
    newMax = np.min([self.max,other.max],axis=0)

    # combine self and other
    combiBox = self.__add__(other) 

    # now we need to look through the contents and only
    # pass through the contents that lie in the new box
    combiBox.clipContents(newMin,newMaxall=all)    
    return combiBox

  def clipContents(self,min,max,all=False):
    '''
    Clip the contents of self to lie in the extent min/max

    Arguments:
      min   : minimum vector
      max   : max vector

    Options:
      all     : when testing contents, we need to decide if an object that is partially
                within scope is allowed. If all is set to False, then if any part of 
                the object bounding box is within scope, we include it.
                Since all=False allows ragged edges, we finally must update
                the box extent to admit this.
    '''
    contents = []
    for c in self.contents:
      if not c.empty:
        # test the minimum point
        delta = (c.min - self.min)/self.extent
        if ((delta <= 1) & (delta >= 0)).all():
          delta = not all or (c.max - self.min)/self.extent
          if delta != True and ((delta <= 1) & (delta >= 0)).all():
            delta = True
        if delta:
          # we want to append these contents according to the bbox
          if type(c) == Box:
            c.clipContents(min,max,all=all)
          else:
            contents.append(c)
    # delete None objects
    contents = list(np.array(contents)[np.where(contents)])
    for c in contents:
      self.min = np.min([self.min,c.min],axis=0)
      self.max = np.max([self.max,c.max],axis=0)
    self.base = self.min
    self.extent = self.getExtent()  
    self.contents = contents

  def getExtent(self,min=None,max=None):
    extent = np.array(max or self.max).astype(float) - (min or self.min)
    if ((extent == 0) == True).all():
      # there is nothing in here
      return None
    extent[(extent == 0)] = PyRatMinExtent
    return extent

  def intersectProjection(self,d,ray,axis):
    '''
    Object interestion test in projection (axis)

    '''
    r = np.arange(3)
    for i in r[r != axis]:
      intersection = ray.origin[i] + d * ray.direction[i]
      if intersection >= self.min[i] and intersection <= self.max[i]:
        return True

  def bbox_shuffle(self,min,max,origin,direction,ray,axis):
    '''
    BBox intercept code using inequality equations to solve for  
    ray length limits from ray_origin in direction ray_direction 
    when intesecting bbox (return(0) if fail)                    

    min       : min distannce
    max       : max distance
    direction : projection scalar
    origin    : origin scalar

    ray       : full ray description
    axis      : ray axis (i.e. 0,1,2)
    '''
    if direction >= -PyRatTol and direction <= PyRatTol:
      if origin<min or origin>max:
        return False
      return True
    t1=(min-origin)/direction
    t2=(max-origin)/direction
    if t1>t2:
      if t2<ray.tnear and self.intersectProjection(t2,ray,axis):
        ray.tnear=t2
      if t1<ray.tfar and self.intersectProjection(t1,ray,axis):
        ray.tfar=t1
    else:
      if t1<ray.tnear: ray.tnear=t1
      if t2<ray.tfar:  ray.tfar=t2
    if ray.tnear>ray.tfar: return False
    if ray.tfar<0.0:    return False
    return True

  def vintersect(self,ray,closest=True):
    '''
    Interesection test for volumetric object
    '''
    if not self.intersect(ray,closest=closest):
      return False
    if not 'lad' in self.info:
      return True 
    if ray.tnear >= PyRatBig or ray.tfar >= PyRatBig:
      return False
    ray.rayLengthThroughObject=ray.tfar-ray.tnear
    travel = -np.log(np.random.random())/(self.g*self.lad)
   
    if travel >= ray.rayLengthThroughObject:
      return False
    ray.tnear += travel
    ray.tfar = ray.tnear
    ray.object = self
    return True

  def intersect(self,ray,closest=True):
    '''
    Interesection test for ray for box object

    Sets: tnear and tfar in ray object

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    if self.empty or self.invisible:
      return False
    ray.tnear = PyRatBig
    ray.tfar = PyRatBig*2
    ok = [False,False,False]

    for i in xrange(3):
      ok[i] = self.bbox_shuffle(self.min[i],self.max[i],ray.origin[i],ray.direction[i],ray,i)
    if (np.array(ok)==True).all():
      if ray.tnear < 0:
        ray.tnear = 0.
      if closest and ray.tnear >= ray.length:
        return False
      ray.length = ray.tnear
      return True
    return False      

  def intersects(self,ray,closest=True):
    '''
    Call intersect but return ray as well
    '''
    hit = False
    for i in self.contents:
      if not i.invisible:
        thisHit,thisRay = i.intersects(ray.copy(),closest=True)
        if thisHit and thisRay.tnear < ray.length:
          #import pdb;pdb.set_trace()
          ray = thisRay
          ray.length = thisRay.tnear
          hit = thisHit
    thisRay = ray.copy()
    thisHit = self.vintersect(thisRay,closest=True)
    if thisHit or hit:
      ray = thisRay
      ray.length = thisRay.tnear
    return thisHit or hit,ray.copy()

  def surfaceNormal(self,ray,length,true=True):
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
    hitPoint = (self.hit(ray,ok=True) - \
                 self.min)/self.extent
    ihitPoint = (hitPoint+0.5).astype(int)
    delta = hitPoint - ihitPoint
    axis = np.where(delta == np.min(np.abs(delta)))[0]
    upDown = 2*ihitPoint[axis]-1
    direction = np.zeros(3)
    direction[axis] = upDown
    return direction

  def hit(self,ray,ok=False):
    '''
    Calculate intersection point of ray with object
    
    ray : the input ray

    Options:
      ok   : set to True is ray.tnear is already calculated
             otherwise self.intersect(ray) is called
    '''
    ok = ok or self.vintersect(ray)
    if ok: return ray.origin + (ray.length-PyRatRayTol) * ray.direction
    return None

  def spherical(self,angles=None):
    '''
    spherical leaf angle distribution G function
    '''
    return 0.5

  def sortContent(self):
    '''
    This method examines the contents of the PyRatBox and sorts
    them to update plane sets (for faster stacking in ray tracing).
  
    It also updates the box limits (base,extent,min,max)
 
    In effect, calling this method converts the representation from 
    a box to a bounding box
    '''
    pass

def test(base,tip,obj=None,type=None,info={}):
  '''
  A simple test of the intersection algorithm
 
  A scan over an object is made and images produced
  in tests with the distances.

  '''

  from PyRatRay import PyRatRay
  import pylab as plt

  try:
    import pp
    isPP = True
  except:
    isPP = False

  import sys
  import os
  import pylab as plt

  type = type or str(globals()['__file__'].split('.')[0])
  if 'verbose' in info:
    sys.stderr.write('Object: %s\n'%type)

  exec('from %s import %s'%(type,type))

  name = type[5:]
  obj = obj or eval('%s(base,tip,info=info)'%type)
  # ray direction
  direction = np.array([0,0,-1])

  # sun direction
  sun = np.array([1,1,1.])
  sun /= np.sqrt(np.dot(sun,sun))

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,4]).astype(float)
  o = origin.copy()
  ray = PyRatRay(o,direction)

  # dimensions of the image in physical units
  dimensions = [2,2]

  result0 = np.zeros(size)
  result1 = np.zeros(size)
  result2 = np.zeros(size)
  sys.stderr.write('from %s in direction %s\n'%(str(origin),str(direction)))

  if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
  else:
    ncpus = -1
 
  # tuple of all parallel python servers to connect with
  if ncpus != 0 and isPP:
    ppservers = ()
    if ncpus > 0:
      # Creates jobserver with ncpus workers
      job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
      # Creates jobserver with automatically detected number of workers
      job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"

    results = []
    l = float(size[0])
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * 2.*(ix-size[0]*0.5)/size[0]
      if 'verbose' in info and int(100.* (ix+1.)/l) % 5 == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100.*(ix+1.)/l))
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * 2.*(iy-size[1]*0.5)/size[1]
        ray.length = PyRatBig
        f = job_server.submit(obj.intersects,(ray,))
        results.append((ix,iy,f))
    if 'verbose' in info:
      sys.stderr.write('\nGathering results\n')
    l = size[0]*size[1]
    for i,(ix,iy,f) in enumerate(results):
      if 'verbose' in info and int(100.*(i+1.)/l) % 5 == 0:
        sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100.*(i+1.)/l))
      val = f()
      try:
        if val[0]:
          ray = val[1]
          n = np.array([0,0,1.])
          result0[ix,iy] = np.dot(n,-ray.direction)
          result1[ix,iy] = ray.tnear
          result2[ix,iy] = ray.tfar
      except:
        pass
  else:
    l = size[0]
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * 2.*(ix-size[0]*0.5)/size[0]
      if 'verbose' in info and int(100.*(ix+1)/l) % 5 == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*(ix+1)/l))
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * 2.*(iy-size[1]*0.5)/size[1]
        ray.length = ray.tnear = ray.tfar =PyRatBig
        ray.origin = o
        #import pdb;pdb.set_trace()
        hit,thisRay = obj.intersects(ray)
        if hit:
          ray = thisRay
          n = np.array([0,0,1.])
          result0[ix,iy] = np.dot(n,-ray.direction)
          result1[ix,iy] = ray.tnear
          result2[ix,iy] = ray.tfar
  if 'verbose' in info:
    sys.stderr.write('\nWriting results\n')
  plt.clf()
  plt.imshow(result0,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRat%s.png'%name)
  plt.clf()
  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRat%s-near.png'%name)
  plt.clf()
  plt.imshow(result2,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRat%s-far.png'%name)


def main():
  '''
  A simple test of the bounding box algorithm
  
  A scan over a cube is made and 2 images produced
  tests/PyRatBox-near.png and tests/PyRatBox-far.png
  with the near and far distances 
  '''
  min = [-0.5,-0.5,2]
  extent = [1,2,1]
  info = {'verbose':True,'lad':3.0}
  test(min,extent,info=info)


if __name__ == "__main__":
    main()

