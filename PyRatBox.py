#!/usr/bin/env python
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
import math
from math import sqrt

# a set of globals 
PyRatTol = 1e-20
PyRatMinExtent = 1e-20
PyRatBig = 1e20
PyRatRayTol = 1.e-10 # 0.00001
PyRatUnityTol = 0.001

def abs(a):
  if a < 0:
    return -a
  return a

def dot(a,b):
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

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

  def __init__(self,base,extent,planeSets=False,contents=None,material=None,info=None):
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
    self.planeSets = planeSets
    self.planes = []

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

  def report(self,level=0):
    '''
    Report on object contents
    '''
    buff = '_'*level
    self.error('%05d %s{'%(level,buff))
    self.error('%s self     \t%s'%(buff,str(self)))
    self.error('%s type     \t%s'%(buff,str(type(self))))
    try:
      self.error('%s offset   \t%s'%(buff,str(self.offset)))
    except:
      pass
    try:
      self.error('%s matrix  \t%s'%(buff,str(self.matrix)))
    except:
      pass
    try:
      self.error('%s tip      \t%s'%(buff,str(self.tip)))
    except:
      pass
    try:
      self.error('%s centre    \t%s'%(buff,str(self.centre)))
    except:
      pass
    try:
      self.error('%s radius    \t%s'%(buff,str(self.radius)))
    except:
      pass
    try:
      self.error('%s normal    \t%s'%(buff,str(self.normal)))
    except:
      pass

    self.error('%s min      \t%s'%(buff,str(self.min)))
    self.error('%s max      \t%s'%(buff,str(self.max)))
    self.error('%s base     \t%s'%(buff,str(self.base)))
    self.error('%s extent   \t%s'%(buff,str(self.extent)))
    self.error('%s invisible\t%s'%(buff,str(self.invisible)))
    self.error('%s empty    \t%s'%(buff,str(self.empty)))
    self.error('%s material \t%s'%(buff,str(self.material)))
    self.error('%s size     \t%s'%(buff,str(self.size)))
    self.error('%s info     \t%s'%(buff,str(self.info)))
    self.error('%s contents:'%(buff))
    for c in self.contents:
      c.report(level=level+1)
    self.error('%05d %s}'%(level,buff))

  def updateContents(self):
    '''
    Resort the contents
    '''
    allMin = np.array([i.min for i in self.contents])
    allMax = np.array([i.max for i in self.contents])
    self.minP = np.argsort(allMin,axis=0)    
    self.maxP = np.argsort(allMax,axis=0)
    self.minN = np.argsort(-allMin,axis=0)           
    self.maxN = np.argsort(-allMax,axis=0)

  def updateBbox(self):
    import numpy as np
    from PyRatClone import PyRatClone
    # look at min and max
    if type(self) != PyRatClone:
      return
    corners = np.array([[self.min[0],self.min[1],self.min[2]]\
                       ,[self.min[0],self.min[1],self.max[2]]\
                       ,[self.min[0],self.max[1],self.min[2]]\
                       ,[self.min[0],self.max[1],self.max[2]]\
                       ,[self.max[0],self.min[1],self.min[2]]\
                       ,[self.max[0],self.min[1],self.max[2]]\
                       ,[self.max[0],self.max[1],self.min[2]]\
                       ,[self.max[0],self.max[1],self.max[2]]])
    # now apply the matrix and then translation
    for n,i in enumerate(corners):
      try:
        corners[n] = np.array(np.matrix(i)*self.matrix).flatten()
      except:
        pass
      try:
        corners[n] += self.offset
      except:
        pass
    self.min = np.min(corners,axis=0)
    self.max = np.max(corners,axis=0)
    try:
      self.scale = np.linalg.det(self.matrix)**(1./3.)
    except:
      pass


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
    intersection = (ray.origin + d * ray.direction - self.min)/self.extent
    if not ((intersection < 0) | (intersection > 1)).all():
      return True
    return False
    #r = np.arange(3)
    #for i in r[r != axis]:
    #  intersection = ray.origin[i] + d * ray.direction[i]
    #  if intersection >= self.min[i] and intersection <= self.max[i]:
    #    return True

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
        return False
    else:
      if t1<ray.tnear and self.intersectProjection(t1,ray,axis):
        ray.tnear=t1
        if t2<ray.tfar and self.intersectProjection(t2,ray,axis):
          ray.tfar=t2
      else:
        return False
    if ray.tnear>ray.tfar: return False
    if ray.tfar<0.0:    return False
    return True

  def vintersect(self,ray,closest=True):
    '''
    Interesection test for volumetric object
    '''
    from PyRatBox import PyRatBox
    from PyRatClone import PyRatClone
    thisRay = ray.copy()
    if not self.intersect(thisRay,closest=closest) or thisRay.length <= 0:
      return False
    # so we intersect the superstructure
    # if its a box, then look into its contents
    if type(self) == PyRatBox or type(self) == PyRatClone:
      if self.invisible:
        return True
    if thisRay.tnear <= 0:
      return False
    ray.tnear = thisRay.tnear
    ray.tfar = thisRay.tfar
    try:
      ray.object = thisRay.object
    except:
      ray.object = None
    if not 'lad' in self.info:
      ray.length = ray.tnear
      ray.object = self
      try:
        ray.object.surfaceNormal(ray)
      except:
        pass
      return True 
    if ray.tnear >= ray.big or ray.tfar >= ray.big:
      return False
    ray.length = ray.tnear
    ray.rayLengthThroughObject=ray.tfar-ray.tnear
    travel = -np.log(np.random.random())/(self.g*self.lad)
   
    if travel >= ray.rayLengthThroughObject:
      return False
    ray.tnear += travel
    ray.tfar = ray.tnear
    ray.object = self
    ray.localNormal = np.array([0.,0.,1.])
    return True

  def intersectOld(self,ray,closest=True):
    '''
    Interesection test for ray for box object

    Sets: tnear and tfar in ray object

    Options:
      closest : set True to return False if the possible
                ray length would be greater than ray.length
    '''
    #if self.empty or self.invisible:
    #  return False
    ray.tnear = ray.big
    ray.tfar = ray.big
    ok = [False,False,False]

    for i in xrange(3):
      ok[i] = self.bbox_shuffle(self.min[i],self.max[i],ray.origin[i],ray.direction[i],ray,i)
    if (np.array(ok)==True).all():
      if ray.tnear < 0:
        ray.tnear = 0.
      if closest and ray.tnear >= ray.length:
        return False
      ray.length = ray.tnear
      #import pdb;pdb.set_trace()
      return True

    return False  

  def intersect(self,r,closest=True):
    '''
    Intersection test for a bounding box
    '''
    direction = r.direction.copy()
    direction[direction==0] = PyRatTol
    tmin = (self.min - r.origin) / direction
    tmax = (self.max - r.origin) / direction
    tt = tmin.copy()
    tt[tmax<tmin] = tmax[tmax<tmin]
    tmax[tmin>tmax] = tmin[tmin>tmax]
    tmin = tt

    if (tmin[0] > tmax[1]) or (tmin[1] > tmax[0]):
      return False
    if tmin[1] > tmin[0]:
      tmin[0] = tmin[1]
    if tmax[1] < tmax[0]:
      tmax[0] = tmax[1]
    if (tmin[0] > tmax[2]) or (tmin[2] > tmax[0]):
      return False

    if tmin[2] > tmin[0]:
      tmin[0] = tmin[2]
    if tmax[2] < tmax[0]:
      tmax[0] = tmax[2]
    r.tnear = tmin[0]
    r.tfar = tmax[0]
    if closest and r.tnear >= r.length:
      return False
    return tmin[0] < r.length

  def intersectsMany(self,rays,closest=True):
    '''
    Call intersects for multiple rays
    '''
    result = []
    for ray in rays:
      result.append(self.intersects(ray,closest=closest))
    return result

  def rayTreeMany(self,rays,closest=True):
    '''
    Call intersects for multiple rays
    '''
    result = []
    for ray in rays:
      result.append(self.rayTree(ray,closest=closest))
    return result

  def isDefined(self):
    from PyRatBox import PyRatBox
    if type(self) != PyRatBox:
      return False 
    try:
      return self.defined
    except:
      return False

  def intersects(self,ray,closest=True):
    '''
    'solid' object intersection tests

    Call intersect but return ray as well
    '''
    from PyRatClone import PyRatClone
    # cheap method if only one entry
    invisible = self.invisible# or type(self) == PyRatClone
    if invisible and len(self.contents) == 1:
      return self.contents[0].intersects(ray,closest=closest)

    # first check to see if we hit this object
    thisRay = ray.copy()
    thisHit = thatHit = self.vintersect(thisRay,closest=True)
    if not thisHit:
      return False,ray
    if invisible:
      # it doesnt count so forget about it
      thisHit = False
      thisRay = ray
    else:
      thisRay.length = thisRay.tnear
    # if not, then don't bother with contents
    # if we do hit, then have a look at the contents
    ray.ccopy(thisRay)
    hit = False
    if len(self.contents) == 0:
      return thisHit or hit,ray.copy()

    axis = np.argsort(np.array([-ray.direction,ray.direction]).flatten())[0]

    if self.planeSets:
      try:
        x = self.minN[axis-3]
      except:
        self.updateContents()
      if axis >3:
        order = self.maxN[:,axis-3]
      else:
        order = self.minP[axis]
      max = ray.origin[axis%3] - np.array([i.max[axis-3] for i in np.array(self.contents)[order]])
      min = ray.origin[axis%3] - np.array([i.min[axis-3] for i in np.array(self.contents)[order]])

      doit = 0 <= max 
      for c,i in enumerate(np.array(self.contents)[order]):
        if doit[c]:
          thatHit,thatRay = i.intersects(ray.copy(),closest=True)
          if thatHit and thatRay.tnear < ray.length:
            projLength = thatRay.tnear * np.abs(thatRay.direction[axis%3])
            doit = projLength >= max
            hit = thatHit
            ray.length = thatRay.tnear
    else:
      for i in self.contents:
        thatHit,thatRay = i.intersects(ray.copy(),closest=True)
        if thatHit and thatRay.tnear < ray.length:
          hit = thatHit
          ray.ccopy(thatRay)
          ray.length = thatRay.tnear

    return thisHit or hit,ray.copy()

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
    #import pdb;pdb.set_trace()
    hitPoint = (self.hit(ray,ok=True) - \
                 self.min)/self.extent
    ihitPoint = (hitPoint+0.5).astype(int)
    delta = hitPoint - ihitPoint
    axis = np.where(np.in1d(np.abs(delta),np.min(np.abs(delta))))[0]
    upDown = 2*ihitPoint[axis]-1
    direction = np.zeros(3)
    direction[axis] = upDown
    ray.localNormal = direction
    return direction

  def hit(self,ray,ok=False,diff=None):
    '''
    Calculate intersection point of ray with object
    
    ray : the input ray

    Options:
      ok   : set to True is ray.tnear is already calculated
             otherwise self.intersect(ray) is called
    '''
    diff = diff or -PyRatRayTol
    ok = ok or self.vintersect(ray)
    if ok: return ray.origin + (ray.length+diff) * ray.direction
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

  def rayTree(self,ray,closest=True):
    '''
    Return a ray tree
    '''
    from PyRatRay import PyRatRay
    sun = ray.sun
    hit = False
    # first try any infinite planes
    for i in self.planes:
      thatHit,thatRay = i.intersects(ray.copy(),closest=True)
      if thatHit and thatRay.tnear < ray.length:
        hit = thatHit
        ray.ccopy(thatRay)
        ray.length = thatRay.tnear
    # then try this object and its contents
    # noting that the infinite planes may have given
    # a ray.length we need to beat
    thisHit,thisRay = self.intersects(ray.copy(),closest=closest)
    if thisHit:
      ray.ccopy(thisRay)
      ray.length = thisRay.tnear
      hit = thisHit
      
    try:
      obj = ray.object
    except:
      obj = None 
    if hit and obj != None:
      try:
        n = ray.localNormal
      except:
        try:
          n = ray.localNormal = ray.object.surfaceNormal(ray)
        except:
          n = np.array([0,0,1.])
      #return True,ray,n
      # now see if we hit the sun
      ray.hitPoint = ray.length * ray.direction + ray.origin
      isTransmittance = dot(n,ray.direction) > 0
      #import pdb;pdb.set_trace()
      if isTransmittance:
        # transmittance
        sunRay = PyRatRay(ray.hitPoint+PyRatRayTol*ray.direction,sun)
      else:
        sunRay = PyRatRay(ray.hitPoint-PyRatRayTol*ray.direction,sun)
      try:
        sunRay.view = ray
        sunHit,sunThisRay = self.intersects(sunRay.copy(),closest=closest)
      except:
        # not sure what to do here 
        sunHit = True
      if not sunHit:
        return True,ray,n
      else:
        return False,ray,np.array([0,0,1.])
    else:
      return False,ray,np.array([0,0,1.])  

def test(base,tip,obj=None,name=None,type=None,file=None,info={},nAtTime=200):
  '''
  A simple test of the intersection algorithm
 
  A scan over an object is made and images produced
  in tests with the distances.

  '''

  from PyRatRay import PyRatRay
  from PyRatEllipsoid import PyRatEllipsoid
  from PyRatRay import PyRatRay
  from PyRatCylinder import PyRatCylinder
  from PyRatFacet import PyRatFacet
  from PyRatSpheroid import PyRatSpheroid
  from PyRatDisk import PyRatDisk
  from PyRatPlane import PyRatPlane
  from PyRatObjParser import PyRatObjParser
  from PyRatClone import PyRatClone
  import pylab as plt
  import matplotlib.cm as cm

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
  type = type.split('/')[-1]
  exec('from %s import %s'%(type,type))

  name = name or type[5:]
  obj = obj or eval('%s(base,tip,info=info)'%type)
  # ray direction
  direction = np.array([-1,-1.,-1])
  direction /= sqrt(dot(direction,direction))
 
  origin = -direction * 6.0 
  focalPoint = origin*1.5
  # sun direction
  sun = np.array([10,-10,20.])
  sun /= sqrt(dot(sun,sun))

  # image size
  size = (200,200)

  # ray origins
  #origin = np.array([0,0,4]).astype(float)
  o = origin.copy()
  ray = PyRatRay(o,direction)

  # dimensions of the image in physical units
  dimensions = [2,2]

  result0 = np.zeros(size)
  result1 = np.zeros(size)
  result2 = np.zeros(size)
  sys.stderr.write('from %s in direction %s\n'%(str(origin),str(direction)))
  sys.stderr.write('Name: %s\n'%name)
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

    sims = []
    l = float(size[0])
    index = []
    ray.sun = sun
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * 2.*(ix-size[0]*0.5)/size[0]
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * 2.*(iy-size[1]*0.5)/size[1]
        ray.length = PyRatBig
        d = (o - focalPoint)
        ray.direction = d/sqrt(dot(d,d))
        index.append((ix,iy))
        sims.append(ray.copy())

    sims = np.array(sims)
    index = np.array(index)
    results = []
    for i in xrange(0,len(sims),nAtTime):
      try:
        f = job_server.submit(obj.rayTreeMany,(sims[i:i+nAtTime],))  
        j = index[i:i+nAtTime]
      except:
        f = job_server.submit(obj.rayTreeMany,(sims[i:],))
        j = index[i:]
      results.append([j[:,0],j[:,1],f])

    if 'verbose' in info:
      sys.stderr.write('\nGathering results\n')
    l = size[0]*size[1]
    
    for c in xrange(len(results)):
      if 'verbose' in info and int(100.*(c+1)/l) % 5 == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*(c+1)/l))
      r = results[c]
      f = r[2]
      thisResult = f()
      # this returns True,ray,n
      try:
       ww = np.where(np.array(thisResult)[:,0])[0]
       iix = r[0][ww]
       iiy = r[1][ww]
       for j,val in enumerate(np.array(thisResult)[ww]):
          if val[0]:
            ix = iix[j]
            iy = iiy[j]
            ray = val[1]
            n = val[2]
            lambert = dot(n,sun)
            if lambert < 0: lambert = 0
            result0[size[0]-1-ix,size[1]-1-iy] = lambert
            result1[size[0]-1-ix,size[1]-1-iy] = ray.tnear
            result2[size[0]-1-ix,size[1]-1-iy] = ray.tfar
      except:
        pass
  else:
    l = size[0]
    ray.sun = sun
    for ix in xrange(size[0]):
      o[0] = origin[0] + dimensions[0] * 2.*(ix-size[0]*0.5)/size[0]
      if 'verbose' in info and int(100.*(ix+1)/l) % 5 == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*(ix+1)/l))
      for iy in xrange(size[1]):
        o[1] = origin[1] + dimensions[1] * 2.*(iy-size[1]*0.5)/size[1]
        ray.length = ray.tnear = ray.tfar =PyRatBig
        ray.origin = o
        d = (o - focalPoint) 
        ray.direction = d/sqrt(dot(d,d))
        try:
          hit,thisRay,n = obj.rayTree(ray.copy())
        except:
          hit = False
          try:
            print obj,ray
            print ray.copy()
          except:
            hit = False
        if hit:
          lambert = dot(n,sun)
          if lambert < 0: lambert = 0
          result0[size[0]-1-ix,size[1]-1-iy] = lambert
          result1[size[0]-1-ix,size[1]-1-iy] = thisRay.tnear
          result2[size[0]-1-ix,size[1]-1-iy] = thisRay.tfar
        else:
          missed = True
  if 'verbose' in info:
    sys.stderr.write('\nWriting results\n')
  plt.clf()
  plt.imshow(result0,interpolation='nearest',cmap=cm.Greys_r)
  if 'verbose' in info: sys.stderr.write('Mean: %f\n'%np.mean(result0))

  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests')
  plt.savefig('tests/PyRat%s.png'%name or file)
  plt.clf()
  plt.imshow(result1,interpolation='nearest',cmap=cm.Greys_r)
  plt.colorbar()
  plt.savefig('tests/PyRat%s-near.png'%name or file)
  plt.clf()
  plt.imshow(result2,interpolation='nearest',cmap=cm.Greys_r)
  plt.colorbar()
  plt.savefig('tests/PyRat%s-far.png'%name or file)


def main():
  '''
  A simple test of the bounding box algorithm
  
  A scan over a cube is made and 2 images produced
  tests/PyRatBox-near.png and tests/PyRatBox-far.png
  with the near and far distances 
  '''
  extent = np.array([1,1,1]).astype(float) * 3
  min = -extent/2.
  
  info = {'verbose':True}# ,'lad':3.0}
  test(min,extent,info=info,name='PyRatBox')

if __name__ == "__main__":
    main()

