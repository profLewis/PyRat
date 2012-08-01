import numpy as np
from PyRatBox import *

class PyRatObjParser(object):
  '''
  Parser for extended wavefront format
  '''
  def __init__(self,filename,verbose=True,reportingFrequency=10):
    self.setupDictionary()

    self.verbose=verbose
    self.reportingFrequency = reportingFrequency

    none = np.zeros(3)
    self.top = PyRatBox(none,none,material=None)
    self.root = self.top
    self.error = self.top.error

    self.group = {}

    self.stack = []
    self.point = []
    self.nPoints = 0
    self.read(filename)
    self.root = self.reconcile(self.root)
    # reset the stack
    self.point = np.array(self.point)
    self.verbose=verbose
    self.reportingFrequency = 10

  def combine(self,other):
    '''
    combine other with self
    '''

  def reconcile(self,bbox):
    '''
    Reconcile world object:
     - update bboxes (min, max)
     - delete Box objects if they only contain an empty box
       (this gets rid of spurious bbox info)
     - updates the surface area (size) info for the box
    '''
    # if there is only one object and its a box and its empty
    # then delete it
    if len(bbox.contents) == 1:
      if type(bbox.contents[0]) == PyRatBox:
        if bbox.contents[0].empty:
          return self.reconcile(bbox.contents[0])
    # scan over the world object
    for c,i in enumerate(bbox.contents):
      bbox.contents[c] = self.reconcile(i)
    try:
      if not bbox.empty:
        min = bbox.min
        max = bbox.max
        size = bbox.size
      bbox.size = np.sum([i.size for i in bbox.contents])
      bbox.min = np.min(np.array([i.min for i in bbox.contents]),axis=0)
      bbox.max = np.max(np.array([i.max for i in bbox.contents]),axis=0)
      if not bbox.empty:
        # if this isnt an empty box then add in its contents
        bbox.min = np.min(np.array([bbox.min,min]),axis=0)
        bbox.max = np.max(np.array([bbox.max,max]),axis=0)
        bbox.size += size
    except:
      pass
    return bbox

  def setupDictionary(self):
    '''
    Set up the parser dictionary
    '''
    self.dict = {\
      '#define':self.define,\
      '!{':self.openBox,\
      '!}':self.closeBox,\
      'usemtl':self.usemtl,\
      'mtllib':self.mtllib,\
      'disk':self.disk,\
      'plane':self.plane,\
      'clone':self.clone,\
      'g':self.g,\
      'v':self.v,\
      'ell':self.ell,\
      'sph':self.sph,\
      'f':self.f}

  def define(self,cmd):
    self.top.invisible=True
    pass

  def clone(self,cmd):
    '''
    clone
    '''
    try:
      pass
    except:
      self.error('could not interpret clone line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def ell(self,cmd):
    '''
    Ellipsoid

    ell base rx ry rz
    '''
    from PyRatEllipsoid import PyRatEllipsoid
    try:
      this = int(cmd[1])
      if this < 0:
        this += self.nPoints
      self.contents.append(\
           PyRatEllipsoid(self.point[this[1]],\
             np.array([this[2],this[3],this[4]]).astype(float),\
             material=self.top.material))
    except:
      self.error('could not interpret "ell" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def cyl(self,cmdi,info={}):
    '''
    Cylinder
  
    cyl base tip radius
    '''
    from PyRatCylinder import PyRatCylinder
    try:
      this = np.array([cmd[1],cmd[2]]).astype(float)
      this[this<0] += self.nPoints
      info.update({'radius':float(this[3])})
      self.top.contents.append(\
           PyRatCylinder(self.point[this[1]],self.point[this[1]],info=info,\
                         material=self.top.material))
    except:
      self.error('could not interpret "cyl" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def ccyl(self,cmd):
    '''
    Capped (closed) Cylinder
 
    ccyl base tip radius
    '''
    self.cyl(cmd,info={'caps':True})

  def sph(self,cmd):
    '''
    Spheroid

    sph centre radius
    '''
    from PyRatSpheroid import PyRatSpheroid
    try:
      this = int(cmd[1])
      if this < 0:
        this += self.nPoints
      self.top.contents.append(\
           PyRatSpheroid(self.point[this[1]],float(this[2]),\
                         material=self.top.material))
    except:
      self.error('could not interpret "sph" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def plane(self,cmd):
    '''
    Infinite plane

      plane normal centre

    '''
    from PyRatPlane import PyRatPlane
    try:
      this = np.array(cmd[1:3]).astype(int)
      this[this<0] += self.nPoints
      self.root.contents.append(\
           PyRatPlane(self.point[this[1]],self.point[this[0]],\
                      material=self.top.material))
    except:
      self.error('could not interpret "plane" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def usemtl(self,cmd):
    '''
    Set the current material
    '''
    try:
      self.top.material = cmd[1]
    except:
      self.error('could not interpret "usemtl" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def read(self,filename):
    '''
    Read extended wavefront file filename
    '''
    import sys
    self.filename = filename
    try:
      this = open(filename).read().split('\n')
      l = float(len(this))
      reportingFrequency = int(l/self.reportingFrequency)
      for i,line in enumerate(this):
        if self.verbose and i%reportingFrequency == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*i/l))
        if len(line):
          self.lineNumber = i+1
          self.parseLine(line.split())
    except:
      self.error('Unable to read wavefront file %s'%filename)

  def parseLine(self,cmd):
    '''
    Parse a single line of input
    '''  
    try:
      self.dict[cmd[0]](cmd)
    except:
      self.error('could not interpret line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def openBox(self,cmd):
    '''
    Open a bounding box

    Start a new container box and push on stack
    '''
    from PyRatBox import PyRatBox
    self.top.contents.append(PyRatBox(np.zeros(3),None))
    self.stack.append(self.top)
    self.top = self.top.contents[-1]

  def closeBox(self,cmd):
    '''
    Close a bounding box
    '''
    self.top = self.stack.pop()
 
  def mtllib(self,cmd):
    '''
    Material library
    '''
    self.root.materialFile = cmd[1] 

  def g(self,cmd):
    '''
    Group:

      g group name
    '''
    try:
      groupName = ' '.join(cmd[1:])
    except:
      self.error('could not interpret "g" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))
    # now we need to reference for this object
    # This should be the last item on the stack but if not its the root
    try:
      this = self.stack[-1]
    except:
      this = top.root
    self.group[groupName] = this

  def v(self,cmd):
    '''
    Vertex:

      v x y z
    '''
    try:
      self.point.append(np.array(cmd[1:4]).astype(float))
      self.nPoints += 1
    except:
      self.error('could not interpret "v" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def disk(self,cmd):
    '''
    Disk
    
      disk centre normal radius
    '''
    from PyRatDisk import PyRatDisk
    try:
      this = np.array(cmd[1:3]).astype(int)
      this[this<0] += self.nPoints
      radius = float(cmd[3])
      self.top.contents.append(\
           PyRatDisk(self.point[this[0]],self.point[this[1]],\
                     material=self.top.material,info={'radius':radius}))
    except:
      self.error('could not interpret "disk" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def f(self,cmd):
    '''
    Facet:

      f n1 n2 n3
    '''
    from PyRatFacet import PyRatFacet
    try:
      this = np.array(cmd[1:4]).astype(int)
      this[this<0] += self.nPoints
      self.top.contents.append(\
           PyRatFacet(np.array([self.point[this[0]],self.point[this[1]],\
                      self.point[this[2]]]),\
                      material=self.top.material))
    except:
      self.error('could not interpret "f" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))


