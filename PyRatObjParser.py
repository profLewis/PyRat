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
    self.top = PyRatBox(none,none)
    self.point = []
    self.nPoints = 0
    self.read(filename)
    self.point = np.array(self.point)
    self.verbose=verbose
    self.reportingFrequency = 10

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
    pass

  def clone(self,cmd):
    '''
    Infinite plane
    '''
    try:
      pass
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

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
      self.top.contents.append(\
           PyRatEllipsoid(self.point[this[1]],\
             np.array([this[2],this[3],this[4]]).astype(float)))
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

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
           PyRatCylinder(self.point[this[1]],self.point[this[1]],info=info))
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

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
           PyRatSpheroid(self.point[this[1]],float(this[2])))
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

  def plane(self,cmd):
    '''
    Infinite plane

      plane normal centre

    '''
    from PyRatPlane import PyRatPlane
    try:
      this = np.array(cmd[1:3]).astype(int)
      this[this<0] += self.nPoints
      self.top.contents.append(\
           PyRatPlane(self.point[this[1]],self.point[this[0]]))
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

  def usemtl(self,cmd):
    '''
    Set the current material
    '''
    try:
      self.top.currentMaterial = cmd[1]
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))


  def read(self,filename):
    '''
    Read extended wavefront file filename
    '''
    import sys
    try:
      this = open(filename).read().split('\n')
      l = float(len(this))
      reportingFrequency = int(l/self.reportingFrequency)
      for i,line in enumerate(this):
        if self.verbose and i%reportingFrequency == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*i/l))
        if len(line):
          self.parseLine(line.split())
    except:
      self.top.error('Unable to read wavefront file %s'%filename)

  def parseLine(self,cmd):
    '''
    Parse a single line of input
    '''  
    try:
      self.dict[cmd[0]](cmd)
    except:
      self.top.error('could not interpret line %s'%' '.join(cmd))

  def openBox(self,cmd):
    '''
    Open a bounding box
    '''
    pass

  def closeBox(self,cmd):
    '''
    Close a bounding box
    '''
    pass

 
  def mtllib(self,cmd):
    '''
    Material library
    '''
    self.top.materialFile = cmd[1] 

  def g(self,cmd):
    '''
    Group:

      g group name
    '''
    try:
      self.top.groupName = ' '.join(cmd[1:])
    except:
      self.top.error('error interpreting g object: %s'%' '.join(cmd)) 


  def v(self,cmd):
    '''
    Vertex:

      v x y z
    '''
    try:
      self.point.append(np.array(cmd[1:4]).astype(float))
      self.nPoints += 1
    except:
      self.top.error('error interpreting v object: %s'%' '.join(cmd))  

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
           PyRatDisk(self.point[this[0]],self.point[this[1]],info={'radius':radius}))
    except:
      self.top.error('error interpreting disk object: %s'%' '.join(cmd))

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
           PyRatFacet(np.array([self.point[this[0]],self.point[this[1]],self.point[this[2]]])))
    except:
      self.top.error('error interpreting f object: %s'%' '.join(cmd))  
