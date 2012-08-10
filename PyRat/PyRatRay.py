#!/usr/bin/env python
import numpy as np
from PyRatBox import PyRatBig

class PyRatRay(object):
  def __init__(self,origin,direction):
    self.origin = origin
    self.direction = direction
    self.tnear = -1e20
    self.tfar = 1e20
    self.length = 1e20
    self.hitPoint = np.zeros(3)
    self.object = None
    self.rays = None
    self.sun = np.array([0,0,1.])
    self.big = PyRatBig

  def ccopy(self,new):
    tnear = new.tnear
    tfar = new.tfar
    length = new.length
    hitPoint = new.hitPoint
    self.origin = new.origin
    self.direction = new.direction
    self.tnear = tnear
    self.tfar = tfar
    self.length = length
    self.hitPoint = hitPoint
    self.object = new.object
    self.rays = new.rays
    self.sun = new.sun
    self.localNormal = new.localNormal

  def copy(self):
    tnear = self.tnear
    tfar = self.tfar
    length = self.length
    hitPoint = self.hitPoint
    new = PyRatRay(np.array([self.origin[0],self.origin[1],self.origin[2]]),\
                   np.array([self.direction[0],self.direction[1],self.direction[2]]))
    new.tnear = tnear
    new.tfar = tfar
    new.length = length
    new.hitPoint = hitPoint
    new.object = self.object
    new.rays = self.rays
    new.sun = self.sun
    try:
      new.localNormal = self.localNormal
    except:
      new.localNormal = np.array([0.,0.,1.])
    return new 

  def error(self,msg):
    '''
    Error reporting
    '''
    import sys
    sys.stderr.write(msg+'\n')

  def report(self,level=0):
    '''
    Report on the ray information
    '''
    buff = '_'*level
    self.error('%s{'%buff)
    self.error('%s origin   \t%s'%(buff,str(self.origin)))
    self.error('%s direction\t%s'%(buff,str(self.direction)))
    self.error('%s tNear    \t%s'%(buff,str(self.tnear)))
    self.error('%s tFar     \t%s'%(buff,str(self.tfar)))
    self.error('%s length   \t%s'%(buff,str(self.length)))
    self.error('%s hit point\t%s'%(buff,str(self.hitPoint)))
    self.error('%s object   \t%s'%(buff,str(self.object)))
    try:
      self.error('%s normal   \t%s'%(buff,str(self.localNormal)))
    except:
      try:
        self.localNormal = self.object.surfaceNormal(self)
        self.error('%s normal   \t%s'%(buff,str(self.localNormal)))
      except:
        pass
    try:
      self.error('%s %s'%(buff,'sun:'))
      self.sun.report(level=level+1)
    except:
      pass
    try:
      self.error('%s %s'%(buff,'view:'))
      self.view.report(level=level+1)
    except:
      pass
    try:
      for ray in self.rays:
        ray.report(level=level+1)
    except:
      pass
    self.error('%s}'%buff)
     
