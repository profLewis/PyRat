import numpy as np

class PyRatRay(object):
  def __init__(self,origin,direction):
    self.origin = origin
    self.direction = direction
    self.tnear = -1e20
    self.tfar = 1e20
    self.length = 1e20
    self.hitPoint = np.zeros(3)

  def copy(self):
    tnear = self.tnear
    tfar = self.tfar
    length = self.length
    hitPoint = self.hitPoint
    new = PyRatRay(self.origin.copy(),self.direction.copy())
    new.tnear = tnear
    new.tfar = tfar
    new.length = length
    new.hitPoint = hitPoint
    return new 
