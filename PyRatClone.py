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

from PyRatBox import *
import numpy as np
import sys

class PyRatClone(PyRatBox):
  '''
  PyRatClone: A PyRat object 

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
    PyRatBox.__init__(self,base,extent,contents=contents,material=material,info=info)

  def intersects(self,ray,closest=True):
    '''
    Call intersect but return ray as well

    and do hierarchical intersections
    '''
    # transform ray
    transformed_ray = ray.copy()

    try:
      transformed_ray.sourceOfRay = ray.sourceOfRay
    except:
      pass

    try:
      transformed_ray.origin -= self.offset
    except:
      pass

    #import pdb;pdb.set_trace()
    try:
      transformed_ray.direction = np.array(np.matrix(ray.direction) * self.matrix.T).flatten()
      mod_ray_direction=np.sqrt(np.dot(transformed_ray.direction,transformed_ray.direction))
      transformed_ray.direction /= mod_ray_direction

      #import pdb;pdb.set_trace()
      transformed_ray.origin = np.array(np.matrix(transformed_ray.origin) * self.matrix.T).flatten()
      # to account for scaling effects
      transformed_ray.origin /= (mod_ray_direction*mod_ray_direction)
      transformed_ray.length /= mod_ray_direction
    except:
      pass

    hit,thisRay  = PyRatBox.intersects(self,transformed_ray,closest=closest)
    if hit:
      try:
        thisRay.object.localNormal = np.array(np.matrix(thisRay.object.localNormal) * self.matrix).flatten()
        mod = np.sqrt(np.dot(thisRay.object.localNormal,thisRay.object.localNormal))
        thisRay.object.localNormal /= mod
        thisRay.length *= mod_ray_direction
        thisRay.origin = ray.origin
        thisRay.direction = ray.direction
      except:
        pass
    return hit,thisRay

def main():
  '''
  A simple test of the clone algorithm
  
  A scan over a cube is made and 2 images produced
  tests/PyRatBox-near.png and tests/PyRatBox-far.png
  with the near and far distances 

  Here we also demonstrate that setting info['lad']
  makes the object volumetric
  '''
  from PyRatBox import PyRatBox
  from PyRatSpheroid import PyRatSpheroid
  from PyRatClone import PyRatClone
  from PyRatEllipsoid import PyRatEllipsoid

  min = [-0.5,-0.5,2]
  extent = [1,2,1]
  info = {'verbose':True,'lad':3.0}
  box = PyRatBox(min,extent,info=info)

  centre = [0,-0.5,2.5]
  radius = 0.25
  info = {'verbose':True}
  sph = PyRatSpheroid(centre,radius,info=info)

  base = np.array(min) + np.array(extent) *0.5
  radius = [0.25,0.5,0.5]
  info = {'verbose':True,'lad':10.}
  ell = PyRatEllipsoid(base,radius,info=info)

  clone = PyRatClone(np.zeros(3),None)
  clone.thisGroup = None
  clone.offset = np.array([0.0,0.0,-0.5])
  clone.matrix = np.eye(3)
  c = np.cos(45*np.pi/180.)
  s = np.sin(45*np.pi/180.)
  clone.matrix[1,1] = clone.matrix[0,0] = c
  clone.matrix[0,1] = -s
  clone.matrix[1,0] = s

  clone.contents = [box,sph,ell]

  name = str(globals()['__file__'].split('.')[0])  
  test(min,extent,obj=clone,info=info,type=name)


if __name__ == "__main__":
    main()
