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
    try:
      #import pdb;pdb.set_trace()
      transformed_ray = ray.copy()
      transformed_ray.origin = np.array(np.matrix(ray.origin) * self.matrix).flatten()
      transformed_ray.direction = np.array(np.matrix(ray.direction) * self.matrix).flatten()
      mod_ray_direction=np.dot(transformed_ray.direction,transformed_ray.direction)
      transformed_ray.direction /= mod_ray_direction
    except:
      transformed_ray = ray
    try:
      transformed_ray.origin += self.offset
    except:
      pass
    try:
      transformed_ray.sourceOfRay = ray.sourceOfRay
    except:
      pass
    return PyRatBox.intersects(self,\
                               transformed_ray,closest=closest)

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

  min = [-0.5,-0.5,2]
  extent = [1,2,1]
  info = {'verbose':True,'lad':3.0}
  box = PyRatBox(min,extent,info=info)

  centre = [0,-0.5,2.5]
  radius = 0.25
  info = {'verbose':True}
  sph = PyRatSpheroid(centre,radius,info=info)

  clone = PyRatClone(np.zeros(3),None)
  clone.thisGroup = None
  clone.offset = np.array([1.0,0,0])
  clone.matrix = np.eye(3)
  clone.contents = [box,sph]

  name = str(globals()['__file__'].split('.')[0])  
  test(min,extent,obj=clone,info=info,type=name)


if __name__ == "__main__":
    main()
