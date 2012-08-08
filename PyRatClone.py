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
    self.planes = []

  def intersects(self,ray,closest=True):
    '''
    Call intersect but return ray as well

    and do hierarchical intersections
    '''
    from PyRatBox import PyRatBox
    import numpy as np
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

    try:
      transformed_ray.direction = np.array(np.matrix(ray.direction) * self.matrix.T).flatten()
      mod_ray_direction=sqrt(dot(transformed_ray.direction,transformed_ray.direction))
      transformed_ray.direction /= mod_ray_direction
      transformed_ray.origin = np.array(np.matrix(transformed_ray.origin) * self.matrix.T).flatten()
      # to account for scaling effects
      transformed_ray.origin /= (mod_ray_direction*mod_ray_direction)
      transformed_ray.length /= mod_ray_direction
      transformed_ray.tnear /= mod_ray_direction
      transformed_ray.tfar /= mod_ray_direction
      transformed_ray.big /= mod_ray_direction
    except:
      pass
    # with a clone, we are not interested in the intersection other
    # than the transformation we apply
    # contents[0] will always be the contents of a clone as it
    # points to a group which must be a bounding box
    #import pdb;pdb.set_trace()
    hit,thisRay  = PyRatBox.intersects(self.contents[0],transformed_ray,closest=closest)
    if hit:
      try:
        #import pdb;pdb.set_trace()
        thisRay.localNormal = np.array(np.matrix(thisRay.localNormal) * self.matrix).flatten()
        mod = sqrt(dot(thisRay.localNormal,thisRay.localNormal))
        thisRay.localNormal /= mod
      except:
        pass

      try:
        thisRay.length *= self.scale
        thisRay.tnear *= self.scale
        thisRay.tfar *= self.scale
        thisRay.big *= self.scale
      except:
        pass

      try:
        #import pdb;pdb.set_trace()
        #thisRay.length *= mod_ray_direction
        thisRay.origin = ray.origin
        thisRay.direction = ray.direction
      except:
        pass
      ray.ccopy(thisRay)
    return hit,ray

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
  from PyRatObjParser import PyRatObjParser
  from PyRatPlane import PyRatPlane
  from tempfile import NamedTemporaryFile
  import os

  data = '''
  !{
  !{
  v 0 0 0
  v 0 0 1
  plane -1 -2
  #define objects
  g group
  box 2 2 0.5 1 1 1
  v 0 0 1
  sph -1 0.5
  v 2 1 2
  sph -1 0.5
  v -1 1 0
  ell -1 0.2 0.3 1
  !}
  clone -0.5 1.5 0 Rz 45 Ry 10 Rx -5 group
  !}
  '''
  f = NamedTemporaryFile(delete=False)
  f.write(data)
  f.close()
  
  p = PyRatObjParser(f.name)
  p.root.planes = p.infinitePlane
 
  name = str(globals()['__file__'].split('/')[-1].split('.')[0])  
  test(np.zeros(3),None,obj=p.root,info={},type=name,nAtTime=100*100/20)
  os.unlink(f.name)

if __name__ == "__main__":
    main()
