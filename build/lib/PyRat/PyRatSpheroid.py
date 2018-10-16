#!/usr/bin/env python
from  PyRat.PyRatEllipsoid import *

class PyRatSpheroid(PyRatEllipsoid):
  '''
  PyRatSpheroid: A PyRat object

  P. Lewis 30/6/2012

  For a sphere
  '''
  def __init__(self,base,radius,contents=None,material=None,info=None):
    '''
    Load the object:

      base : vector that define the spheroid centre
      radius : vector (for ellipsoid) or scalar for sphere

      OPTIONAL:
      info['coords']:
        3 x 2 vectors that define a coordinate for each vertex

      Options:
        contents   : a list that may contain other objects
                     or None
        material   : a dictionary defining materials
        info       : placeholder

      Note the the centre of the object is stored as self.base

    '''
    base[2] -= radius
    PyRatEllipsoid.__init__(self,base,radius,contents=contents,material=material,info=info)

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
    hitPoint = self.hit(ray,ok=True)
    v1 = (hitPoint - self.centre)/self.radius
    ray.localNormal = v1/sqrt(dot(v1,v1))
    return v1

  def draw(self,matrix=None,offset=None,scale=1.0):
    '''
    mayavi/tvtk drawing method
    '''
    try:
      from enthought.tvtk.tools import visual
    except:
      return None
    sph = visual.Sphere(pos=tuple(modify(self.centre,matrix,offset)),\
         radius=self.radius[0]*scale)
    return sph



  def tesselate(self,N=8):
    '''
    This method forms a new (triangulated) version of the object
    with N segments.

    It returns a PyRatBox object which contains the set of N PyRatFacet
    objects.

    Options:
      N   : number of segments
    '''
    phi, theta = np.mgrid[0:np.pi:complex(0,N), 0:2*np.pi:complex(0,N)]
    x = np.sin(phi) * np.cos(theta) * self.radius[0] * 0.5 + self.centre[0]
    y = np.sin(phi) * np.sin(theta) * self.radius[1] * 0.5 + self.centre[0]
    z = np.cos(phi) * self.radius[2] * 0.5 + self.centre[0]
    return np.array([x,y,z])

  def tesselate1(self,N=8):
    '''
    This method forms a new (triangulated) version of the object
    with N segments.

    It returns a PyRatBox object which contains the set of N PyRatFacet
    objects.

    Options:
      N   : number of segments

    Algorithm from:
      http://paulbourke.net/miscellaneous/sphere_cylinder/

    '''
    from PyRatBox import PyRatBox
    from PyRatFacet import PyRatFacet
    centre = self.centre
    radius = self.radius

    p1 = np.array([1.0,1.0,1.0])
    p2 = np.array([-1.0,-1.0,1.0])
    p3 = np.array([1.0,-1.0,-1.0])
    p4 = np.array([-1.0,1.0,-1.0])
    p1 = normalise(p1)
    p2 = normalise(p2)
    p3 = normalise(p3)
    p4 = normalise(p4)

    facets = [[p1,p2,p3],[p2,p1,p4],[p2,p4,p3],[p1,p3,p4]]
    n = 4
    for i in range(N):
      nstart = n
      for j in range(nstart):
        # Create initially copies for the new facets
        facets.append([i.copy() for i in facets[j]])
        facets.append([i.copy() for i in facets[j]])
        facets.append([i.copy() for i in facets[j]])
        # Calculate the midpoints
        p1 = (facets[j][0]+facets[j][1])*0.5
        p2 = (facets[j][1]+facets[j][2])*0.5
        p3 = (facets[j][2]+facets[j][0])*0.5
        # Replace the current facet
        facets[j][1] = p1
        facets[j][2] = p3
        # Create the changed vertices in the new facets
        facets[n][0]   = p1
        facets[n][2]   = p2
        facets[n+1][0] = p3
        facets[n+1][1] = p2
        facets[n+2] = [p1,p2,p3]
        n += 3

    output = []
    for j in range(n):
      facets[j] = [normalise(facets[j][0])*radius+centre,\
                   normalise(facets[j][1])*radius+centre,\
                   normalise(facets[j][2])*radius+centre] 
      output.append(PyRatFacet(facets[j],None)) 
    return output

def main():
  '''
  A simple test of the Spheroid algorithm
 
  A scan over a Spheroid is made and an image produced
  tests/PyRatSpheroid-near.png with the distances.

  It should be ~1 in the centre (since the camera is located at z=4)
  but this is a volumetric example.

  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  from PyRatBox import test

  base = np.array([0,0,0.])
  radius = 1.0
  info = {'verbose':True}#,'lad':3.0}

  name = str(globals()['__file__'].split('/')[-1].split('.')[0])
  test(base,radius,info=info,type=name,name=name[5:])

if __name__ == "__main__":
    main()


