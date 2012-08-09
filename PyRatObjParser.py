#!/usr/bin/env python 
import numpy as np
from PyRatBox import *

try:
  from mayavi import mlab
  from tvtk.api import tvtk
  #from mayavi.scripts import mayavi2
  hasGL = True
except:
  hasGL = False

class PyRatObjParser(object):
  '''
  Parser for extended wavefront format
  '''
  def __init__(self,filename,GL=False,verbose=True,reportingFrequency=10):
    self.setupDictionary()
    self.GL = hasGL and GL
    self.verbose=verbose
    self.reportingFrequency = reportingFrequency

    none = np.zeros(3)
    self.top = PyRatBox(none,none,material=None)
    self.top.invisible = True
    self.root = self.top
    self.root.defined = False
    self.root.invisible = True
    self.root.empty = False

    self.infinitePlane = []

    self.error = self.top.error

    # this is where we keep groups
    # we have to reconcile all of these as well as
    # the main object in case some of them
    # are #define'd
    self.group = {}

    self.stack = []
    self.point = []
    self.nPoints = 0
    if self.verbose:
      sys.stderr.write('Reading data from %s\n'%filename)
    self.read(filename)
    if self.verbose:
      sys.stderr.write('\n ... sorting bbox contents\n')
    if not self.root.isDefined():
      self.root = self.reconcile(self.root,1)
    for i in self.group.keys():
      self.group[i] = self.reconcile(self.group[i],1,defined=True) 
    # reset the stack
    self.point = np.array(self.point)
    self.verbose=verbose
    self.reportingFrequency = 10
    if self.GL:
      self.triangleN = []
      self.triangleData = []
      if self.verbose:
        sys.stderr.write('\n ... sorting GL representation\n')
      self.GLobjects = []
      self.loadGL(self.root,np.eye(3),np.matrix(np.zeros(3)))

  def loadGL(self,bbox,Mmatrix,Ooffset):
    '''
    Load a GL representation
    '''
    from PyRatClone import PyRatClone
    from PyRatSpheroid import PyRatSpheroid
    from PyRatCylinder import PyRatCylinder
    from PyRatFacet import PyRatFacet
    matrix = Mmatrix.copy()
    offset = Ooffset.copy()
    if bbox == self.root:
      pass
    if type(bbox) == PyRatClone:
      try:
        matrix = matrix * bbox.matrix
        offset = offset * bbox.matrix
      except:
        pass
      try:
        offset += bbox.offset
      except:
        pass

    try:
      if type(bbox) != PyRatFacet:
        self.GLobjects.append(bbox.draw(matrix=matrix,offset=offset,scale=1.0))
      else:
        data,triangles = bbox.draw(matrix=matrix,offset=offset,scale=1.0)
        [self.triangleData.append(data[j]) for j in xrange(3)]
        [self.triangleN.append(triangles[j]+len(self.triangleN)) \
                             for j in xrange(3)]
    except:
      pass
    for i in bbox.contents:
      self.loadGL(i,matrix,offset)

    return None

  def dump(self,filename):
    '''
    Dump a numpy representation
    '''
    np.savez(filename,self=self.root)

  def load(self,filename):
    '''
    unDump a numpy representation
    '''
    return np.load(filename)


  def combine(self,other):
    '''
    combine other with self
    '''

  def reconcile(self,bbox,level,defined=False,ignoreVisit=False):
    '''
    Reconcile world object:
     - update bboxes (min, max)
     - delete Box objects if they only contain an empty box
       (this gets rid of spurious bbox info)
     - updates the surface area (size) info for the box
    '''
    from PyRatClone import PyRatClone
    try:
      if not ignoreVisit and bbox.visited:
        return bbox
    except:
      bbox.visited = True
    if defined and bbox.isDefined():
      # if it has a #define, remove from list
      return None
    if type(self) == PyRatClone:
      print 'clone'
    # if there is only one object and its a box and its empty
    # then delete it
    if len(bbox.contents) == 1:
      if type(self) == PyRatBox and type(bbox.contents[0]) == PyRatBox:
        if bbox.contents[0].invisible and self.invisible:
          return self.reconcile(bbox.contents[0],level,ignoreVisit=ignoreVisit)
    # scan over the world object
    for c,i in enumerate(bbox.contents):
      bbox.contents[c] = self.reconcile(i,level+1,ignoreVisit=ignoreVisit)
    #  if type(i) == PyRatBox and bbox.contents[c].invisible:
    #    # remove it
    #    bbox.contents.pop(c)
    try:
      if self.verbose:
        arrowClear = ['\b']*35
        arrowClear = ''.join(arrowClear)
        arrow = ['-']*level
        arrow = '%20s>(%03d)'%(''.join(arrow),level)
        sys.stderr.write(arrowClear + arrow)
      if not bbox.empty and type(bbox) != PyRatClone:
        min = bbox.min
        max = bbox.max
        size = bbox.size
      bbox.size = np.sum([i.size for i in bbox.contents])
      bbox.min = np.min(np.atleast_2d([i.min for i in bbox.contents]),axis=0)
      bbox.max = np.max(np.atleast_2d([i.max for i in bbox.contents]),axis=0)
      if len(bbox.min) == 0:
        bbox.min = np.zeros(3) + PyRatBig
      if len(bbox.max) == 0:
        bbox.max = np.zeros(3) + -PyRatBig
      if not bbox.empty and type(bbox) != PyRatClone:
        # if this isnt an empty box then add in its contents
        bbox.min = np.min(np.atleast_2d([bbox.min,min]),axis=0)
        bbox.max = np.max(np.atleast_2d([bbox.max,max]),axis=0)
        bbox.size += size
      bbox.updateBbox()
    except:
      pass
    if bbox.size > 0:
      bbox.empty = False
    #bbox.base = bbox.min
    bbox.extent = bbox.max - bbox.min
    # get rid of any None entries and zero sized objects
    for c,i in enumerate(bbox.contents):
      if type(i) == PyRatBox and i.empty:
        bbox.contents[c] = None
    if len(bbox.contents):
      bbox.contents = [j for j in np.array(bbox.contents)[np.array([i != None for i in bbox.contents])]]
    return bbox

  def setupDictionary(self):
    '''
    Set up the parser dictionary
    '''
    self.dict = {\
      '#':self.comment,\
      '#define':self.define,\
      'box':self.box,\
      '!{':self.openBox,\
      '!}':self.closeBox,\
      'usemtl':self.usemtl,\
      'mtllib':self.mtllib,\
      'disk':self.disk,\
      'plane':self.plane,\
      'clone':self.clone,\
      'transformation_matrix':self.transformation_matrix,\
      'g':self.g,\
      'v':self.v,\
      'ell':self.ell,\
      'sph':self.sph,\
      'cyl':self.cyl,\
      'ccyl':self.ccyl,\
      'f':self.f}

  def comment(self,cmd):
    pass

  def define(self,cmd):
    self.top.invisible=True
    self.top.defined = True
    pass

  def clone(self,cmd):
    '''
    clone

    clone dx dy dz [Rz] name
    '''
    from PyRatBox import PyRatBox
    from PyRatClone import PyRatClone
    try:
      try:
        this = cmd.split[1:]
      except:
        this = cmd[1:]
      # look up groups
      try:
        offset = np.array(this[:3]).astype(float)
        this = cmd[4:]
      except:
        offset = np.zeros(3)
        
      thisGroup = None
      thisGroupName = None
      matrix = np.matrix(np.eye(3))
      done = np.zeros_like(np.array(this)).astype(bool)
      for i in xrange(len(this)):
        if not done[i]:
          if this[i] == 'Rx':
            theta = float(this[i+1])*np.pi/180.
            done[i:i+2] = True
            if theta != 0:
              c = np.cos(theta)
              s = np.sin(theta)
              m = np.eye(3)
              m[1,1] = m[2,2] = c
              m[2,1] = s
              m[1,2] = -s
              matrix = matrix * np.matrix(m)
          elif this[i] == 'Ry':
            theta = float(this[i+1])*np.pi/180.
            done[i:i+2] = True
            if theta != 0:
              c = np.cos(theta)
              s = np.sin(theta)
              m = np.eye(3)
              m[2,2] = m[0,0] = c
              m[0,2] = s
              m[2,0] = -s
              matrix = matrix * np.matrix(m)
          elif this[i] == 'Rz':
            theta = float(this[i+1])*np.pi/180.
            done[i:i+2] = True
            if theta != 0:
              c = np.cos(theta)
              s = np.sin(theta)
              m = np.eye(3)
              m[0,0] = m[1,1] = c
              m[0,1] = -s
              m[1,0] = s
              matrix = matrix * np.matrix(m)
          elif this[i] == 'Transform':
          # arbitrary 3 x 3 rotation
            m = \
               np.matrix(np.array(this[i+1:i+9+1]).astype(float).reshape((3,3)))
            done[i+1:i+9+1] = True  
            matrix = matrix * np.matrix(m)
          else:
            # maybe its a grp name
            possible = ' '.join(this[i:])
            try:
              thisGroup = self.group[possible] 
              thisGroupName = possible
              done[:] = True
            except:
              # maybe its a z rotation
              try:
                theta = float(this[i])*np.pi/180.
                done[i] = True
                if theta != 0:
                  c = np.cos(theta)
                  s = np.sin(theta)
                  m = np.eye(3)
                  m[0,0] = m[1,1] = c
                  m[0,1] = -s
                  m[1,0] = s
                  matrix = matrix * np.matrix(m)
              except:
                pass
      # NB dummy extent at this point        
      clone = PyRatClone(np.zeros(3),None)
      clone.empty = False
      clone.invisible = True
      if (matrix != np.matrix(np.eye(3))).any():
        clone.matrix = matrix
      clone.thisGroup = thisGroup
      clone.thisGroupName = thisGroupName
      if np.abs(offset).sum() > 0:
        clone.offset = offset
      clone.contents = [clone.thisGroup]
      self.top.contents.append(clone)
      self.verbose == 2 and sys.stderr.write('<Clone %s>'%thisGroupName)      
    except:
      self.error('could not interpret clone line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def transformation_matrix(self,cmd):
    '''
    transformation_matrix

    '''
    from PyRatBox import PyRatBox
    from PyRatClone import PyRatClone
    try:
      matrix = self.matrix
    except:
      matrix = np.matrix(np.eye(3))
    try:
      offset = self.offset
    except:
      offset = np.zeros(3)

    try:
      try:
        this = cmd.split[1:]
      except:
        this = cmd[1:]
      clone = PyRatClone(np.ones(3),np.ones(3))
      for i,n in enumerate(this):
        noGood = False
        if n == 'scale':
          scale = float(this[i+1])
          matrix2,offset2 = self.load_scaling_matrix4(scale)        
        elif n == 'transformation_matrix':
          matrix2,offset2 = np.matrix(np.array(this[i+1:i+1+9])\
                                      .reshape((3,3)).astype(float))
        elif n == 'scale_differential':
          fix_point = np.array(this[i+1:i+1+3]).astype(float)
          matrix2,offset2 = self.load_differential_scaling_matrix4(fix_point)  
        elif n == 'translate':
          offset3 = np.array(this[i+1:i+1+3]).astype(float)
          matrix2,offset2 = self.load_translation_matrix4(offset3)
        elif n == 'rotate_about_x_axis':
          rot = float(this[i+1])
          matrix2,offset2 = self.load_x_axis_rotation_matrix4(rot)
        elif n == 'rotate_about_y_axis':
          rot = float(this[i+1])
          matrix2,offset2 = self.load_y_axis_rotation_matrix4(rot)
        elif n == 'rotate_about_z_axis':
          rot = float(this[i+1])
          matrix2,offset2 = self.load_z_axis_rotation_matrix4(rot)
        elif n == 'rotate_about_arbitrary_axis':
          axis = np.array(this[i+1:i+1+3]).astype(float)
          rot = float(this[i+3])
          fix_point = np.array(this[i+4:i+4+3]).astype(float)
          matrix2,offset2 = self.rotate_about_arbitrary_axis(axis,rot,fix_point)
        elif n == 'scale_fix_point':
          scale = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_scaling_fix_point_matrix4(scale,fix_point)
        elif n == 'scale_differential_fix_point':
          scale = np.array(this[i+1:i+1+3]).astype(float)
          fix_point = np.array(this[i+3:i+3+3]).astype(float)
          matrix2,offset2 = self.load_differential_scaling_fix_point_matrix4(\
                                    scale,fix_point)
        elif n == 'rotate_about_x_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_x_axis_rotation_fix_point_matrix4(\
                                     theta,fix_point)
        elif n == 'rotate_about_y_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_y_axis_rotation_fix_point_matrix4(\
                                     theta,fix_point)
        elif n == 'rotate_about_z_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_z_axis_rotation_fix_point_matrix4(\
                                     theta,fix_point)
        else:
          noGood = True
        if not noGood:
          clone.matrix = matrix2*matrix
          clone.matrixI = clone.matrix.I
          
      if (matrix != np.matrix(np.eye(3))).all():
        clone.matrix = matrix
      if np.abs(offset).sum() > 0:
        clone.offset = offset
      clone.contents = []

      self.top.contents.append(clone)
      self.verbose == 2 and sys.stderr.write('<M>')
    except:
      self.error('could not interpret clone line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))


  def load_z_axis_rotation_fix_point_matrix4(self,theta,fix_point):
    '''
    rotation of theta in z axis about fix_point
    '''
    

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
           PyRatEllipsoid(self.point[this],\
             np.array([cmd[2],cmd[3],cmd[4]]).astype(float),\
             material=self.top.material))
      self.verbose == 2 and sys.stderr.write('e')
    except:
      self.error('could not interpret "ell" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def cyl(self,cmd,info={}):
    '''
    Cylinder
  
    cyl base tip radius
    '''
    from PyRatCylinder import PyRatCylinder
    try:
      this = np.array([cmd[1],cmd[2]]).astype(int)
      this[this<0] += self.nPoints
      info.update({'radius':float(cmd[3])})
      self.top.contents.append(\
           PyRatCylinder(self.point[this[0]],self.point[this[1]],info=info,\
                         material=self.top.material))
      self.verbose == 2 and sys.stderr.write('c')
    except:
      self.error('could not interpret "cyl" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def ccyl(self,cmd):
    '''
    Capped (closed) Cylinder
 
    ccyl base tip radius
    '''
    self.cyl(cmd,info={'caps':True})
    self.verbose == 2 and sys.stderr.write('C')


  def sph(self,cmd,N=11):
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
           PyRatSpheroid(self.point[int(this)],float(cmd[2]),\
                         material=self.top.material))
      self.verbose == 2 and sys.stderr.write('s')

      #if self.GL and hasGL:
      #  this = self.top.contents[-1]
      #  this.GLRepresentation = this.tesselate(N=N)

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
      self.infinitePlane.append(\
           PyRatPlane(self.point[this[1]],self.point[this[0]],\
                      material=self.top.material))
      self.verbose == 2 and sys.stderr.write('<plane>')
    except:
      self.error('could not interpret "plane" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def usemtl(self,cmd):
    '''
    Set the current material
    '''
    try:
      self.top.material = cmd[1]
      self.verbose == 2 and sys.stderr.write('(%s)'%cmd[1])
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
      reportingFrequency = np.max([1,int(l/self.reportingFrequency)])
      for i,line in enumerate(this):
        if self.verbose and i%reportingFrequency == 0:
          sys.stderr.write('\b\b\b\b\b\b\b\b%.2f%%'%(100*i/l))
        if len(line):
          self.lineNumber = i+1
          self.parseLine(line.split())
      try:
        self.closeBox('')
      except:
        pass
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

    '''
    from PyRatBox import PyRatBox
    box = PyRatBox(np.zeros(3),None)
    box.invisible = True
    box.empty = True
    self.top.contents.append(box)
    self.stack.append(self.top)
    self.top = self.top.contents[-1]

    self.verbose == 2 and sys.stderr.write('[%d}'%len(self.stack))

  def box(self,cmd):
    '''
    Open a solid box

    Start a new container box and push on stack
   
    box minx miny minz maxx maxy maxz 
    '''
    from PyRatBox import PyRatBox
    box = PyRatBox(np.array(cmd[1:1+3]).astype(float),np.array(cmd[4:4+3]).astype(float))
    box.invisible = False
    box.empty = False
    self.top.contents.append(box)
    
    self.verbose == 2 and sys.stderr.write('B'%len(self.stack))

  def closeBox(self,cmd):
    '''
    Close a bounding box
    '''
    if len(self.top.contents) == 0:
      self.top.empty = True
    else:
      self.top.empty = False

    #if len(self.top.contents) == 1:
    #  # only one item, so put it in the upper box
    #  if len(self.stack[-1].contents) == 1 and self.stack[-1].contents[0] == self.top:
    #    # replace item above with this one
    #    older = self.stack[-1].contents[0]
    #    self.stack[-1].contents[0] = self.top.contents[0]
    #    older.contents = []
    #    older.empty = True
    #    del older
    #  else: 
    #    self.stack[-1].contents.append(self.top.contents[0])
    if len(self.top.contents) > 0:
       for c,i in enumerate(self.top.contents):
         if type(i) == PyRatBox and (i.isDefined() or i.empty):
           self.top.contents[c] = None  
       self.top.contents = [j for j in np.array(self.top.contents)[np.array([i != None for i in self.top.contents])]]
    try:
      self.top = self.stack.pop()
      self.verbose == 2 and sys.stderr.write('{%d]'%len(self.stack))
    except:
      pass



  def mtllib(self,cmd):
    '''
    Material library
    '''
    self.root.materialFile = cmd[1] 
    self.verbose == 2 and sys.stderr.write('(Materials: %s)'%cmd[1])

  def g(self,cmd):
    '''
    Group:

      g group name
    '''
    try:
      groupName = ' '.join(cmd[1:])
      self.verbose == 2 and sys.stderr.write('(g %s)'%groupName)
    except:
      self.error('could not interpret "g" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))
    # now we need to reference for this object
    # This should be the last item on the stack but if not its the root
    self.group[groupName] = self.top

  def v(self,cmd):
    '''
    Vertex:

      v x y z
    '''
    try:
      self.point.append(np.array(cmd[1:4]).astype(float))
      self.nPoints += 1
      self.verbose == 2 and sys.stderr.write('v')
    except:
      self.error('could not interpret "v" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def disk(self,cmd,N=8):
    '''
    Disk
    
      disk centre normal radius
    '''
    from PyRatDisk import PyRatDisk
    isOK = True
    try:
      this = np.array(cmd[1:3]).astype(int)
      this[this<0] += self.nPoints
      radius = float(cmd[3])
      self.top.contents.append(\
           PyRatDisk(self.point[this[0]],self.point[this[1]],\
                     material=self.top.material,info={'radius':radius}))
      self.verbose == 2 and sys.stderr.write('d')
    except:
      isOK = False
      self.error('could not interpret "disk" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))
    
    #if self.GL and isOK:
    #  this = self.top.contents[-1]
    #  this.GLRepresentation = this.tesselate(N=N)

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
                      self.point[this[2]]]),None,\
                      material=self.top.material))
      self.verbose == 2 and sys.stderr.write('f')
    except:
      self.error('could not interpret "f" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

  def write(self,filename,type='compatible'):
    '''
    Write out a wavefront format representation
    '''
    M = np.eye(3)
    C = np.matrix(np.zeros(3))
    file = open(filename,'w')
    try:
      self.root.write(file,M=M.copy(),C=C.copy(),type=type)      
    except:
      pass


def main():
  '''
  Test to demonstrate reading an obj file
  and using clones

  '''
  from PyRatObjParser import PyRatObjParser
  from PyRatClone import PyRatClone
  from PyRatBox import test
  filename = 'tests/clone2.obj'
  filename = 'tests/new_plant.obj'
  if hasGL:
    tvtk.Property(representation='wireframe') 
  world = PyRatObjParser(filename,verbose=True,GL=True)
  if world.GL:
    if len(world.triangleData):
      triangleData = np.array(world.triangleData)
      triangleN = np.array(world.triangleN)
      mesh = tvtk.PolyData()
      mesh.points = triangleData
      mesh.polys = triangleN
      # add mesh.point_data.scalars = temperature
      # at some point
    mlab.show()
  if world.root.size == 0:
    world.error('Zero size in world root')
    return False
  world.root.planes = world.infinitePlane

  info = {'verbose':True}
  name = str(globals()['__file__'].split('/')[-1].split('.')[0])
  test(np.zeros(3),np.zeros(3),obj=world.root,info=info,type=name,nAtTime=100*100/20,name=name[5:])
  return True  

if __name__ == "__main__":
    main()

