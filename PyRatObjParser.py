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
    self.infinitePlane = []

    self.error = self.top.error

    self.group = {}

    self.stack = []
    self.point = []
    self.nPoints = 0
    if self.verbose:
      sys.stderr.write('Reading data from %s\n'%filename)
    self.read(filename)
    if self.verbose:
      sys.stderr.write('\n ... sorting bbox contents\n')
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
      if type(i) == PyRatBox and bbox.contents[c].invisible:
        # remove it
        bbox.contents.pop(c)
    try:
      if not bbox.empty:
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
      if not bbox.empty:
        # if this isnt an empty box then add in its contents
        bbox.min = np.min(np.atleast_2d([bbox.min,min]),axis=0)
        bbox.max = np.max(np.atleast_2d([bbox.max,max]),axis=0)
        bbox.size += size
    except:
      pass
    if bbox.size > 0:
      bbox.empty = False
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
      'transformation_matrix':self.transformation_matrix,\
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
      offset = np.array(this[:3]).astype(float)
      thisGroup = None
      thisGroupName = None
      matrix = np.matrix(np.eye(3))
      done = np.zeros_like(np.array(this)).astype(bool)
      for i in xrange(3,len(this)):
        if not done[i]:
          if this[i] == 'Rx':
            theta = float(this[i+1])
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
            theta = float(this[i+1])
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
              
      clone = PyRatClone(np.ones(3),np.ones(3))
      if (matrix != np.matrix(np.eye(3))).all():
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
        elif n == 'scale_differential_fix_point':
          scale = np.array(this[i+1:i+1+3]).astype(float)
          fix_point = np.array(this[i+3:i+3+3]).astype(float)
          matrix2,offset2 = self.load_differential_scaling_fix_point_matrix4(\
                                    scale,fix_point)
        elif n == 'rotate_about_x_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_x_axis_rotation_fix_point_matrix4(\                                     theta,fix_point)
        elif n == 'rotate_about_y_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_y_axis_rotation_fix_point_matrix4(\
                                     theta,fix_point)
        elif n == 'rotate_about_z_axis_fix_point':
          rot = float(this[i+1])
          fix_point = np.array(this[i+2:i+2+3]).astype(float)
          matrix2,offset2 = self.load_z_axis_rotation_fix_point_matrix4(\                                     theta,fix_point)
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
           PyRatEllipsoid(self.point[this[1]],\
             np.array([this[2],this[3],this[4]]).astype(float),\
             material=self.top.material))
      self.verbose == 2 and sys.stderr.write('e')
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
      self.verbose == 2 and sys.stderr.write('s')
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
    
    self.verbose == 2 and sys.stderr.write('[%d}'%len(self.stack))

  def closeBox(self,cmd):
    '''
    Close a bounding box
    '''
    self.top = self.stack.pop()
    self.verbose == 2 and sys.stderr.write('{%d]'%len(self.stack))
 
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
      self.verbose == 2 and sys.stderr.write('v')
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
      self.verbose == 2 and sys.stderr.write('d')
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
      self.verbose == 2 and sys.stderr.write('f')
    except:
      self.error('could not interpret "f" line %d of %s: %s'%\
                (self.lineNumber,self.filename,' '.join(cmd)))

def main():
  '''
  Test
  '''
  from PyRatObjParser import PyRatObjParser
  filename = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
  world = PyRatObjParser(filename,verbose=True)


if __name__ == "__main__":
    main()

