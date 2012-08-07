

import OpenGL
OpenGL.ERROR_CHECKING = False
OpenGL.ERROR_LOGGING = False
from OpenGL.GL import *
from OpenGL.GLU import *

def main():
  from PyRatObjParser import PyRatObjParser
  from PyRatClone import PyRatClone
  from PyRatBox import test

  filename = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
  filename = 'tests/clone3.obj'
  world = PyRatObjParser(filename,verbose=True,GL=True)
  import pdb;pdb.set_trace()
  print 'done'

if __name__ == "__main__":
    main()


  
