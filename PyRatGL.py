#! /usr/bin/env python

from OpenGLContext import testingcontext
BaseContext = testingcontext.getInteractive()

import OpenGL
OpenGL.ERROR_CHECKING = False
OpenGL.ERROR_LOGGING = False
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGLContext.arrays import *
from OpenGL.GL import shaders


class TestContext( BaseContext ):
  """Creates a simple vertex shader."""
  def OnInit( self ):
    VERTEX_SHADER = shaders.compileShader("""#version 330
        void main() {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
        }""", GL_VERTEX_SHADER)
    FRAGMENT_SHADER = shaders.compileShader("""#version 330
        void main() {
            gl_FragColor = vec4( 0, 1, 0, 1 );
        }""", GL_FRAGMENT_SHADER)
    self.shader = shaders.compileProgram(VERTEX_SHADER,FRAGMENT_SHADER)

    self.vbo = vbo.VBO(
            array( [
                [  0, 1, 0 ],
                [ -1,-1, 0 ],
                [  1,-1, 0 ],
                [  2,-1, 0 ],
                [  4,-1, 0 ],
                [  4, 1, 0 ],
                [  2,-1, 0 ],
                [  4, 1, 0 ],
                [  2, 1, 0 ],
            ],'f')
        )

  def Render( self, mode):
    """Render the geometry for the scene."""
    shaders.glUseProgram(self.shader)
    try:
        self.vbo.bind()
        try:
          glEnableClientState(GL_VERTEX_ARRAY);
          glVertexPointerf( self.vbo )
          glDrawArrays(GL_TRIANGLES, 0, 9)
        finally:
          self.vbo.unbind()
          glDisableClientState(GL_VERTEX_ARRAY);
    finally:
      shaders.glUseProgram( 0 )

def test():
  from PyRatClone import PyRatClone
  from PyRatBox import test

  filename = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
  filename = 'tests/clone3.obj'
  world = PyRatObjParser(filename,verbose=True,GL=True)
  import pdb;pdb.set_trace()
  print 'done'

if __name__ == "__main__":
    TestContext.ContextMainLoop()

  
