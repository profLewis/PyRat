from PyRatBox import *

import pyglet
from pyglet.window import mouse
from pyglet.window import key
from pyglet.gl import *


platform = pyglet.window.get_platform()
display = platform.get_default_display()
screen = display.get_default_screen()

template = pyglet.gl.Config(alpha_size=8,sample_buffers=1, samples=4)
try:
  config = screen.get_best_config(template)
except pyglet.window.NoSuchConfigException:
    template = gl.Config()
    config = screen.get_best_config(template)

context = config.create_context(None)
window = pyglet.window.Window(caption='PyRat',context=context,config=config)

@window.event
def on_key_press(symbol, modifiers):
    if symbol == key.A:
        print 'The "A" key was pressed.'
    elif symbol == key.Q:
        pyglet.app.exit()
    elif symbol == key.LEFT:
        print 'The left arrow key was pressed.'
    elif symbol == key.ENTER:
        print 'The enter key was pressed.'
    else:
        print 'A key was pressed'

@window.event
def on_mouse_press(x, y, button, modifiers):
    if button == mouse.LEFT:
        print 'The left mouse button was pressed.'

vertices = [
    0, 0,
    window.width, 0,
    window.width, window.height]
vertices_gl = (GLfloat * len(vertices))(*vertices)

glEnableClientState(GL_VERTEX_ARRAY)
glVertexPointer(2, GL_FLOAT, 0, vertices_gl)

@window.event
def on_resize(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(65, width / float(height), .1, 1000)
    glMatrixMode(GL_MODELVIEW)
    return pyglet.event.EVENT_HANDLED

@window.event
def on_draw():
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    glDrawArrays(GL_TRIANGLES, 0, len(vertices) // 2)


pyglet.app.run()



def main():
  '''
  A simple test of the bounding box algorithm
  
  A scan over a cube is made and 2 images produced
  tests/PyRatBox-near.png and tests/PyRatBox-far.png
  with the near and far distances 
  '''
  import sys
  import os

  from PyRatRay import PyRatRay
  import pylab as plt

  # set up a test object: a cube
  min = [-0.5,-0.5,0]
  extent = [1,1,1]

  box = PyRatBox(min,extent)

  # ray direction
  direction = np.array([0,0,-1])

  # image size
  size = (100,100)

  # ray origins
  origin = np.array([0,0,2]).astype(float)
  # dimensions of the image in physical units
  dimensions = [4,4]

  result1 = np.zeros(size) 
  result2 = np.zeros(size)

  o = origin.copy()
  ray = PyRatRay(o,direction)
  for ix in xrange(size[0]):
    o[0] = origin[0] + dimensions[0] * (ix-size[0]*0.5)/size[0]   
    for iy in xrange(size[1]):
      o[1] = origin[1] + dimensions[1] * (iy-size[1]*0.5)/size[1] 
      if box.intersect(ray):
        distance = ray.tnear * ray.direction
        result1[ix,iy] = np.sqrt(np.dot(distance,distance))
        distance = ray.tfar * ray.direction
        result2[ix,iy] = np.sqrt(np.dot(distance,distance))

  plt.imshow(result1,interpolation='nearest')
  plt.colorbar()
  if not os.path.exists('tests'):
    os.makedirs('tests') 
  plt.savefig('tests/PyRatBox-near.png')
  plt.clf()
  plt.imshow(result2,interpolation='nearest')
  plt.colorbar()
  plt.savefig('tests/PyRatBox-far.png')


if __name__ == "__main__":
    pyglet.app.run()
