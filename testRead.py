#!/usr/bin/env python

'''
A test PyRay scripy

Here, we read data from an extended wavefront format file (file)
to form a world object. The root of the contents of this 
is in world.root and we could use this directly to do ray tracing etc.

We can get a report on the structure with:

world.root.report()

but here, instead, we manually warp a clone object
around the world

where we have shrunk by scaling by 0.025
and rotated 15 degrees.

When we run this, we can specify the number of processors to use

e.g. 

testRead.py 10

Specifying 0 means that it is done on a single processor.
Not specifying anything, it will use all core on the
current processor.

There is an overhead to running the parallel jobs, but it will generally be faster.

The granularity can be changed by using nAtTime in calling test() which specifies how many
primary rays to pass to each parallel job. You'd generally want that to be around the number
of total primary rays divided by the number of processors (or perhaps half of that)

Process time is around 1.5 minutes on a machine with e.g. 12 cores.

The result in in tests/PyRat-RAMI-near.png

'''
import PyRat
from PyRat import test, PyRatObjParser, PyRatClone
import numpy as np
import time
import pdb
file = 'PyRat/spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
world = PyRatObjParser(file,verbose=True)
clone = PyRatClone(np.zeros(3),None)
clone.thisGroup = None
clone.offset = np.array([0,0,0.])
clone.matrix = np.eye(3)

c = np.cos(-15*np.pi/180.)
s = np.sin(-15*np.pi/180.)
clone.matrix[1,1] = clone.matrix[0,0] = c
clone.matrix[0,1] = -s
clone.matrix[1,0] = s
clone.matrix *= 0.025
clone.empty = False
clone.invisible = True
clone.contents = [world.root]
world.reconcile(clone,0)

clone.updateBbox()
info = {'verbose':True}

t0 = time.clock()
t0w = time.time()
test(np.zeros(3),np.zeros(3),obj=clone,info=info,type=None,file='RAMI',nAtTime=100*100/20)
print time.clock() - t0, "seconds process time"
print time.time() - t0w, "seconds wall time"
