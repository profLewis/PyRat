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
import numpy as np
import time
import pdb
file = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_20.obj'
#file = 'spheresTest/HET01_DIS_UNI_NIR_20/HET01_DIS_UNI_NIR_201.obj'
world = PyRat.PyRatObjParser.PyRatObjParser(file,verbose=True)
info = {'verbose':True}

t0 = time.clock()
t0w = time.time()
PyRat.test(np.zeros(3),np.zeros(3),obj=world.root,info=info,type=None,file='RAMI',nAtTime=100*100/50,name='RAMI')
print (time.clock() - t0, "seconds process time")
print (time.time() - t0w, "seconds wall time")
