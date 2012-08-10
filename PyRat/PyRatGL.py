#! /usr/bin/env python

'''
'''


import sys
import wx
import numpy
import time
from math import sin, cos, pi, sqrt, acos, asin
from vtk.util import colors as color
# Enthought library imports.
from enthought.traits.api import HasTraits, Trait, Instance, Tuple, Int, \
     Float, Range, Button, Array, Color, Bool, Any, List, Enum, ListFloat, ListInt
from enthought.traits.ui.api import View, Item, Group, RGBColorEditor, RangeEditor
from enthought.traits.ui.message import message
from enthought.tvtk.api import tvtk
from enthought.tvtk.tools import ivtk
from enthought.pyface.api import GUI
from enthought.pyface.timer.api import Timer
from enthought.tvtk.tvtk_base import TVTKBase, vtk_color_trait


from enthought.tvtk.tools.visual import *

class Facet(Sphere):
    """Facet class creates Facet from triangle list, follows the
    usual VTK pipeline and creates a Facet actor, which is passed to
    the show_actor() function as an argument.
    """
    #####################################################################
    # Traits definitions
    coords = List(value = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],0.0, 0.0, 1.0], desc = 'the vertex list')
    indices = List(value = [[0,1,2]],desc = 'the triangle index list')
    color = vtk_color_trait((1.0, 1.0, 1.0))
    representation = Enum('s', 'w', 'p')
    visibility = Bool(True)
    viewer = Any
    polydata = Instance(tvtk.PolyData, ())
    property = Instance(tvtk.Property)
    actor = Instance(tvtk.Actor, ()) # tvtk Actor, for the usual pipeline architecture.
    ######################################################################
    # User interface view     
    traits_view = View(Group(Item(name = 'coords', style = 'simple', label = 'Vertices'),
                             Item(name = 'indices', label = 'Triangle indices'),
                             Item(name = 'color'),
                             Item(name = 'visibility'),
                             Item(name = 'representation'), 
                             label = 'Facet Properties',
                             show_border = True))
    def __init__(self, **traits):
        self.property = self.actor.property
       
        HasTraits.__init__(self, **traits)
        self._create_points(self.coords,self.indices)
        self._color_changed(self.color)
        self._visibility_changed(self.visibility)
        normals = tvtk.PolyDataNormals(input = self.polydata)
        m = tvtk.PolyDataMapper(input = normals.output) # the usual vtk pipleine countinuation
        self.actor.mapper = m
        self.property = self.actor.property
        self.property.representation = self.representation
        show_actor(self.actor) # passing the actors function for rendering
        self.viewer = get_viewer() # getting the ivtk viewer
        self.property.on_trait_change(self.viewer.scene.render)
        self.actor.on_trait_change(self.viewer.scene.render)
    ######################################################################
    # Non-public methods, Event handlers   
    def _create_points(self, c, i):
        self.polydata = tvtk.PolyData(points=numpy.array(c),polys=numpy.array(i))
        self.points = c
        return self.points, self.polydata.polys
       

