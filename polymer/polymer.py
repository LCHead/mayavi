# Example guide on using Mayavi for 3D visualisations
# Example 1: Polymer
# By Louise Head
# 1st April 2024

################################################################################
# Import libraries
import numpy as np
import math
from mayavi import mlab
from pprint import pprint

################################################################################
# Routines

def examplepolymer():
	# Example on how to make a polymer

	# System size (not required)
	xyz = [10,10,10]

	# Points that make up the polymer
	xpoints = [0.28*xyz[0],0.4*xyz[0],0.6*xyz[0],0.67*xyz[0],0.57*xyz[0],0.4*xyz[0]]
	ypoints = [0.86*xyz[1],0.68*xyz[1],0.65*xyz[1],0.38*xyz[1],0.21*xyz[1],0.1*xyz[1]]
	zpoints = [0.1*xyz[2], 0.4*xyz[2], 0.45*xyz[2], 0.65*xyz[2], 0.8*xyz[2], 0.9*xyz[2]]
	diameter = [1.,1.,1.,1.,1.,1.]

	# Colour and opacity
	polycolour = (157./255., 217./255., 243./255.)
	polyopacity = 0.6
	bondcolour = (157./255., 217./255., 243./255.)
	bondopacity = 0.3

	# Visualise monomers
	spheres = mlab.points3d(xpoints,ypoints,zpoints,diameter,resolution=60,scale_factor=1,color=polycolour,opacity=polyopacity)

	# Set visualisation preferences
	spheres.actor.property.ambient = 0.25
	spheres.actor.property.diffuse = 0.4
	spheres.actor.property.specular = 0.35
	spheres.actor.property.specular_power = 3
	spheres.actor.property.backface_culling = True
	spheres.actor.property.lighting = True

	# Visualise bonds
	bonds = mlab.plot3d(xpoints,ypoints,zpoints,color=bondcolour,opacity=bondopacity,tube_radius = 0.2)

	# Visualise box
	box = mlab.outline(extent=[0,xyz[0],0,xyz[1],0,xyz[2]],line_width=2)

################################################################################
saveshow = 'show'
savename = 'mayaviexamplepolymer.png'

# Setting up the Mayavi scene
if saveshow == 'save':
	# Needs to be before mlab.figure
	# Stops the pipeline from appearing.
	# Useful if you want to save multiple frames from a simulation
	mlab.options.offscreen = True

# Setting up the Mayavi scene
mlab.figure(size=(1024,768),bgcolor=(1.,1.,1.),fgcolor=(160./255.,160./255.,160./255.))

# Polymer example
examplepolymer()

# Save
if saveshow == 'save':
	mlab.savefig(savename)

if saveshow == 'show':
	# Show the scene (can interact with pipeline)
	mlab.show()
