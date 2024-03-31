# Example guide on using Mayavi for 3D visualisations
# Example 1: Defect schematic
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

def exampledefect():
	# Example on how to make a nematic defect

	# Defect charge
	charge = 0.5

	# System size
	xyz = [5,5,1]

	# Make director field
	dir = np.zeros(shape=(3,xyz[0],xyz[1],xyz[2]))			# orientation
	pos = np.zeros(shape=(3,xyz[0],xyz[1],xyz[2]))			# position (don't need)
	for x in range(xyz[0]):
		for y in range(xyz[1]):
			phi = math.atan2((y-0.5*xyz[1]+0.5),(x-0.5*xyz[0]+0.5))		# so phi defined with defect core at origin

			# Remove director at defect core
			if x == int(0.5*xyz[0]) and y == int(0.5*xyz[1]):
				dir[0][x][y][0] = 0.
				dir[1][x][y][0] = 0.

			# Nematic orientation
			else:
				dir[0][x][y][0] = math.cos(charge*phi)
				dir[1][x][y][0] = math.sin(charge*phi)

			# Position field
			pos[0][x][y][0] = x
			pos[1][x][y][0] = y
			pos[2][x][y][0] = 0.

	# Plot director field
	directorcolour = (194./255., 211./255., 223./255.)
	DF = mlab.quiver3d(pos[0],pos[1],pos[2],dir[0],dir[1],dir[2],color=directorcolour,scale_factor = 0.7,mode='cylinder',opacity=0.8)
	# DF = mlab.quiver3d(dir[0],dir[1],dir[2],color=directorcolour,scale_factor = 0.7,mode='cylinder',opacity=0.8)
	DF.glyph.glyph_source.glyph_source.center = np.array([ 0. , 0.,  0. ])			# remove shift
	DF.actor.property.specular = 0.2
	DF.actor.property.specular_power = 10
	DF.actor.property.lighting = True
	DF.actor.property.edge_visibility = 0
	DF.actor.property.backface_culling = True

	# Plot defect core
	corecolour = (0./255.,50./255.,95./255)
	spheres = mlab.points3d([int(0.5*xyz[0])],[int(0.5*xyz[1])],[0.],[0.5],resolution=60,scale_factor=1,color=corecolour,opacity=0.8)
	spheres.glyph.glyph_source.glyph_source.center = np.array([ 0. , 0.,  0. ])		# remove shift

################################################################################
saveshow = 'show'
savename = 'mayaviexampledefect.png'

# Setting up the Mayavi scene
if saveshow == 'save':
	# Needs to be before mlab.figure
	# Stops the pipeline from appearing.
	# Useful if you want to save multiple frames from a simulation
	mlab.options.offscreen = True

# Setting up the Mayavi scene
mlab.figure(size=(1024,768),bgcolor=(1.,1.,1.),fgcolor=(160./255.,160./255.,160./255.))

# Polymer example
exampledefect()

# Save
if saveshow == 'save':
	mlab.savefig(savename)

if saveshow == 'show':
	# Show the scene (can interact with pipeline)
	mlab.show()
