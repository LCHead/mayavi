# Example guide on using Mayavi for 3D visualisations
# Example 3: Vortex lattice
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

def readVTKvfield(f,xyz):
	# Reading a vector field in from vtk format

	# Initialise field
	field = np.zeros(shape=(3,xyz[0],xyz[1],xyz[2]),dtype=float)

	for z in range(xyz[2]):
		for y in range(xyz[1]):
			for x in range(xyz[0]):

				# read line
				l=f.readline().split()

				# Stop the iteration through timesteps if at end
				if len(l)!=3:
					infile = False
					break

				# Find field
				else:
					vx,vy,vz = [float(l[0]),float(l[1]),float(l[2])]
					field[0][x][y][z] = vx
					field[1][x][y][z] = vy
					field[2][x][y][z] = vz

	return field
def readVTKsfield(f,xyz):
	# Reading a scalar field in from vtk format

	# Initialise field
	field = np.zeros(shape=(xyz[0],xyz[1],xyz[2]),dtype=float)

	for z in range(xyz[2]):
		for y in range(xyz[1]):
			for x in range(xyz[0]):

				# read line
				l=f.readline().split('\t')

				# Stop the iteration through timesteps if at end
				if len(l)!=1:
					infile = False
					break

				# Find field
				else:
					s = float(l[0])
					field[x][y][z] = s
	return field
def tossHeader(f, n):
	for i in range(n):
		line = f.readline()
def systemSizeVTK(filename):
	# Finding the system size from a .vtk file

	# Initialise xyz array
	xyz=np.zeros(shape=(3),dtype=int)

	# Extract dimensions from header of an input file
	f=open(filename,'r')
	_=tossHeader(f,4)
	l=f.readline().split()
	for i in range(3):
		xyz[i]=int(l[i+1])
	f.close()

	return xyz

def calcQtensor(xyz,dir,s):
	# Calculating the nematic Q tensor

	# Initialise Q tensor
	Q = np.zeros(shape=(3, 3, xyz[0], xyz[1], xyz[2]),dtype=float)

	for x in range(xyz[0]):
		for y in range(xyz[1]):
			for z in range(xyz[2]):

				# Iterate rows (i) and columns (j)
				for i in range(3):
					for j in range(3):
						Q[i][j][x][y][z]=3.0*dir[i][x][y][z]*dir[j][x][y][z]

						# Diagonal elements
						if(i==j):
							Q[i][j][x][y][z]-=1.0

						Q[i][j][x][y][z]*=0.5*s[x][y][z]

	return Q
def calcD(xyz,Q,walltypes):
	# Calculating disclination density tensor

	# Initialise array
	DT = np.zeros(shape=(3, 3, xyz[0],xyz[1],xyz[2]))

	for x in range(xyz[0]):
		for y in range(xyz[1]):
			for z in range(xyz[2]):

				# Neighbouring cells
				U = y+1			# up
				D = y-1			# down
				R = x+1			# right
				L = x-1			# left
				T = z+1			# top
				B = z-1			# bottom

				# Applying PBC
				if walltypes[0] == 'p':
					if R >= xyz[0]:
						R -= xyz[0]
					if L < 0:
						L += xyz[0]
				if walltypes[1] == 'p':
					if U >= xyz[1]:
						U -= xyz[1]
					if D < 0:
						D += xyz[1]
				if walltypes[2] == 'p':
					if T >= xyz[2]:
						T -= xyz[2]
					if B < 0:
						B += xyz[2]

				# Calculate derivatives
				if walltypes[0] == 'h':
					# x- derivatives
					if (x == 0):
						dxQxx=( Q[0][0][R][y][z]-Q[0][0][x][y][z] )
						dxQyy=( Q[1][1][R][y][z]-Q[1][1][x][y][z] )
						dxQzz=( Q[2][2][R][y][z]-Q[2][2][x][y][z] )
						dxQxy=( Q[0][1][R][y][z]-Q[0][1][x][y][z] )
						dxQyx=dxQxy
						dxQyz=( Q[1][2][R][y][z]-Q[1][2][x][y][z] )
						dxQzy=dxQyz
						dxQxz=( Q[0][2][R][y][z]-Q[0][2][x][y][z] )
						dxQzx=dxQxz
					elif (x == xyz[0]-1):
						dxQxx=( Q[0][0][x][y][z]-Q[0][0][L][y][z] )
						dxQyy=( Q[1][1][x][y][z]-Q[1][1][L][y][z] )
						dxQzz=( Q[2][2][x][y][z]-Q[2][2][L][y][z] )
						dxQxy=( Q[0][1][x][y][z]-Q[0][1][L][y][z] )
						dxQyx=dxQxy
						dxQyz=( Q[1][2][x][y][z]-Q[1][2][L][y][z] )
						dxQzy=dxQyz
						dxQxz=( Q[0][2][x][y][z]-Q[0][2][L][y][z] )
						dxQzx=dxQxz
					else:
						dxQxx=0.5*( Q[0][0][R][y][z]-Q[0][0][L][y][z] )
						dxQyy=0.5*( Q[1][1][R][y][z]-Q[1][1][L][y][z] )
						dxQzz=0.5*( Q[2][2][R][y][z]-Q[2][2][L][y][z] )
						dxQxy=0.5*( Q[0][1][R][y][z]-Q[0][1][L][y][z] )
						dxQyx=dxQxy
						dxQyz=0.5*( Q[1][2][R][y][z]-Q[1][2][L][y][z] )
						dxQzy=dxQyz
						dxQxz=0.5*( Q[0][2][R][y][z]-Q[0][2][L][y][z] )
						dxQzx=dxQxz
				elif walltypes[0] == 'p':
					dxQxx=0.5*( Q[0][0][R][y][z]-Q[0][0][L][y][z] )
					dxQyy=0.5*( Q[1][1][R][y][z]-Q[1][1][L][y][z] )
					dxQzz=0.5*( Q[2][2][R][y][z]-Q[2][2][L][y][z] )
					dxQxy=0.5*( Q[0][1][R][y][z]-Q[0][1][L][y][z] )
					dxQyx=dxQxy
					dxQyz=0.5*( Q[1][2][R][y][z]-Q[1][2][L][y][z] )
					dxQzy=dxQyz
					dxQxz=0.5*( Q[0][2][R][y][z]-Q[0][2][L][y][z] )
					dxQzx=dxQxz

				if walltypes[1] == 'h':
					# y- derivatives
					if (y == 0):
						dyQxx=( Q[0][0][x][U][z]-Q[0][0][x][y][z] )
						dyQyy=( Q[1][1][x][U][z]-Q[1][1][x][y][z] )
						dyQzz=( Q[2][2][x][U][z]-Q[2][2][x][y][z] )
						dyQxy=( Q[0][1][x][U][z]-Q[0][1][x][y][z] )
						dyQyx=dyQxy
						dyQyz=( Q[1][2][x][U][z]-Q[1][2][x][y][z] )
						dyQzy=dyQyz
						dyQxz=( Q[0][2][x][U][z]-Q[0][2][x][y][z] )
						dyQzx=dyQxz
					elif (y == xyz[1]-1):
						dyQxx=( Q[0][0][x][y][z]-Q[0][0][x][D][z] )
						dyQyy=( Q[1][1][x][y][z]-Q[1][1][x][D][z] )
						dyQzz=( Q[2][2][x][y][z]-Q[2][2][x][D][z] )
						dyQxy=( Q[0][1][x][y][z]-Q[0][1][x][D][z] )
						dyQyx=dyQxy
						dyQyz=( Q[1][2][x][y][z]-Q[1][2][x][D][z] )
						dyQzy=dyQyz
						dyQxz=( Q[0][2][x][y][z]-Q[0][2][x][D][z] )
						dyQzx=dyQxz
					else:
						dyQxx=0.5*( Q[0][0][x][U][z]-Q[0][0][x][D][z] )
						dyQyy=0.5*( Q[1][1][x][U][z]-Q[1][1][x][D][z] )
						dyQzz=0.5*( Q[2][2][x][U][z]-Q[2][2][x][D][z] )
						dyQxy=0.5*( Q[0][1][x][U][z]-Q[0][1][x][D][z] )
						dyQyx=dyQxy
						dyQyz=0.5*( Q[1][2][x][U][z]-Q[1][2][x][D][z] )
						dyQzy=dyQyz
						dyQxz=0.5*( Q[0][2][x][U][z]-Q[0][2][x][D][z] )
						dyQzx=dyQxz
				elif walltypes[1] == 'p':
						dyQxx=0.5*( Q[0][0][x][U][z]-Q[0][0][x][D][z] )
						dyQyy=0.5*( Q[1][1][x][U][z]-Q[1][1][x][D][z] )
						dyQzz=0.5*( Q[2][2][x][U][z]-Q[2][2][x][D][z] )
						dyQxy=0.5*( Q[0][1][x][U][z]-Q[0][1][x][D][z] )
						dyQyx=dyQxy
						dyQyz=0.5*( Q[1][2][x][U][z]-Q[1][2][x][D][z] )
						dyQzy=dyQyz
						dyQxz=0.5*( Q[0][2][x][U][z]-Q[0][2][x][D][z] )
						dyQzx=dyQxz

				if walltypes[2] == 'h':
					# z- derivatives
					if (z == 0):
						dzQxx=( Q[0][0][x][y][T]-Q[0][0][x][y][z] )
						dzQyy=( Q[1][1][x][y][T]-Q[1][1][x][y][z] )
						dzQzz=( Q[2][2][x][y][T]-Q[2][2][x][y][z] )
						dzQxy=( Q[0][1][x][y][T]-Q[0][1][x][y][z] )
						dzQyx=dzQxy
						dzQyz=( Q[1][2][x][y][T]-Q[1][2][x][y][z] )
						dzQzy=dzQyz
						dzQxz=( Q[0][2][x][y][T]-Q[0][2][x][y][z] )
						dzQzx=dzQxz
					elif (z == xyz[2] - 1):
						dzQxx=( Q[0][0][x][y][z]-Q[0][0][x][y][B] )
						dzQyy=( Q[1][1][x][y][z]-Q[1][1][x][y][B] )
						dzQzz=( Q[2][2][x][y][z]-Q[2][2][x][y][B] )
						dzQxy=( Q[0][1][x][y][z]-Q[0][1][x][y][B] )
						dzQyx=dzQxy
						dzQyz=( Q[1][2][x][y][z]-Q[1][2][x][y][B] )
						dzQzy=dzQyz
						dzQxz=( Q[0][2][x][y][z]-Q[0][2][x][y][B] )
						dzQzx=dzQxz
					else:
						dzQxx=0.5*( Q[0][0][x][y][T]-Q[0][0][x][y][B] )
						dzQyy=0.5*( Q[1][1][x][y][T]-Q[1][1][x][y][B] )
						dzQzz=0.5*( Q[2][2][x][y][T]-Q[2][2][x][y][B] )
						dzQxy=0.5*( Q[0][1][x][y][T]-Q[0][1][x][y][B] )
						dzQyx=dzQxy
						dzQyz=0.5*( Q[1][2][x][y][T]-Q[1][2][x][y][B] )
						dzQzy=dzQyz
						dzQxz=0.5*( Q[0][2][x][y][T]-Q[0][2][x][y][B] )
						dzQzx=dzQxz
				elif walltypes[2] == 'p':
					dzQxx=0.5*( Q[0][0][x][y][T]-Q[0][0][x][y][B] )
					dzQyy=0.5*( Q[1][1][x][y][T]-Q[1][1][x][y][B] )
					dzQzz=0.5*( Q[2][2][x][y][T]-Q[2][2][x][y][B] )
					dzQxy=0.5*( Q[0][1][x][y][T]-Q[0][1][x][y][B] )
					dzQyx=dzQxy
					dzQyz=0.5*( Q[1][2][x][y][T]-Q[1][2][x][y][B] )
					dzQzy=dzQyz
					dzQxz=0.5*( Q[0][2][x][y][T]-Q[0][2][x][y][B] )
					dzQzx=dzQxz

				DT[0][0][x][y][z] = dyQyx*dzQzx + dyQyy*dzQzy + dyQyz*dzQzz - dyQzx*dzQyx - dyQzy*dzQyy - dyQzz*dzQyz + dzQzx*dyQyx + dzQzy*dyQyy + dzQzz*dyQyz - dzQyx*dyQzx - dzQyy*dyQzy - dzQyz*dyQzz
				DT[0][1][x][y][z] = -dxQyx*dzQzx - dxQyy*dzQzy - dxQyz*dzQzz + dxQzx*dzQyx + dxQzy*dzQyy + dxQzz*dzQyz - dzQzx*dxQyx - dzQzy*dxQyy - dzQzz*dxQyz + dzQyx*dxQzx + dzQyy*dxQzy + dzQyz*dxQzz
				DT[0][2][x][y][z] = dxQyx*dyQzx + dxQyy*dyQzy + dxQyz*dyQzz - dxQzx*dyQyx - dxQzy*dyQyy - dxQzz*dyQyz + dyQzx*dxQyx + dyQzy*dxQyy + dyQzz*dxQyz - dyQyx*dxQzx - dyQyy*dxQzy - dyQyz*dxQzz

				DT[1][0][x][y][z] = -dyQxx*dzQzx - dyQxy*dzQzy - dyQxz*dzQzz + dyQzx*dzQxx + dyQzy*dzQxy + dyQzz*dzQxz - dzQzx*dyQxx - dzQzy*dyQxy - dzQzz*dyQxz + dzQxx*dyQzx + dzQxy*dyQzy + dzQxz*dyQzz
				DT[1][1][x][y][z] = dxQxx*dzQzx + dxQxy*dzQzy + dxQxz*dzQzz - dxQzx*dzQxx - dxQzy*dzQxy - dxQzz*dzQxz + dzQzx*dxQxx + dzQzy*dxQxy + dzQzz*dxQxz - dzQxx*dxQzx - dzQxy*dxQzy - dzQxz*dxQzz
				DT[1][2][x][y][z] = -dxQxx*dyQzx - dxQxy*dyQzy - dxQxz*dyQzz + dxQzx*dyQxx + dxQzy*dyQxy + dxQzz*dyQxz - dyQzx*dxQxx - dyQzy*dxQxy - dyQzz*dxQxz + dyQxx*dxQzx + dyQxy*dxQzy + dyQxz*dxQzz

				DT[2][0][x][y][z] = dyQxx*dzQyx + dyQxy*dzQyy + dyQxz*dzQyz - dyQyx*dzQxx - dyQyy*dzQxy - dyQyz*dzQxz + dzQyx*dyQxx + dzQyy*dyQxy + dzQyz*dyQxz - dzQxx*dyQyx - dzQxy*dyQyy - dzQxz*dyQyz
				DT[2][1][x][y][z] = -dxQxx*dzQyx - dxQxy*dzQyy - dxQxz*dzQyz + dxQyx*dzQxx + dxQyy*dzQxy + dxQyz*dzQxz - dzQyx*dxQxx - dzQyy*dxQxy - dzQyz*dxQxz + dzQxx*dxQyx + dzQxy*dxQyy + dzQxz*dxQyz
				DT[2][2][x][y][z] = dxQxx*dyQyx + dxQxy*dyQyy + dxQxz*dyQyz - dxQyx*dyQxx - dxQyy*dyQxy - dxQyz*dyQxz + dyQyx*dxQxx + dyQyy*dxQxy + dyQyz*dxQxz - dyQxx*dxQyx - dyQxy*dxQyy - dyQxz*dxQyz

	return DT
def calcFnorm(xyz,D):
	# Calculating the Frobenius norm of tensor field

	norm = np.zeros(shape=(xyz[0],xyz[1],xyz[2]))

	# loop over field
	for x in range(xyz[0]):
		for y in range(xyz[1]):
			for z in range(xyz[2]):

				# initialise norm element for each cell to 0
				nf = 0

				# loop over elements of D and compute sum of squares
				for i in range(3):
					for j in range(3):
						nf += D[i][j][x][y][z]*D[i][j][x][y][z]

				# square root the sum of squares
				norm[x][y][z] = math.sqrt(nf)

	return norm
def calcTr(xyz,T):
	# Calculating the trace

	# Initialise invariant field
	trT = np.zeros(shape=(xyz[0],xyz[1],xyz[2]))

	for x in range(xyz[0]):
		for y in range(xyz[1]):
			for z in range(xyz[2]):

				trT[x][y][z] = T[0][0][x][y][z] + T[1][1][x][y][z] + T[2][2][x][y][z]

	return trT
def dividefields(xyz,A,B):
	# Divide scalar field A by field B

	newfield = np.zeros(shape=(xyz[0],xyz[1],xyz[2]))

	for x in range(xyz[0]):
		for y in range(xyz[1]):
			for z in range(xyz[2]):
				newfield[x][y][z] = A[x][y][z]/B[x][y][z]

	return newfield
def plSmooth(contour):
	# smoothing contour

	PL = mlab.pipeline.poly_data_normals(contour)
	PL.filter.feature_angle = 95
	PL.filter.auto_orient_normals = 1
	PL.filter.compute_cell_normals = 1
	PL.filter.compute_point_normals = 1

	return PL
################################################################################
# Extract data and analyse fields

# Load field data
dirF = open("dir-q-00195000.vtk", "r")		# open director field
scaF = open("op-q-00195000.vtk", "r")		# open scalar order parameter
velF = open("vel-00195000.vtk", "r")		# open velocity field
tossHeader(dirF, 9)				# skip header
tossHeader(scaF, 10)
tossHeader(velF, 9)
xyz = systemSizeVTK("dir-q-00195000.vtk")
dir = readVTKvfield(dirF,xyz)	# read director
sca = readVTKsfield(scaF,xyz)	# read scalar order parameter
vel = readVTKvfield(velF,xyz)	# read velocity field

# Nematic defect analysis
boundarylist = ['p','h','h']	# periodic ('p') or hard wall ('h') (I use this for derivatives)
Q = calcQtensor(xyz,dir,sca)	# nematic tensor order parameter
D = calcD(xyz,Q,boundarylist)	# disclination density tensor
s = calcFnorm(xyz,D)			# scalar field s
trD = calcTr(xyz,D)		# trace of D
cosBeta = dividefields(xyz,trD,s)	# cosBeta from Tr(D)/s

################################################################################
# Visualise

saveshow = 'show'
savename = 'mayaviexamplevortexlattice.png'

if saveshow == 'save':
	# Needs to be before mlab.figure
	# Stops the pipeline from appearing.
	# Useful if you want to save multiple frames from a simulation
	mlab.options.offscreen = True

# Setting up the Mayavi scene
mlab.figure(size=(1024,768),bgcolor=(1.,1.,1.),fgcolor=(160./255.,160./255.,160./255.))

# Visualise velocity field
optionVel = 3
if optionVel == 1:
	# Option 1 -- directly visualise field as a quiver
	VF1 = mlab.quiver3d(vel[0],vel[1],vel[2],mask_points=5,mode='arrow',colormap='plasma',opacity=0.2)
	VF1.glyph.mask_points.maximum_number_of_points = 50000
if optionVel == 2:
	# Load vector field
	VF2 = mlab.pipeline.vector_field(vel[0],vel[1],vel[2])
	# Visualise as vector slice
	VF2 = mlab.pipeline.vector_cut_plane(VF2,mask_points=5,scale_factor=10,mode='arrow',colormap='plasma',opacity=0.4)
	VF2.implicit_plane.widget.normal = [0,1,0]
	VF2.implicit_plane.widget.origin = [0,12.5,0]
	VF2.implicit_plane.visible = False
if optionVel == 3:
	# Option 2 part 2 -- full vector field
	# Load vector field
	VF2 = mlab.pipeline.vector_field(vel[0],vel[1],vel[2])
	# Visualise as glyph
	VF3 = mlab.pipeline.glyph(VF2,mask_points=5,scale_factor=10,mode='arrow',colormap='plasma',opacity=0.3)
	VF3.glyph.mask_points.maximum_number_of_points = 2000
	VF3.glyph.color_mode = 'color_by_vector'
	VF3.glyph.scale_mode = 'scale_by_vector'

# Visualise director field
optionDir = 1
if optionDir == 1:
	DF2 = mlab.pipeline.vector_field(dir[0],dir[1],dir[2])
	DF2 = mlab.pipeline.vector_cut_plane(DF2,mask_points=11,scale_factor=3,mode='2ddash',color=(140./255., 140./255., 140./255.),opacity=0.65)
	DF2.glyph.mask_points.random_mode = False
	DF2.implicit_plane.widget.normal = [0,1,0]
	DF2.implicit_plane.widget.origin = [0,12.5,0]
	DF2.implicit_plane.visible = False

# Visualise defect
optionDefect = 1
if optionDefect == 1:
	dcontour=0.1
	# Create scalar field
	SF = mlab.pipeline.scalar_field(s)
	# Add the colouring field as an additional field
	SF.image_data.point_data.add_array(cosBeta.T.ravel())
	SF.image_data.point_data.get_array(1).name = 'colour'
	SF.update()
	# Make contour of scalar field
	CP = mlab.pipeline.contour(SF)
	CP.filter.contours = [dcontour]
	# Select the other scalar field as active attribute for colouring surface
	CP = mlab.pipeline.set_active_attribute(CP,point_scalars='colour')
	CP = plSmooth(CP)
	# Create surface of the contour coloured by the colouring field
	SP = mlab.pipeline.surface(CP,vmax=1,vmin=-1,colormap='viridis',opacity=1)

# Plot outline
box = mlab.outline(extent=[0,xyz[0],0,xyz[1],0,xyz[2]],line_width=1)

# Save
if saveshow == 'save':
	mlab.savefig(savename)

if saveshow == 'show':
	# Show the scene (can interact with pipeline)
	mlab.show()
