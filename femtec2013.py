# -*- coding: utf-8 -*-

from lar import *
from scipy import *
import scipy
import numpy as np
from time import time
from pngstack2array3d import pngstack2array3d

DEBUG = False
MAX_CHUNKS = 25

# ------------------------------------------------------------
# configuration parameters
# ------------------------------------------------------------

beginImageStack = 430
colors = 2
endImage = beginImageStack
nx = ny = nz = 50
imageDx = imageDy = imageDz = 50
imageWidth, imageHeight, imageDepth = 250, 250, 250
xEnd = yEnd = zEnd = 0,0,0
count = 0

# ------------------------------------------------------------
# Utility toolbox
# ------------------------------------------------------------

def ind(x,y,z): return x + (nx+1) * (y + (ny+1) * (z))

def invertIndex(nx,ny,nz):
	nx,ny,nz = nx+1,ny+1,nz+1
	def invertIndex0(offset):
		a0, b0 = offset / nx, offset % nx
		a1, b1 = a0 / ny, a0 % ny
		a2, b2 = a1 / nz, a1 % nz
		return b0,b1,b2
	return invertIndex0


# ------------------------------------------------------------
# computation of d-chain generators (d-cells)
# ------------------------------------------------------------

# cubic cell complex
# ------------------------------------------------------------
def the3Dcell(coords):
	x,y,z = coords
	return [ind(x,y,z),ind(x+1,y,z),ind(x,y+1,z),ind(x,y,z+1),ind(x+1,y+1,z),
			ind(x+1,y,z+1),ind(x,y+1,z+1),ind(x+1,y+1,z+1)]

# construction of vertex coordinates (nx x ny x nz)
# ------------------------------------------------------------
V = [[x,y,z] for z in range(nz+1) for y in range(ny+1) for x in range(nx+1) ]


if __name__=="__main__" and DEBUG == True:
	print "\nV =", V

# construction of CV relation (nx x ny x nz)
# ------------------------------------------------------------
CV = [the3Dcell([x,y,z]) for z in range(nz) for y in range(ny) for x in range(nx)]

if __name__=="__main__" and DEBUG == True:
	print "\nCV =", CV
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,CV[:500]+CV[-500:])))
	box = SKELETON(1)(BOX([1,2,3])(hpc))
	#VIEW(STRUCT([box,hpc]))

# construction of FV relation (nx x ny x nz)
# ------------------------------------------------------------
FV = []
v2coords = invertIndex(nx,ny,nz)
for h in range(len(V)):
	x,y,z = v2coords(h)
	if (x < nx) and (y < ny): FV.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])
	if (x < nx) and (z < nz): FV.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])
	if (y < ny) and (z < nz): FV.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])


if __name__=="__main__" and DEBUG == True:
	print "\nFV =",FV
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV[:500]+FV[-500:])))
	box = SKELETON(1)(BOX([1,2,3])(hpc))
	#VIEW(STRUCT([box,hpc]))

# construction of EV relation (nx x ny x nz)
# ------------------------------------------------------------
EV = []
v2coords = invertIndex(nx,ny,nz)
for h in range(len(V)):
	x,y,z = v2coords(h)
	if x < nx: EV.append([h,ind(x+1,y,z)])
	if y < ny: EV.append([h,ind(x,y+1,z)])
	if z < nz: EV.append([h,ind(x,y,z+1)])

if __name__=="__main__" and DEBUG == True:
	print "\nEV =",EV
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EV[:500]+EV[-500:])))
	box = SKELETON(1)(BOX([1,2,3])(hpc))
	#VIEW(STRUCT([box,hpc]))

# ------------------------------------------------------------
# computation of boundary operators (∂3 and ∂2s)
# ------------------------------------------------------------
"""
# computation of the 2D boundary complex of the image space 
# ------------------------------------------------------------
Fx0V, Ex0V = [],[]  # x == 0
Fx1V, Ex1V = [],[]  # x == nx-1
Fy0V, Ey0V = [],[]  # y == 0
Fy1V, Ey1V = [],[]  # y == ny-1
Fz0V, Ez0V = [],[]  # z == 0
Fz1V, Ez1V = [],[]  # z == nz-1
v2coords = invertIndex(nx,ny,nz)
for h in range(len(V)):
	x,y,z = v2coords(h)
	if (z == 0):
		# if x < nx: Ez0V.append([h,ind(x+1,y,z)])
		# if y < ny: Ez0V.append([h,ind(x,y+1,z)])
		if (x < nx) and (y < ny):
			Fz0V.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])
	elif (z == nz):
		# if x < nx: Ez1V.append([h,ind(x+1,y,z)])
		# if y < ny: Ez1V.append([h,ind(x,y+1,z)])
		if (x < nx)  and (y < ny):
			Fz1V.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])

	if (y == 0):
		# if x < nx: Ey0V.append([h,ind(x+1,y,z)])
		# if z < nz: Ey0V.append([h,ind(x,y,z+1)])
		if (x < nx) and (z < nz):
			Fy0V.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])
	elif (y == ny):
		# if x < nx: Ey1V.append([h,ind(x+1,y,z)])
		# if z < nz: Ey1V.append([h,ind(x,y,z+1)])
		if (x < nx) and (z < nz):
			Fy1V.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])

	if (x == 0):
		# if y < ny: Ex0V.append([h,ind(x,y+1,z)])
		# if z < nz: Ex0V.append([h,ind(x,y,z+1)])
		if (y < ny) and (z < nz):
			Fx0V.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])
	elif (x == nx):
		# if y < ny: Ex1V.append([h,ind(x,y+1,z)])
		# if z < nz: Ex1V.append([h,ind(x,y,z+1)])
		if (y < ny) and (z < nz):
			Fx1V.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])

# FbV = Fz0V+Fz1V+Fy0V+Fy1V+Fx0V+Fx1V
# EbV = Ez0V+Ez1V+Ey0V+Ey1V+Ex0V+Ex1V

if __name__=="__main__" and DEBUG == True:
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FbV)))
	VIEW(hpc)
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EbV)))
	VIEW(hpc)

# computation of the ∂2 operator on the boundary space
# ------------------------------------------------------------
# print "start partial_2_b computation"
# partial_2_b = larBoundary(EbV,FbV)
# print "end partial_2_b computation"
"""

# computation of ∂3 operator on the image space
# ------------------------------------------------------------
print "start partial_3 computation"
partial_3 = larBoundary(FV,CV)
print "end partial_3 computation"


# ------------------------------------------------------------
# input from volume image (test: 250 x 250 x 250)
# ------------------------------------------------------------
out = []

for inputIteration in range(imageDepth/imageDz):
	startImage = endImage
	endImage = startImage + imageDz
	xEnd, yEnd = 0,0
	theImage,colors = pngstack2array3d('SLICES2/', startImage, endImage, colors)

	if __name__=="__main__" and DEBUG == True:
		print "\nstartImage, endImage =", (startImage, endImage)
	
	for i in range(imageWidth/imageDx):
		
		for j in range(imageWidth/imageDx):
			
			xStart, yStart = i * imageDx, j * imageDy
			xEnd, yEnd = xStart+imageDx, yStart+imageDy
			
			count += 1
			image = theImage[:, xStart:xEnd, yStart:yEnd]
			nz,nx,ny = image.shape

			if __name__=="__main__" and DEBUG == True:
				print "\n\tsubimage count =",count
				print "\txStart, yStart =", (xStart, yStart)
				print "\txEnd, yEnd =", (xEnd, yEnd)
				print "\timage.shape",image.shape
			
			# ------------------------------------------------------------
			# image elaboration  (chunck: 50 x 50 x 50)
			# ------------------------------------------------------------
			
			
			# compute a quotient complex of chains with constant field
			# ------------------------------------------------------------
			chains3D = [[] for k in range(colors)]
			theColors,k = {},-1
			def addr(x,y,z): return x + (nx) * (y + (ny) * (z))
			for x in range(nx):
				for y in range(ny):
					for z in range(nz):
						if not (image[x,y,z] in theColors):
							k += 1
							theColors[image[x,y,z]] = k
						chains3D[theColors[image[x,y,z]]].append(addr(x,y,z))

			if __name__=="__main__" and DEBUG == True:
				print "\nchains3D =\n", chains3D


			# compute the boundary complex of the quotient cell
			# ------------------------------------------------------------
			objectBoundaryChain = larBoundaryChain(partial_3,chains3D[1])
			b2cells = csrChainToCellList(objectBoundaryChain)
			sup_cell_boundary = MKPOLS((V,[FV[f] for f in b2cells]))
		

			# ------------------------------------------------------------
			# visualize the generated model  
			# ------------------------------------------------------------
			
			zStart = startImage - beginImageStack
			print "xStart, yStart, zStart =", xStart, yStart, zStart

			if __name__=="__main__":
				sup_cell_boundary = MKPOLS((V,[FV[f] for f in b2cells]))
				if sup_cell_boundary != []:
					out += [T([1,2,3])([zStart,xStart,yStart])(STRUCT(sup_cell_boundary))]
				if count == MAX_CHUNKS:
					VIEW(STRUCT(out))

			
			
			
			# ------------------------------------------------------------
			# interrupt the cycle of image elaboration  
			# ------------------------------------------------------------
			if count == MAX_CHUNKS: break
		if count == MAX_CHUNKS: break
	if count == MAX_CHUNKS: break
