# -*- coding: utf-8 -*-

from lar import *
from scipy import *
import scipy
import numpy as np
from time import time
from pngstack2array3d import pngstack2array3d

colors = 2
theColors = []
DEBUG = False
MAX_CHAINS = colors

# It is VERY important that the below parameter values 
# correspond exactly to each other !!
# ------------------------------------------------------------
MAX_CHUNKS = 75
imageHeight, imageWidth, imageDepth = 250,250,150  # Dx, Dy, Dz


# configuration parameters
# ------------------------------------------------------------

beginImageStack = 430
endImage = beginImageStack
nx = ny = nz = 50
imageDx = imageDy = imageDz = 50
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

def invertPiece(nx,ny,nz):
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

# construction of vertex coordinates (nx * ny * nz)
# ------------------------------------------------------------
V = [[x,y,z] for z in range(nz+1) for y in range(ny+1) for x in range(nx+1) ]


if __name__=="__main__" and DEBUG == True:
	print "\nV =", V

# construction of CV relation (nx * ny * nz)
# ------------------------------------------------------------
CV = [the3Dcell([x,y,z]) for z in range(nz) for y in range(ny) for x in range(nx)]

if __name__=="__main__" and DEBUG == True:
	print "\nCV =", CV
	#hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,CV[:500]+CV[-500:])))
	#box = SKELETON(1)(BOX([1,2,3])(hpc))
	#VIEW(STRUCT([box,hpc]))

# construction of FV relation (nx * ny * nz)
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
	#hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV[:500]+FV[-500:])))
	#box = SKELETON(1)(BOX([1,2,3])(hpc))
	#VIEW(STRUCT([box,hpc]))

# construction of EV relation (nx * ny * nz)
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
	#hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EV[:500]+EV[-500:])))
	#box = SKELETON(1)(BOX([1,2,3])(hpc))
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
		if x < nx: Ez0V.append([h,ind(x+1,y,z)])
		if y < ny: Ez0V.append([h,ind(x,y+1,z)])
		if (x < nx) and (y < ny):
			Fz0V.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])
	elif (z == nz):
		if x < nx: Ez1V.append([h,ind(x+1,y,z)])
		if y < ny: Ez1V.append([h,ind(x,y+1,z)])
		if (x < nx)  and (y < ny):
			Fz1V.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])

	if (y == 0):
		if x < nx: Ey0V.append([h,ind(x+1,y,z)])
		if z < nz: Ey0V.append([h,ind(x,y,z+1)])
		if (x < nx) and (z < nz):
			Fy0V.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])
	elif (y == ny):
		if x < nx: Ey1V.append([h,ind(x+1,y,z)])
		if z < nz: Ey1V.append([h,ind(x,y,z+1)])
		if (x < nx) and (z < nz):
			Fy1V.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])

	if (x == 0):
		if y < ny: Ex0V.append([h,ind(x,y+1,z)])
		if z < nz: Ex0V.append([h,ind(x,y,z+1)])
		if (y < ny) and (z < nz):
			Fx0V.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])
	elif (x == nx):
		if y < ny: Ex1V.append([h,ind(x,y+1,z)])
		if z < nz: Ex1V.append([h,ind(x,y,z+1)])
		if (y < ny) and (z < nz):
			Fx1V.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])

FbV = Fz0V+Fz1V+Fy0V+Fy1V+Fx0V+Fx1V
EbV = Ez0V+Ez1V+Ey0V+Ey1V+Ex0V+Ex1V
"""

"""
if __name__=="__main__" and DEBUG == True:
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FbV)))
	VIEW(hpc)
	hpc = EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EbV)))
	VIEW(hpc)
"""
	
# computation of the ∂2 operator on the boundary space
# ------------------------------------------------------------
print "start partial_2_b computation"
#partial_2_b = larBoundary(EbV,FbV)
print "end partial_2_b computation"



# computation of ∂3 operator on the image space
# ------------------------------------------------------------
print "start partial_3 computation"
partial_3 = larBoundary(FV,CV)
print "end partial_3 computation"



# ------------------------------------------------------------
# input from volume image (test: 250 x 250 x 250)
# ------------------------------------------------------------
out = []
Nx,Ny,Nz = imageHeight/imageDx, imageWidth/imageDx, imageDepth/imageDz
segFaces = set(["Fz0V","Fz1V","Fy0V","Fy1V","Fx0V","Fx1V"])

for inputIteration in range(imageDepth/imageDz):
	startImage = endImage
	endImage = startImage + imageDz
	xEnd, yEnd = 0,0
	theImage,colors,theColors = pngstack2array3d('SLICES2/', startImage, endImage, colors)
	print "\ntheColors =",theColors
	theColors = theColors.reshape(1,2)
	background = max(theColors[0])
	foreground = min(theColors[0])
	print "\n(background,foreground) =",(background,foreground)

	if __name__=="__main__" and DEBUG == True:
		print "\nstartImage, endImage =", (startImage, endImage)
	
	for i in range(imageHeight/imageDx):
		
		for j in range(imageWidth/imageDy):
			
			xStart, yStart = i * imageDx, j * imageDy
			xEnd, yEnd = xStart+imageDx, yStart+imageDy
			
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
			
			"""
			# Computation of (local) boundary to be removed by pieces
			# ------------------------------------------------------------
			pieceCoords = invertPiece(Nx,Ny,Nz)(count)
			if pieceCoords[2] == Nz-1: boundaryPlanes = ["Fz1V"]
			else:  boundaryPlanes = []

			if pieceCoords[0] == 0:  boundaryPlanes += ["Fx0V"]
			elif pieceCoords[0] == Nx-1:  boundaryPlanes += ["Fx1V"]
			if pieceCoords[1] == 0:  boundaryPlanes += ["Fy0V"]
			elif pieceCoords[1] == Ny-1:  boundaryPlanes += ["Fy1V"]
			if pieceCoords[2] == 0:  boundaryPlanes += ["Fz0V"]
			elif pieceCoords[2] == Nz-1:  boundaryPlanes += ["Fz1V"]
			
			"""
			#if __name__=="__main__" and DEBUG == True:
				#planesToRemove = list(segFaces.difference(boundaryPlanes))
				#FVtoRemove = CAT(map(eval,planesToRemove))


			count += 1
			

			# compute a quotient complex of chains with constant field
			# ------------------------------------------------------------
			chains3D = [[] for k in range(colors)]
			def addr(x,y,z): return x + (nx) * (y + (ny) * (z))
			for x in range(nx):
				for y in range(ny):
					for z in range(nz):
						if (image[x,y,z] == background):
							chains3D[1].append(addr(x,y,z))
						else:
							chains3D[0].append(addr(x,y,z))

            #if __name__=="__main__" and DEBUG == True:
                #print "\nchains3D =\n", chains3D


			# compute the boundary complex of the quotient cell
			# ------------------------------------------------------------
			objectBoundaryChain = larBoundaryChain(partial_3,chains3D[1])
			b2cells = csrChainToCellList(objectBoundaryChain)
			sup_cell_boundary = MKPOLS((V,[FV[f] for f in b2cells]))

			# remove the (local) boundary (shared with the piece boundary) from the quotient cell
			# ------------------------------------------------------------
			
			"""
			cellIntersection = matrixProduct(csrCreate([FV[f] for f in b2cells]),csrCreate(FVtoRemove).T)
			#print "\ncellIntersection =", cellIntersection
			cooCellInt = cellIntersection.tocoo()
			b2cells = [cooCellInt.row[k] for k,val in enumerate(cooCellInt.data) if val >= 4]
			"""
				

			# ------------------------------------------------------------
			# visualize the generated model  
			# ------------------------------------------------------------
			
			zStart = startImage - beginImageStack
			print "xStart, yStart, zStart =", xStart, yStart, zStart

			if __name__=="__main__":
				sup_cell_boundary = MKPOLS((V,[FV[f]  for f in b2cells]))
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
