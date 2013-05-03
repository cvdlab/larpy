from lar import *
from scipy import *
import scipy
import numpy as np

#-------------------------------------------------------------
# data generation
#-------------------------------------------------------------
nx,ny,nz = 6,6,6
colors = 6

image = ndarray(shape=(nx,ny,nz),dtype=np.dtype(np.int8))
for x in range(nx):
	for y in range(ny):
		for z in range(nz):
			# random generation (no filtering)
			#-------------------------------------------------------------
			image[x,y,z] = trunc(random.random() * colors)
			# random generation (filtering)
			#-------------------------------------------------------------
			done = False
			while not done:
				value = trunc(random.random() * colors)
				if image[x-1,y,z]==value or image[x,y-1,z]==value or image[x,y,z-1]==value or \
				x==0 or y==0 or z==0:
					image[x,y,z] = value
					done = True
print "\nimage =\n", image

#-------------------------------------------------------------
# LAR building from data
#-------------------------------------------------------------
def ind(x,y,z): return x + (nx+1) * (y + (ny+1) * (z))
def invertIndex(nx,ny,nz):
	def invertIndex0(offset):
		a0, b0 = offset / nx, offset % nx
		a1, b1 = a0 / ny, a0 % ny
		a2, b2 = a1 / nz, a1 % nz
		return b0,b1,b2
	return invertIndex0
	

def cell(coords):
	x,y,z = coords
	return [ind(x,y,z),ind(x+1,y,z),ind(x,y+1,z),ind(x,y,z+1),ind(x+1,y+1,z),
			ind(x+1,y,z+1),ind(x,y+1,z+1),ind(x+1,y+1,z+1)]
	
V = [[x,y,z] for z in range(nz+1) for y in range(ny+1) for x in range(nx+1) ]
CV = [cell([x,y,z]) for z in range(nz) for y in range(ny) for x in range(nx)]

print "\nV =", V
print "\nCV =", CV

model = V,CV
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))

xmin = [ind(0,y,z)  for z in range(nz+1) for y in range(ny+1)]
xmax = [ind(nx,y,z)  for z in range(nz+1) for y in range(ny+1)]
ymin = [ind(x,0,z)  for x in range(nx+1) for z in range(nz+1)]
ymax = [ind(x,ny,z)  for x in range(nx+1) for z in range(nz+1)]
zmin = [ind(x,y,0)  for x in range(nx+1) for y in range(ny+1)]
zmax = [ind(x,y,nz)  for x in range(nx+1) for y in range(ny+1)]

CV = CV+[xmin,xmax,ymin,ymax,zmin,zmax]
model = V,CV
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))

V,FV = larFacets(model,dim=3)
FV = [facet for facet in FV if len(facet)==4]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

CV = CV[:-6]
model = V,CV
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))

#-------------------------------------------------------------
# 3D chains extraction (C3 matrix, by rows)
#-------------------------------------------------------------
chains3D = [[] for k in range(colors)]
def k(x,y,z): return x + (nx) * (y + (ny) * (z))
for x in range(nx):
	for y in range(ny):
		for z in range(nz):
			chain = image[x,y,z]
			chains3D[chain].append(k(x,y,z))

print "\nchains3D =\n", chains3D

# 3D chains visualization
#-------------------------------------------------------------
theColor = [RED,GREEN,BLUE,CYAN,MAGENTA,YELLOW]
VIEW(
	#STRUCT(
	EXPLODE(1.4,1.4,1.4)(
[
	COLOR(RED)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[0]])),
	COLOR(GREEN)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[1]])),
	COLOR(BLUE)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[2]])),
	COLOR(CYAN)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[3]])),
	COLOR(MAGENTA)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[4]])),
	COLOR(YELLOW)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[5]]))
]))


#-------------------------------------------------------------
# Boundary operator
#-------------------------------------------------------------
csrBoundary_3 = larBoundary(FV,CV)
chain_2 = larBoundaryChain(csrBoundary_3, chains3D[0])
boundary_2_cells = csrChainToCellList( chain_2 )
print "\nboundary_2_cells =",boundary_2_cells
boundary2D = AA(JOIN)([ AA(MK)([V[v] for v in FV[f]]) for f in boundary_2_cells ])
VIEW(EXPLODE(1.2,1.2,1.2)(boundary2D))

# -------------------------------------------------------------
# Extraction of vertices of "supercells" (vertex method)
# -------------------------------------------------------------


F3V = [[7,22,23,24,13,14,16,17],
[22,5,24,25,14,15,17,18],
[23,24,6,26,16,17,19,20],
[24,25,26,3,17,18,20,21],
[13,14,16,17,4,8,9,10],
[14,15,17,18,8,1,10,11],
[16,17,19,20,9,10,2,12],
[17,18,20,21,10,11,12,0]#,
# # exterior 3-cells
# [7,22,5,13,14,15,4,8,1],
# [7,23,6,13,16,19,4,9,2],
# [7,22,5,23,24,25,6,26,3],
# [6,26,3,19,20,21,2,12,0],
# [5,25,3,15,18,21,1,11,0],
# [4,8,1,9,10,11,2,12,0]
]

F2V = [[0, 10, 11, 12],
[0, 11, 18, 21],
[0, 12, 20, 21],
[1, 8, 10, 11],
[1, 8, 14, 15],
[1, 11, 15, 18],
[2, 9, 10, 12],
[2, 9, 16, 19],
[2, 12, 19, 20],
[3, 18, 21, 25],
[3, 20, 21, 26],
[3, 24, 25, 26],
[4, 8, 9, 10],
[4, 8, 13, 14],
[4, 9, 13, 16],
[5, 14, 15, 22],
[5, 15, 18, 25],
[5, 22, 24, 25],
[6, 16, 19, 23],
[6, 19, 20, 26],
[6, 23, 24, 26],
[7, 13, 14, 22],
[7, 13, 16, 23],
[7, 22, 23, 24],
# facets ortho to x axis (24:28)
[9, 10, 16, 17],
[10, 11, 17, 18],
[16, 17, 23, 24],
[17, 18, 24, 25],
# facets ortho to y axis (28:32)
[8, 10, 14, 17],
[10, 12, 17, 20],
[14, 17, 22, 24],
[17, 20, 24, 26],
# facets ortho to z axis (32:36)
[13, 14, 16, 17],
[14, 15, 17, 18],
[16, 17, 19, 20],
[17, 18, 20, 21]]

def boundaryVertex(bx,by,bz):
	vcoords = invertIndex(bx+1,by+1,bz+1)
	def boundaryVertex0(v):
		x,y,z = vcoords(v)
		constraints = AA(lambda test: (0,1)[test])([x==0, x==bx, y==0, y==by, z==0, z==bz])
		return constraints,x,y,z
	return boundaryVertex0

def count(theSet):
	def count0(subset):
		return [sum([1 for k in theSet if k==elem]) for elem in subset ]
	return count0


def single(v,cells,supcells):
	nodes = {
		0:[1,2,3], 1:[0,4,5], 2:[0,4,6], 3:[0,5,6], 
		4:[1,2,7], 5:[1,3,7], 6:[2,3,7], 7:[4,5,6],
	}
	supsets = [[] for color in range(colors)]
	for cell,col in enumerate(supcells):
		supsets[col] = supsets[col] + [cell] 

	supnodes = dict([ ( tuple(key), list(set(CAT([nodes[k] for k in key])).difference(key) ) ) 
					for key in supsets if [nodes[k] for k in key] != []])
	m = len(supnodes)
	insert = OR([True for key in supnodes if len(key)==1])
	return insert


def vertexTest(supcells,vsc):
	for color in vsc:
		chain_on_v = [cell for cell,col in enumerate(supcells) if col==color]
		chain_2 = larBoundaryChain(csrBoundary_3,[0,1,2,3])
		chain_on_v_boundary = [k for k,cell in enumerate(chain_2) if cell!=[]]
		normal2x = sum([1 for facet in chain_on_v_boundary if facet in range(24,28)])
		normal2y = sum([1 for facet in chain_on_v_boundary if facet in range(28,32)])
		normal2z = sum([1 for facet in chain_on_v_boundary if facet in range(32,36)])
		if normal2x*normal2y*normal2z != 0:  return True


verts = len(V)
csrVF = csrCreate([],shape=(verts,colors))
vcoords = boundaryVertex(nx,ny,nz)
VSC = [] 
for v in range(verts):
	constraints,x,y,z = vcoords(v)
	vtype = sum(constraints)
	label = "".join(AA(str)(constraints))
	if   vtype == 3:
		# corner vertex ----------------------------------
		cells = {
		"010101": [ind(nx,ny,nz)],
		"010110": [ind(nx,ny,0)],
		"011001": [ind(nx,0,nz)],
		"011010": [ind(nx,0,0)],
		"100101": [ind(0,ny,nz)],
		"100110": [ind(0,ny,0)],
		"101001": [ind(0,0,nz)],
		"101010": [ind(0,0,0)]
		}[label]
	elif vtype == 2:
		# edge vertex ----------------------------------
		cells = {
		"000101": [ind(x-1,ny,nz),ind(x,ny,nz)],
		"000110": [ind(x-1,ny,0),ind(x,ny,0)],
		"001001": [ind(x-1,0,nz),ind(x,0,nz)],
		"001010": [ind(x-1,0,0),ind(x,0,0)],
		"010001": [ind(nx,y-1,nz),ind(nx,y,nz)],
		"010010": [ind(nx,y-1,0),ind(nx,y,0)],
		"010100": [ind(nx,ny,z-1),ind(nx,ny,z)],
		"011000": [ind(nx,0,z-1),ind(nx,0,z)],
		"100001": [ind(0,y-1,nz),ind(0,y,nz)],
		"100010": [ind(0,y-1,0),ind(0,y,0)],
		"100100": [ind(0,ny,z-1),ind(0,ny,z)],
		"101000": [ind(0,0,z-1),ind(0,0,z)]
		}[label]
	elif vtype == 1:
		# face vertex ----------------------------------
		cells = {
		"000001": [ind(x,y,nz),ind(x-1,y,nz),ind(x,y-1,nz),ind(x-1,y-1,nz)],
		"000010": [ind(x,y,0),ind(x-1,y,0),ind(x,y-1,0),ind(x-1,y-1,0)],
		"000100": [ind(x,ny,z),ind(x-1,ny,z),ind(x,ny,z-1),ind(x-1,ny,z-1)],
		"001000": [ind(x,0,z),ind(x-1,0,z),ind(x,0,z-1),ind(x-1,0,z-1)],
		"010000": [ind(nx,y,z),ind(nx,y-1,z),ind(nx,y,z-1),ind(nx,y-1,z-1)],
		"100000": [ind(0,y,z),ind(0,y-1,z),ind(0,y,z-1),ind(0,y-1,z-1)]
		}[label]
	elif vtype == 0:
		# interior vertex ----------------------------------
		cells = {
		"000000": [ind(x-1,y-1,z-1),ind(x-1,y,z-1),ind(x,y-1,z-1),ind(x,y,z-1),
					ind(x-1,y-1,z),ind(x-1,y,z),ind(x,y-1,z),ind(x,y,z)]
		}[label]
	supcells = [image[invertIndex(nx,ny,nz)(cell)].tolist() for cell in cells]
	vsc = [supcell for supcell in range(colors) if supcell in supcells]
	if vtype == 3: theVSC = vsc
	elif vtype==2 and len(vsc)==1: theVSC = []
	elif vtype==2 and len(vsc)==2: theVSC = vsc
	elif vtype==1 and len(vsc)==1: theVSC = []
	elif vtype==1 and len(vsc)==2 and PROD(count(supcells)(vsc))%2==1: theVSC = vsc
	elif vtype==1 and len(vsc)==2 and cells[1]==cells[2]: theVSC = vsc
	elif vtype==1 and len(vsc)==2 and cells[1]!=cells[2]: theVSC = []
	elif vtype==1 and len(vsc)==3: theVSC = vsc
	elif vtype==1 and len(vsc)==4: theVSC = vsc
	elif vtype==0 and len(vsc)==1: theVSC = []
	elif vtype==0 and vertexTest(supcells,vsc): theVSC = vsc
	else: theVSC = []
	
	VSC += [theVSC]
	print v,vtype,label,cells,supcells,vsc,theVSC

# 
# # -------------------------------------------------------------
# # Extraction of vertices of "supercells" (vector method)
# # -------------------------------------------------------------
# def csrVertFilter(CSRm):
# 	for k,row in enumerate(CSRm):
# 		theRow = row[0].tocoo().data
# 		theCols = row[0].tocoo().col
# 		odds = sum([val%2 for val in theRow])
# 		if odds > 0:
# 			for j in theCols: CSRm[k,j] = 1 
# 		else:
# 			for j in theCols: CSRm[k,j] = 0
# 	return CSRm
# 
# faces,verts = len(FV),len(V)
# csrFV = csrCreate(FV,shape=(faces,verts))
# csrSCV = csrCreate([],shape=(verts, 1)).tocsc()
# for supercell in range(colors):
# 	chain_2 = larBoundaryChain(csrBoundary_3, chains3D[supercell])	# 2-cells in each supercell
# 	csrSFV = matrixProduct(csrFV.T,chain_2) 						# vertices in each supercell 2-boundary
# 	csrSCV = csrAppendByColumn(csrSCV,csrSFV)						# vertices in all supercell 2-boundaryies
# 
# csrSCV = csrSCV.tocsr()[:,1:]
# csrSCV = csrVertFilter(csrSCV)
# SCV = [supercell.tocoo().col.tolist() for supercell in csrSCV.T]
# print "\nSCV =",SCV
# 
# for supercell in range(colors):
# 	VIEW(STRUCT(
# 		[COLOR(theColor[supercell])(STRUCT([ JOIN([MK(V[vert]) for vert in CV[cell] ]) for cell in chains3D[supercell]]))
# 			] + MKPOLS((V,AA(LIST)(csrSCV.T[supercell].tocoo().col)))
# 	))
# 
# print len([sum(row.tocoo().data) for row in csrSCV if sum(row.tocoo().data)==0 ]) # number of filtered vertices
# 
# #-------------------------------------------------------------
# # Building of LAR "supercomplex"
# #-------------------------------------------------------------
# model = V,SCV
# VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
# 
# #csrSCV = csrCreate(SCV,shape=(colors,verts))
# 
# SCV = SCV+[xmin,xmax,ymin,ymax,zmin,zmax]
# model = V,SCV
# V,faces = larSkeletons(model,dim=3)
# F0V, F1V, F2V, F3V = faces
# VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V[:-6])) ))
# VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V)) ))
# VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
# #VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))
# 
