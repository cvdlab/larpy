from lar import *
from scipy import *
import scipy
import numpy as np
from time import time
from pngstack2array3d import pngstack2array3d

# semplice esempio di array 3D
# ------------------------------------------------------------

image = ndarray(shape=(6,6,6),dtype=uint8)
image[0]=[
[0,0,0,3,2,2],
[0,3,0,3,3,2],
[0,3,3,3,3,2],
[1,1,3,3,2,2],
[1,1,1,3,3,3],
[1,1,1,3,3,3]]
image[1]=[
[0,0,0,2,2,2],
[0,0,0,3,3,2],
[3,3,3,3,3,2],
[3,3,3,3,2,2],
[1,1,1,1,3,3],
[1,1,1,1,3,3]]
image[2]=[
[0,0,0,0,2,2],
[0,3,3,3,3,2],
[3,3,3,2,2,2],
[3,3,3,2,2,1],
[1,1,1,4,4,1],
[1,1,1,4,4,1]]
image[3]=[
[2,2,0,0,0,0],
[2,2,3,3,0,0],
[3,3,3,3,1,1],
[1,4,4,2,2,1],
[1,1,1,2,2,1],
[1,1,1,2,2,1]]
image[4]=[
[2,2,2,2,0,0],
[2,2,3,3,0,0],
[0,0,3,3,1,1],
[0,0,4,4,2,2],
[0,0,0,2,2,2],
[0,0,0,2,2,2]]
image[5]=[
[2,2,2,1,1,1],
[3,3,3,4,1,1],
[3,3,3,4,1,1],
[0,0,0,4,2,2],
[0,0,0,2,2,2],
[0,0,0,2,2,2]]

nx,ny,nz = image.shape

# ------------------------------------------------------------
# input da immagine di volume
# ------------------------------------------------------------

colors = 2
image,colors = pngstack2array3d('SLICES2/', 430, 480, colors)
image = image[:,:50,:50]

nx,ny,nz = image.shape

# ------------------------------------------------------------
# complesso a celle cubiche
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

def theCell(coords):
    x,y,z = coords
    return [ind(x,y,z),ind(x+1,y,z),ind(x,y+1,z),ind(x,y,z+1),ind(x+1,y+1,z),
            ind(x+1,y,z+1),ind(x,y+1,z+1),ind(x+1,y+1,z+1)]

V = [[x,y,z] for z in range(nz+1) for y in range(ny+1) for x in range(nx+1) ]
CV = [theCell([x,y,z]) for z in range(nz) for y in range(ny) for x in range(nx)]
print "\nV =", V
print "\nCV =", CV
model = (V,CV)
#VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))

# ------------------------------------------------------------
# complesso di catene di campo costante
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

print "\nchains3D =\n", chains3D

"""
theColor = [RED,GREEN,BLUE,CYAN,MAGENTA,YELLOW]
VIEW(
    STRUCT(
[
    COLOR(RED)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[0]])),
    COLOR(GREEN)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[1]])),
    COLOR(BLUE)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[2]])),
    COLOR(CYAN)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[3]])),
    COLOR(MAGENTA)(STRUCT([ JOIN([MK(V[v]) for v in CV[cell] ]) for cell in chains3D[4]])),
]))
"""

# ------------------------------------------------------------
# complesso di catene di campo costante
# ------------------------------------------------------------
"""
csrSCC = csrCreate(chains3D)
#csrToMatrixRepresentation(csrSCC)[0]
csrCV = csrCreate(CV)

SCV = matrixProduct(csrSCC,csrCV)
for k in range(colors):
    print csrToMatrixRepresentation(SCV)[k]

vinfo = [(SCV.T[h].sum(), SCV.T[h].nnz) for h in range(len(V))]

for h in range(len(V)):
	print "v,sum,nnz =", (h, vinfo[h])
"""
# ------------------------------------------------------------
# costruzione della relazione VE (complesso originario)
# ------------------------------------------------------------
"""
EV = []
v2coords = invertIndex(nx,ny,nz)
for h in range(len(V)):
	x,y,z = v2coords(h)
	if x < nx: EV.append([h,ind(x+1,y,z)])
	if y < ny: EV.append([h,ind(x,y+1,z)])
	if z < nz: EV.append([h,ind(x,y,z+1)])

print "\nEV =",EV
#VIEW(STRUCT(MKPOLS((V,EV))))
"""
# ------------------------------------------------------------
# costruzione della relazione VF (complesso originario)
# ------------------------------------------------------------
FV = []
v2coords = invertIndex(nx,ny,nz)
for h in range(len(V)):
	x,y,z = v2coords(h)
	if (x < nx) and (y < ny): FV.append([h,ind(x+1,y,z),ind(x,y+1,z),ind(x+1,y+1,z)])
	if (x < nx) and (z < nz): FV.append([h,ind(x+1,y,z),ind(x,y,z+1),ind(x+1,y,z+1)])
	if (y < ny) and (z < nz): FV.append([h,ind(x,y+1,z),ind(x,y,z+1),ind(x,y+1,z+1)])

print "\nFV =",FV
#VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))
#VIEW(STRUCT(MKPOLS((V,FV))))

# ------------------------------------------------------------
# Drawing of supervertices
# ------------------------------------------------------------
'''
superverts = [k  for k, (valence,incidences) in enumerate(vinfo) 
	if (valence==1) or (valence==2 and incidences>1) or (valence==4 and incidences>2) or (valence==8 and incidences>1)]
sphere = SPHERE(0.1)([8,16])
supnodes = STRUCT(CONS(AA(T([1,2,3]))([V[v] for v in superverts]))(sphere))
VIEW(supnodes)
'''

# ------------------------------------------------------------
# Bordi delle supercelle
# ------------------------------------------------------------
print "\ncsrboundary_3 =","inizio"
csrboundary_3 = larBoundary(FV,CV)  # <<---------- BUG
print "\ncsrboundary_3 =",csrboundary_3


sup_cell_boundary = []
"""
for h in range(len(chains3D)):
	print "\nh =",h
	boundaryChain = larBoundaryChain(csrboundary_3,chains3D[h])
	b2cells = csrChainToCellList(boundaryChain)
	sup_cell_boundary += [MKPOLS((V,[FV[f] for f in b2cells]))]
"""
boundaryChain = larBoundaryChain(csrboundary_3,chains3D[1])
b2cells = csrChainToCellList(boundaryChain)
sup_cell_boundary += [MKPOLS((V,[FV[f] for f in b2cells]))]

print "visualizza"

mycolors = [RED,GREEN,BLUE,CYAN,MAGENTA,YELLOW]

VIEW(COLOR(BLUE)(S([1,2,3])([-1,-1,-1])(STRUCT(sup_cell_boundary[0]))))

"""
VIEW(COLOR(WHITE)(STRUCT(sup_cell_boundary[1])))

VIEW(STRUCT(
	[COLOR(mycolors[k])(supcell) for k in range(colors) for supcell in sup_cell_boundary[k]  if sup_cell_boundary[k] != []]
	# + [supnodes]
))
"""
