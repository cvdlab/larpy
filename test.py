from lar import *
from voxels2lar import *
from random import random
from blocks2solid import *

#-------------------------------------------------------------------
# data input of 2D image

solids = [[283,388,12,7],[288,395,10,16],[291,411,26,8],[309,401,12,10],[317,411,4,3],[298,405,3,6],[301,409,8,2],[296,419,10,3],[306,419,5,2],[286,395,2,6],[295,391,2,4],[298,399,2,6],[285,386,7,2],[284,395,2,2],[289,411,2,4],[317,414,3,2],[305,408,4,1],[321,402,1,8],[287,401,1,5],[298,396,1,3],[297,393,1,2],[282,391,1,2],[292,387,2,1],[294,419,2,2],[311,419,4,1],[317,416,2,1],[284,387,1,1],[288,385,1,1],[285,397,1,1],[308,407,1,1],[293,419,1,1],[317,417,1,1],[319,399,1,1],[315,397,1,1],[311,399,1,1],[310,400,2,1],[319,400,2,1],[312,398,7,3] ]
voids = [[299,381,27,16],[276,422,50,9],[276,398,10,24],[301,397,8,10],[322,397,4,25],[276,381,6,17],[282,381,17,4],[295,385,4,5],[286,415,4,7],[317,418,5,4],[311,420,6,2],[301,407,4,2],[286,411,3,4],[286,406,2,5],[290,419,3,3],[282,395,2,3],[282,385,2,3],[292,385,3,2],[297,390,2,3],[309,397,3,2],[319,397,3,2],[319,416,3,2],[320,414,2,2],[305,407,3,1],[299,397,2,2],[300,399,1,6],[298,393,1,3],[284,385,4,1],[306,421,5,1],[289,385,3,1],[312,397,3,1],[316,397,3,1],[282,388,1,3],[309,399,1,2],[321,399,1,3],[321,410,1,4],[293,421,3,1],[286,401,1,5],[282,393,1,2],[290,415,1,4],[284,386,1,1],[294,387,1,1],[295,390,2,1],[284,397,1,1],[293,420,1,1],[318,417,1,1],[315,419,2,1],[320,399,1,1],[310,399,1,1] ]
rects = solids + voids

#-------------------------------------------------------------------
# input testing

blocks = [[[rect[0],rect[1]],[rect[0]+rect[2], rect[1]+rect[3]]] for rect in rects]
boundary = [ [[276,381],[326,381]], [[276,381],[276,431]], [[276,431],[326,431]], [[326,326],[326,431]] ]

VIEW(STRUCT([COLOR(BLUE)(STRUCT([T([1,2])(rect[:2])(CUBOID(rect[2:])) for rect in solids])),
COLOR(RED)(STRUCT([T([1,2])(rect[:2])(CUBOID(rect[2:])) for rect in voids]))]))



#-------------------------------------------------------------------
# random generation of 3D image data

def block2DtoBlocks3D(block):
    p0,p1 = block
    z0 = -10
    z1 = int(18*random())-9
    z2 = int(18*random())-9
    z3 = 10
    if z1 == z2: z2 = z1+1
    elif z1 > z2: z1,z2 = z2,z1
    block0 = [p0+[z0],p1+[z1]]
    block1 = [p0+[z1],p1+[z2]]
    block2 = [p0+[z2],p1+[z3]]
    return [block0,block1,block2]

def block2DtoSingleBlock3D(block):
    p0,p1 = block
    z0 = -10
    z1 = 10
    block = [p0+[z0],p1+[z1]]
    return block

solid2D = [[[rect[0],rect[1]],[rect[0]+rect[2], rect[1]+rect[3]]] for rect in solids]
voids2D = [[[rect[0],rect[1]],[rect[0]+rect[2], rect[1]+rect[3]]] for rect in voids]


voids1,solids,voids2 = TRANS([block2DtoBlocks3D(block) for block in solid2D])
voids0 = [block2DtoSingleBlock3D(block) for block in voids2D]

voxels0 = AA(COLOR(RED))([T([1,2,3])(s[0])(CUBOID(VECTDIFF([s[1],s[0]])))  for s in voids1])
voxels1 = AA(COLOR(GREEN))([T([1,2,3])(s[0])(CUBOID(VECTDIFF([s[1],s[0]])))  for s in solids])
voxels2 = AA(COLOR(BLUE))([T([1,2,3])(s[0])(CUBOID(VECTDIFF([s[1],s[0]])))  for s in voids2])
voxels3 = AA(COLOR(WHITE))([T([1,2,3])(s[0])(CUBOID(VECTDIFF([s[1],s[0]])))  for s in voids0])

View(EXPLODE(1.2,1.2,1.2)(voxels0 + voxels1 + voxels2 + voxels3))
View(STRUCT(voxels0 + voxels1 + voxels2 + voxels3))


voxels0 = AA(COLOR(RED))([T([1,2,3])(s[0])(SKELETON(1)(CUBOID(VECTDIFF([s[1],s[0]]))))  for s in voids1])
voxels1 = AA(COLOR(GREEN))([T([1,2,3])(s[0])(SKELETON(1)(CUBOID(VECTDIFF([s[1],s[0]]))))  for s in solids])
voxels2 = AA(COLOR(BLUE))([T([1,2,3])(s[0])(SKELETON(1)(CUBOID(VECTDIFF([s[1],s[0]]))))  for s in voids2])
voxels3 = AA(COLOR(WHITE))([T([1,2,3])(s[0])(SKELETON(1)(CUBOID(VECTDIFF([s[1],s[0]]))))  for s in voids0])

View(EXPLODE(1.2,1.2,1.2)(voxels0 + voxels1 + voxels2 + voxels3))
View(STRUCT(voxels0 + voxels1 + voxels2 + voxels3))

boundary = [ [[276,381,-10],[326,381,10]], [[276,381,-10],[276,431,10]], [[276,431,-10],[326,431,10]], [[326,326,-10],[326,431,10]], [[276,381,-10],[326,431,-10]], [[276,381,10],[326,431,10]] ]

voxels = solids + voids0 + voids1 + voids2 + boundary
print "\nvoxels =\n", voxels

#-------------------------------------------------------------------
# cellular complex generation

def voxel2block(voxel):
    p0,p1 = voxel
    x,y,z = p0
    dx,dy,dz = AA(abs)(VECTDIFF([p1,p0]))
    return x,y,z,dx,dy,dz

def computeMins(blocks):
    return AA(min)(TRANS([[x+dx,y+dy,z+dz] for [x,y,z,dx,dy,dz] in blocks]))

def translateImageOnOrigin(voxels):
    firstPoints = [voxel[0] for voxel in voxels]
    tvect = AA(min)(TRANS(firstPoints))
    blocks = [ CAT([ VECTDIFF([voxel[0],tvect]), VECTDIFF(REVERSE(voxel)) ]) for voxel in voxels ]
    return blocks

blocks = translateImageOnOrigin(voxels)

#model = larFromImageBlocks(voxels)
model = lar3DFromImageBlocks(blocks)
#V,cells = model
V,faces = larSkeletons(model,dim=3)
F0V, F1V, F2V, F3V = faces
#V = model[0]
VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V[:-6])) ))
#VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V)) ))
VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V)) ))
VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))


#-------------------------------------------------------------------
# boundaries extraction
    
csrBoundary_3 = larBoundary(F2V,F3V)
print "\ncsrBoundary_3.shape =", csrBoundary_3.shape
chain_2 = larBoundaryChain(csrBoundary_3,range(len(solids)))
print "\nchain_2.T =", csrToMatrixRepresentation(chain_2.T)
_2cells = csrExtractAllGenerators(csrBoundary_3)
print "\n_2cells =", _2cells        # list of 3-cells, given as lists of (boundary) 2-cells

# boundary 2-chain computation
boundary_2_cells = csrChainToCellList( chain_2 )
print "\nboundary_2_cells =\n",boundary_2_cells
# boundary 2-chain visualization
boundary2D = AA(JOIN)([ AA(MK)([V[v] for v in F2V[f]]) for f in boundary_2_cells ])
VIEW(EXPLODE(1.2,1.2,1.2)(boundary2D))
VIEW(STRUCT(boundary2D))

# boundary 2-chain of the empty space
emptyChain = len(voxels[:-6])
chain_2_emptyCells = larBoundaryChain(csrBoundary_3,range(len(solids),len(voxels[:-6])))
boundary_2_emptyCells = csrChainToCellList( chain_2_emptyCells )
boundary2D_emptyCells = AA(JOIN)([ AA(MK)([V[v] for v in F2V[f]]) for f in boundary_2_emptyCells ])
VIEW(EXPLODE(1.2,1.2,1.2)(boundary2D_emptyCells))
VIEW(STRUCT(boundary2D_emptyCells))

# boundary 2-chain solidification
VIEW(SOLIDIFY(T([1,2])([-301.0, -406.0])(STRUCT(boundary2D))))

# 1-chain of 2-boundary
csrBoundary_2 = larBoundary(F1V,F2V)
n = csrGetNumberOfColumns(csrBoundary_2)
csrChain_1 = sum([csrProduct(csrBoundary_2, csrCreateUnitChain(n,f)) for f in boundary_2_cells ])
chain_1 = csrChainToCellList(csrPredFilter(csrChain_1,C(EQ)(2)))
model = (V,[F1V[e] for e in chain_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))


