from lar import *
import collections
from time import time
from image2blocks import *
"""
Fast conversion from 2D blocks (rectangles: blocks = [[x,y,dx,dy]]) to LAR 2D input (model: V,FV)
"""
#-- API ---------------------------------------------------------------------------------------------

def write2DBlock(cooStore):
    def write2DBlock0(block):
        cooStore.append([0,0,2])
        x,y,dx,dy = block
        i = x
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        i = x+dx
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        j = y
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        j = y+dy
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        return cooStore
    return write2DBlock0

def write3DBlock(store):
    def write3DBlock0(block):
        x,y,z,dx,dy,dz = block
        for j in range(y,y+dy+1):
            for k in range(z,z+dz+1):
                store[x,j,k] += 1
                store[x+dx,j,k] += 1
        for k in range(z,z+dz+1):
            for i in range(x,x+dx+1):
                store[i,y,k] += 1
                store[i,y+dy,k] += 1
        for i in range(x,x+dx+1):
            for j in range(y,y+dy+1):
                store[i,j,z] += 1
                store[i,j,z+dz] += 1
        return store
    return write3DBlock0


def read2DBlock(lilStore):
    def read2DBlock0(block):
        x,y,dx,dy = block
        outBlock = [[(i,y) for i in range(x,x+dx) if lilStore[i,y] > 2]]
        outBlock.append([(x+dx,j) for j in range(y,y+dy) if lilStore[x+dx,j] > 2])
        outBlock.append([(i,y+dy) for i in range(x+dx,x,-1) if lilStore[i,y+dy] > 2])
        outBlock.append([(x,j) for j in range(y+dy,y,-1) if lilStore[x,j] > 2])
        return CAT(outBlock)
    return read2DBlock0

def read3DBlock(store):
    def read3DBlock0(block):
        x,y,z,dx,dy,dz = block
        outBlock = [[(x,j,k) for j in range(y,y+dy+1) for k in range(z,z+dz+1) if store[x,j,k] > 3]]
        outBlock.append([(x+dx,j,k) for j in range(y,y+dy+1) for k in range(z,z+dz+1) if store[x+dx,j,k] > 3])
        outBlock.append([(i,y,k) for k in range(z,z+dz+1) for i in range(x,x+dx+1) if store[i,y,k] > 3])
        outBlock.append([(i,y+dy,k) for k in range(z,z+dz+1) for i in range(x,x+dx+1) if store[i,y+dy,k] > 3])
        outBlock.append([(i,j,z) for i in range(x,x+dx+1) for j in range(y,y+dy+1) if store[i,j,z] > 3])
        outBlock.append([(i,j,z+dz) for i in range(x,x+dx+1) for j in range(y,y+dy+1) if store[i,j,z+dz] > 3])
        return CAT(outBlock)
    return read3DBlock0


def lar2DFromImageBlocks(blocks):
    cooStore = []
    for block in blocks:
        write2DBlock(cooStore)(block)
    lilStore = csrCreateFromCoo(cooStore).tolil()
    updatedBlock = [read2DBlock(lilStore)(block) for block in blocks]
    verts = collections.OrderedDict(); k = 0
    for block in updatedBlock:
        for vert in block:
            if vert not in verts:
                verts[vert] = k
                k += 1
    V = AA(list)(AA(AA(float))(verts.keys()))
    F2V = [[verts[vert] for vert in block] for block in updatedBlock]
    model = V,F2V
    return model

def lar3DFromImageBlocks(blocks):
    def computeShape(blocks):
        return AA(max)(TRANS([[x+dx,y+dy,z+dz] for [x,y,z,dx,dy,dz] in blocks]))
    ax,ay,az = computeShape(blocks)
    store = zeros(shape=(ax+1,ay+1,az+1),dtype=int)
    [ write3DBlock(store)(block) for block in blocks ]
    updatedBlocks = [ read3DBlock(store)(block) for block in blocks ]
    verts = collections.OrderedDict(); k = 0
    for block in updatedBlocks:
        for vert in block:
            if vert not in verts:
                verts[vert] = k
                k += 1
    V = verts.keys()
    F3V = [[verts[vert] for vert in block] for block in updatedBlocks]
    model = V,F3V
    return model

"""
def lar3DFromImageBlocks(blocks):
    outBlocks = []
    verts = collections.OrderedDict()
    k = 0
    for block in blocks:
        x0,y0,z0,dx,dy,dz = block
        x1,y1,z1 = x0+dx, y0+dy, z0+dz
        for vert in [(x0,y0,z0), (x1,y0,z0), (x0,y1,z0), (x0,y0,z1), (x0,y1,z1), (x1,y0,z1), (x1,y1,z0), (x1,y1,z1)]:
            if vert not in verts:
                verts[vert] = k
                k += 1
    for block in blocks:
        x0,y0,z0,dx,dy,dz = block
        x1,y1,z1 = x0+dx, y0+dy, z0+dz
        outBlock = []
        outBlock.append([verts[x0,y,z] for y in range(y0,y1+1) for z in range(z0,z1+1) if (x0,y,z) in verts])
        outBlock.append([verts[x,y0,z] for z in range(z0,z1+1) for x in range(x0,x1+1) if (x,y0,z) in verts])
        outBlock.append([verts[x,y,z0] for x in range(x0,x1+1) for y in range(y0,y1+1) if (x,y,z0) in verts])
        outBlock.append([verts[x1,y,z] for y in range(y0,y1+1) for z in range(z0,z1+1) if (x1,y,z) in verts])
        outBlock.append([verts[x,y1,z] for z in range(z0,z1+1) for x in range(x0,x1+1) if (x,y1,z) in verts])
        outBlock.append([verts[x,y,z1] for x in range(x0,x1+1) for y in range(y0,y1+1) if (x,y,z1) in verts])
        outBlock = list(set(CAT(outBlock)))
        outBlocks.append(outBlock)

    model = verts.keys(), outBlocks
    return model
"""

if __name__ == "__main__":

    #-- Test data ---------------------------------------------------------------------------------------------

    image = mahotas.imread('test1.bmp')
    pylab.imshow(image)
    pylab.show()

    solid,empty = BMPimage2blocks(image)
    
    nsolid = len(solid)
    nempty = len(empty)
    blocks = solid + empty

    start = time()
    V,F2V = lar2DFromImageBlocks(blocks)
    end = time()

    print "\ntotal time (sec) =", end-start
    print "\nV = ", V
    print "\nF2V = ", F2V, "\n"

    solids = AA(COLOR(BLUE))(MKPOLS((V,F2V[:nsolid])))
    voids = AA(COLOR(RED))(MKPOLS((V,F2V[nsolid:])))

    VIEW(EXPLODE(1.2,1.2,1.2)(solids + voids))


    Vmy = [k for k,v in enumerate(V) if v[0]==0.0]
    VMy = [k for k,v in enumerate(V) if v[0]==50.0]
    Vxm = [k for k,v in enumerate(V) if v[1]==0.0]
    VxM = [k for k,v in enumerate(V) if v[1]==50.0]

    F2V = F2V + [Vmy,VMy,Vxm,VxM]

    model = (V,F2V)
    V,faces = larSkeletons(model,dim=2,grid=False)
    F0V, F1V, F2V = faces
    print "AA(LEN)([F0V, F1V, F2V]) =", AA(LEN)([F0V, F1V, F2V])
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))





    csrBoundary_2 = larBoundary(F1V,F2V)
    print "\ncsrBoundary_2.shape =", csrBoundary_2.shape
    chain_1 = larBoundaryChain(csrBoundary_2,range(len(solids)))
    print "\nchain_1.T =", csrToMatrixRepresentation(chain_1.T)
    _1cells = csrExtractAllGenerators(csrBoundary_2)
    print "\n_1cells =", _1cells        # list of 2-cells, given as lists of (boundary) 1-cells

    # boundary 1-chain computation
    boundary_1_cells = csrChainToCellList( chain_1 )
    print "\nboundary_1_cells =\n",boundary_1_cells

    # boundary 1-chain visualization
    boundary1D = AA(POLYLINE)([[V[v] for v in F1V[e]] for e in boundary_1_cells ])
    VIEW(EXPLODE(1.2,1.2,1.2)(boundary1D))
    VIEW(STRUCT(boundary1D))
