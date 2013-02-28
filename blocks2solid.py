from lar import *
import collections
from time import time
from image2blocks import *
"""
Fast conversion from 2D blocks (rectangles: blocks = [[x,y,dx,dy]]) to LAR 2D input (model: V,FV)
"""
#-- API ---------------------------------------------------------------------------------------------

def writeBlock(cooStore):
    def writeBlock0(block):
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
    return writeBlock0


def readBlock(lilStore):
    def readBlock0(block):
        x,y,dx,dy = block
        outBlock = [[(i,y) for i in range(x,x+dx) if lilStore[i,y] > 2]]
        outBlock.append([(x+dx,j) for j in range(y,y+dy) if lilStore[x+dx,j] > 2])
        outBlock.append([(i,y+dy) for i in range(x+dx,x,-1) if lilStore[i,y+dy] > 2])
        outBlock.append([(x,j) for j in range(y+dy,y,-1) if lilStore[x,j] > 2])
        return CAT(outBlock)
    return readBlock0


def lar2DFromImageBlocks(blocks):
    cooStore = []
    for block in blocks:
        writeBlock(cooStore)(block)
    lilStore = csrCreateFromCoo(cooStore).tolil()
    updatedBlock = [readBlock(lilStore)(block) for block in blocks]
    verts = collections.OrderedDict(); k = 0
    index = 0
    for block in updatedBlock:
        for vert in block:
            if vert not in verts:
                verts[vert] = k
                k += 1
    V = AA(list)(AA(AA(float))(verts.keys()))
    FV = [[verts[vert] for vert in block] for block in updatedBlock]
    model = V,FV
    return model

if __name__ == "__main__":

    #-- Test data ---------------------------------------------------------------------------------------------
    
    solid,empty = BMPimage2blocks('test1.bmp')
    
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
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))





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
