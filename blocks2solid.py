from lar import *
import collections
from time import time
from image2blocks import *
"""
Fast conversion from 2D blocks (rectangles: blocks = [[x,y,dx,dy]]) to LAR 2D input (model: V,FV)
"""


if __name__ == "__main__":

    #-- Test data ---------------------------------------------------------------------------------------------

    image = mahotas.imread('test1.bmp')
    pylab.imshow(image)
    pylab.show()

    start = time()
    solid,empty = BMPimage2blocks(image)
    end = time()
    print "\ntotal time (sec) BMPimage2blocks:", end-start
    
    nsolid = len(solid)
    nempty = len(empty)
    blocks = solid + empty

    start = time()
    V,F2V = lar2DFromImageBlocks(blocks)
    end = time()
    print "\ntotal time (sec) lar2DFromImageBlocks:", end-start
    print "\nV = ", V
    print "\nF2V = ", F2V, "\n"

    solids = AA(COLOR(BLUE))(MKPOLS((V,F2V[:nsolid])))
    voids = AA(COLOR(RED))(MKPOLS((V,F2V[nsolid:-4])))

    VIEW(EXPLODE(1.2,1.2,1.2)(solids + voids))

    model = (V,F2V)
    V,faces = larSkeletons(model,dim=2)
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
