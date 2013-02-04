from lar import *



#-------------------------------------------------------------------
# dimension-independent conversion from discrete image to LAR

def larFromImageBlocks(blocks):
    dim = len(blocks[0][0])
    blocks = AA(CAT)(blocks)
    bounds = TRANS(blocks)
    cells = [[] for k in range(len(blocks))]
    minMaxCoords = zip(AA(min)(bounds[:dim]), AA(max)(bounds[dim:]))
    gridPoints = CART(AA(FROMTO)(AA(list)(minMaxCoords)))
    counter,V,i = [],[],0
    for point in gridPoints:
        for k,cell in enumerate(blocks):
            pmin,pmax = cell[:dim],cell[dim:]
            classify = AND(CAT([AA(ISGE)(TRANS([pmin,point])),AA(ISLE)(TRANS([pmax,point]))]))
            if classify: counter.append(k)
        if len(counter) >= dim+1:
            [cells[k].append(i) for k in counter]
            V += [point]
            i += 1
        counter = []
    return V,cells

#-------------------------------------------------------------------
# 2D image example

if __name__=="__main__":
    blocks = [ [[0,0],[5,10]], [[5,0],[9,3]], [[9,0],[13,3]], [[5,3],[8,10]],  [[8,3],[13,10]], [[0,10],[9,12]], [[9,10],[13,12]], [[0,0],[0,12]], [[0,0],[13,0]], [[13,0],[13,12]], [[0,12],[13,12]] ]
    
    model = larFromImageBlocks(blocks)
    V,cells = model
    V,faces = larSkeletons(model,dim=2)
    F0V, F1V, F2V = faces
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))

#-------------------------------------------------------------------
# 3D image example

if __name__=="__main__":
    blocks = [ [[0,0,0],[5,10,3]], [[5,0,0],[9,3,3]], [[9,0,0],[13,3,3]], [[5,3,0],[8,10,3]],  [[8,3,0],[13,10,3]], [[0,10,0],[9,12,3]], [[9,10,0],[13,12,3]], [[0,0,0],[0,12,3]], [[0,0,0],[13,0,3]], [[13,0,0],[13,12,3]], [[0,12,0],[13,12,3]], [[0,0,0],[13,12,0]], [[0,0,3],[13,12,3]] ]
    
    model = larFromImageBlocks(blocks)
    V,cells = model
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V[:-6])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))


#-------------------------------------------------------------------
# bigger 2D image example

if __name__=="__main__":

    rects = [
    [283,388,12,7],[288,395,10,16],[291,411,26,8],[309,401,12,10],[317,411,4,3],[298,405,3,6],[301,409,8,2],[296,419,10,3],[306,419,5,2],[286,395,2,6],[295,391,2,4],[298,399,2,6],[285,386,7,2],[284,395,2,2],[289,411,2,4],[317,414,3,2],[305,408,4,1],[321,402,1,8],[287,401,1,5],[298,396,1,3],[297,393,1,2],[282,391,1,2],[292,387,2,1],[294,419,2,2],[311,419,4,1],[317,416,2,1],[284,387,1,1],[288,385,1,1],[285,397,1,1],[308,407,1,1],[293,419,1,1],[317,417,1,1],[319,399,1,1],[315,397,1,1],[311,399,1,1],[310,400,2,1],[319,400,2,1],[312,398,7,3],
    [299,381,27,16],[276,422,50,9],[276,398,10,24],[301,397,8,10],[322,397,4,25],[276,381,6,17],[282,381,17,4],[295,385,4,5],[286,415,4,7],[317,418,5,4],[311,420,6,2],[301,407,4,2],[286,411,3,4],[286,406,2,5],[290,419,3,3],[282,395,2,3],[282,385,2,3],[292,385,3,2],[297,390,2,3],[309,397,3,2],[319,397,3,2],[319,416,3,2],[320,414,2,2],[305,407,3,1],[299,397,2,2],[300,399,1,6],[298,393,1,3],[284,385,4,1],[306,421,5,1],[289,385,3,1],[312,397,3,1],[316,397,3,1],[282,388,1,3],[309,399,1,2],[321,399,1,3],[321,410,1,4],[293,421,3,1],[286,401,1,5],[282,393,1,2],[290,415,1,4],[284,386,1,1],[294,387,1,1],[295,390,2,1],[284,397,1,1],[293,420,1,1],[318,417,1,1],[315,419,2,1],[320,399,1,1],[310,399,1,1],
     ]

    blocks = [[[rect[0],rect[1]],[rect[0]+rect[2], rect[1]+rect[3]]] for rect in rects] + [ [[276,381],[326,381]], [[276,381],[276,431]], [[276,431],[326,431]], [[326,326],[326,431]] ]

    model = larFromImageBlocks(blocks)
    V,cells = model
    V,faces = larSkeletons(model,dim=2)
    F0V, F1V, F2V = faces
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))


#-------------------------------------------------------------------
# boundaries extraction

if __name__=="__main__":

    csrBoundary_2 = larBoundary(F1V,F2V)
    print "\ncsrBoundary_2.shape =", csrBoundary_2.shape
    chain_1 = larBoundaryChain(csrBoundary_2,range(38))
    print "\nchain_1.T =", csrToMatrixRepresentation(chain_1.T)
    _1cells = csrExtractAllGenerators(csrBoundary_2)
    print "\n_1cells =", _1cells        # list of 2-cells, given as lists of (boundary) 1-cells

    # boundary 1-chain computation
    boundary_1_cells = csrChainToCellList( chain_1 )
    print "\nboundary_1_cells =\n",boundary_1_cells
    # boundary 1-chain visualization
    boundary1D = AA(POLYLINE)([[V[v] for v in F1V[e]] for e in boundary_1_cells ])
    VIEW(EXPLODE(1.2,1.2,1.2)(boundary1D))
    View(SOLIDIFY(T([1,2])([-301.0, -406.0])(STRUCT(boundary1D))))

    chain_2 = larBoundaryChain(csrBoundary_2,range(38,len(F2V[:-4])))
    # boundary 1-chain computation
    boundary_1_void = csrChainToCellList( chain_2 )
    print "\nboundary_1_cells =\n",boundary_1_void
    # boundary 1-chain visualization
    boundary1D = AA(POLYLINE)([[V[v] for v in F1V[e]] for e in boundary_1_void ])
    VIEW(EXPLODE(1.2,1.2,1.2)(boundary1D))
    MED([1,2])(STRUCT(boundary1D))
    VIEW(SOLIDIFY(T([1,2])([-301.0, -406.0])(STRUCT(boundary1D))))

