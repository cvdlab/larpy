# -*- coding: utf-8 -*-
import collections
import scipy.sparse
from scipy import zeros,arange,mat
from scipy.sparse import vstack,hstack,csr_matrix,lil_matrix,triu
from scipy.spatial import Delaunay
from scipy.linalg import *

from pyplasm import *
#------------------------------------------------------------------
#--geometry layer (using PyPlasm)----------------------------------
#------------------------------------------------------------------

def bezier(points):
    return MAP(BEZIERCURVE(points))(INTERVALS(1)(20))

# convex combination operator
CCOMB = COMP([ SCALARVECTPROD,
              CONS([ COMP([ DIV, CONS([K(1),LEN]) ]), VECTSUM ]) ])

def EXPLODE (sx,sy,sz):
    def explode0 (scene):
        centers = [CCOMB(S1(UKPOL(obj))) for obj in scene]
        scalings = len(centers) * [S([1,2,3])([sx,sy,sz])]
        scaledCenters = [UK(APPLY(pair)) for pair in
                         zip(scalings, [MK(p) for p in centers])]
        translVectors = [ VECTDIFF((p,q)) for (p,q) in zip(scaledCenters, centers) ]
        translations = [ T([1,2,3])(v) for v in translVectors ]
        return STRUCT([ APPLY((t,obj)) for (t,obj) in zip(translations,scene) ])
    return explode0

def MKPOLS (model):
    V, FV = model
    pols = [MKPOL([[V[v] for v in f],[range(1,len(f)+1)], None]) for f in FV]
    return pols

def LAR2PLASM (topology):
    return AA(AA(lambda k: k+1))(topology)

def VERTS(geoms):
    return COMP([AA(REVERSE),CART,REVERSE])(geoms)

def VERTEXTRUDE((V,coords)):
    return CAT(AA(COMP([AA(AR),DISTR]))(DISTL([V,coords])))


#------------------------------------------------------------------
#--coo is the standard rep using non-ordered triples of numbers----
#--coo := (row::integer, column::integer, value::float)------------
#------------------------------------------------------------------

def format(cmat,shape="csr"):
    """ Transform from plasm's compressed format (list)
        of triples to scipy.sparse corresponding formats.
        """
    n = len(cmat)
    data = arange(n)
    ij = arange(2*n).reshape(2,n)
    for k,item in enumerate(cmat):
        ij[0][k],ij[1][k],data[k] = item
    return scipy.sparse.coo_matrix((data, ij)).asformat(shape)

#------------------------------------------------------------------
def cooCreateFromBrc(ListOfListOfInt):
    COOm = sorted([[k,col,1] for k,row in enumerate(ListOfListOfInt)
                   for col in row ])
    return COOm

if __name__ == "__main__" and True:
    print "\n>>> cooCreateFromBrc"
    V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
    FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
    EV = [[0,1],[0,3],[1,2],[1,3],[1,4],[2,4],[2,5],[3,4],[4,5]]
    cooFV = cooCreateFromBrc(FV)
    cooEV = cooCreateFromBrc(EV)
    print "\ncooCreateFromBrc(FV) =\n", cooFV
    print "\ncooCreateFromBrc(EV) =\n", cooEV


#------------------------------------------------------------------
def csrCreateFromCoo(COOm):
    CSRm = format(COOm,"csr")
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrCreateFromCoo"
    csrFV = csrCreateFromCoo(cooFV)
    csrEV = csrCreateFromCoo(cooEV)
    print "\ncsr(FV) =\n", repr(csrFV)
    print "\ncsr(EV) =\n", repr(csrEV)



#------------------------------------------------------------------
def csrCreate(BRCm,shape=(0,0)):
    if shape == (0,0):
        return csrCreateFromCoo(cooCreateFromBrc(BRCm))
    else:
        CSRm = scipy.sparse.csr_matrix(shape)
        for i,j,v in cooCreateFromBrc(BRCm):
            CSRm[i,j] = v
        return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrCreateFromCoo"
    V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
    FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
    csrFV = csrCreate(FV)
    print "\ncsrCreate(FV) =\n", csrFV


#------------------------------------------------------------------
def csrGetNumberOfRows(CSRm):
    Int = CSRm.shape[0]
    return Int

if __name__ == "__main__" and True:
    print "\n>>> csrGetNumberOfRows"
    print "\ncsrGetNumberOfRows(csrFV) =", csrGetNumberOfRows(csrFV)
    print "\ncsrGetNumberOfRows(csrEV) =", csrGetNumberOfRows(csrEV)

#------------------------------------------------------------------
def csrGetNumberOfColumns(CSRm):
    Int = CSRm.shape[1]
    return Int

if __name__ == "__main__" and True:
    print "\n>>> csrGetNumberOfColumns"
    print "\ncsrGetNumberOfColumns(csrFV) =", csrGetNumberOfColumns(csrFV)
    print "\ncsrGetNumberOfColumns(csrEV) =", csrGetNumberOfColumns(csrEV)

#------------------------------------------------------------------
def csrToMatrixRepresentation(CSRm):
    nrows = csrGetNumberOfRows(CSRm)
    ncolumns = csrGetNumberOfColumns(CSRm)
    ScipyMat = zeros((nrows,ncolumns),int)
    C = CSRm.tocoo()
    for triple in zip(C.row,C.col,C.data):
        ScipyMat[triple[0],triple[1]] = triple[2]
    return ScipyMat

if __name__ == "__main__" and True:
    print "\n>>> csrToMatrixRepresentation"
    print "\nFV =\n", csrToMatrixRepresentation(csrFV)
    print "\nEV =\n", csrToMatrixRepresentation(csrEV)

#------------------------------------------------------------------
def csrToBrc(CSRm):
    nrows = csrGetNumberOfRows(CSRm)
    C = CSRm.tocoo()
    out = [[] for i in range (nrows)]
    [out[i].append(j) for i,j in zip(C.row,C.col)]
    return out

if __name__ == "__main__" and True:
    print "\n>>> csrToBrc"
    print "\nFV =\n", csrToBrc(csrFV)
    print "\nEV =\n", csrToBrc(csrEV)

#------------------------------------------------------------------
#--matrix utility layer--------------------------------------------
#------------------------------------------------------------------

#------------------------------------------------------------------
def csrIsA(CSRm):
    test = CSRm.check_format(True)
    return test==None

if __name__ == "__main__" and True:
    print "\n>>> csrIsA"
    print "\ncsrIsA(csrFV) =",csrIsA(csrFV)

#------------------------------------------------------------------
def csrGet(CSRm,row,column):
    Num = CSRm[row,column]
    return Num

if __name__ == "__main__" and True:
    print "\n>>> csrGet"
    print "\ncsrGet(csrFV,2,3) =",csrGet(csrFV,2,3)
    print "\ncsrGet(csrFV,1,3) =",csrGet(csrFV,1,3)

#------------------------------------------------------------------
def csrSet(CSRm,row,column,value):
    CSRm[row,column] = value
    return None

if __name__ == "__main__" and True:
    print "\n>>> csrSet"
    csrSet(csrFV,2,3,10)
    print "\ncsrSet(csrFV,2,3,10) =",csrGet(csrFV,2,3)
    csrSet(csrFV,2,3,1)
    print "\ncsrSet(csrFV,2,3,1) =",csrGet(csrFV,2,3)

#------------------------------------------------------------------
def csrAppendByRow(CSRm1,CSRm2):
    CSRm = vstack([CSRm1,CSRm2])
    return CSRm

if __name__ == "__main__" and True:
    
    print "\n>>> csrAppendByRow"
    CSRm = csrAppendByRow(csrFV,csrEV)
    print "\ncsrAppendByRow(csrFV,csrEV) =\n", \
        csrToMatrixRepresentation(CSRm)


#------------------------------------------------------------------
def csrAppendByColumn(CSRm1,CSRm2):
    CSRm = hstack([CSRm1,CSRm2])
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrAppendByColumn"
    CSRm = csrAppendByColumn(csrFV,csrFV)
    print "\ncsrAppendByColumn(csrFV,csrFV) =\n", \
        csrToMatrixRepresentation(CSRm)


#------------------------------------------------------------------
def csrSplitByRow(CSRm,k):
    CSRm1 = CSRm[:k]
    CSRm2 = CSRm[k:]
    return CSRm1,CSRm2

if __name__ == "__main__" and True:
    print "\n>>> csrSplitByRow"
    CSRm1,CSRm2 = csrSplitByRow(csrFV,2)
    print "\ncsrSplitByRow(csrFV,2) ="
    print csrToMatrixRepresentation(CSRm1)
    print csrToMatrixRepresentation(CSRm2)


#------------------------------------------------------------------
def csrSplitByColumn(CSRm,k):
    CSRm1 = CSRm.T[:k]
    CSRm2 = CSRm.T[k:]
    return CSRm1.T,CSRm2.T

if __name__ == "__main__" and True:
    print "\n>>> csrSpliceByColumn"
    CSRm1,CSRm2 = csrSplitByColumn(csrFV,4)
    print "\ncsrSplitByColumn(csrFV,4) ="
    print csrToMatrixRepresentation(CSRm1)
    print csrToMatrixRepresentation(CSRm2)


#------------------------------------------------------------------
#--sparse matrix operations layer----------------------------------
#------------------------------------------------------------------

#------------------------------------------------------------------
def csrTranspose(CSRm):
    CSRm = CSRm.T
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrTranspose"
    CSRm = csrTranspose(csrFV)
    print "\ncsrTranspose(csrFV) =\n", csrToMatrixRepresentation(CSRm)


#------------------------------------------------------------------
def csrProduct(CSRm1,CSRm2):
    CSRm = CSRm1 * CSRm2
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrProduct"
    CSRm = csrProduct(csrFV, csrTranspose(csrEV))
    print "\ncsrFE =\n", csrToMatrixRepresentation(CSRm)

#------------------------------------------------------------------
def csrMaxFilter(CSRm):
    # can be done in parallel (by rows)
    nrows = csrGetNumberOfRows(CSRm)
    maxs = [max(CSRm[k].data) for k in range(nrows)]
    coo = CSRm.tocoo()
    triples = [[row,col,1] for row,col,val in zip(coo.row,coo.col,coo.data)
               if maxs[row]==val]
    CSRm = csrCreateFromCoo(triples)
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrMaxFilter"
    CSRm = csrMaxFilter(csrProduct(csrFV, csrTranspose(csrEV)).T).T
    print "\ncsrMaxFilter(csrFE) =\n", csrToMatrixRepresentation(CSRm)

#------------------------------------------------------------------
def csrBinFilter(CSRm):
    # can be done in parallel (by rows)
    coo = CSRm.tocoo()
    triples = [[row,col,1] for row,col,val in zip(coo.row,coo.col,coo.data)
               if val % 2 == 1]
    i, j, data = TRANS(triples)
    CSRm = scipy.sparse.coo_matrix((data,(i,j)),CSRm.shape).tocsr()
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrBinFilter"
    CSRm = csrBinFilter(csrProduct(csrFV, csrTranspose(csrEV)).T).T
    print "\nccsrBinFilter(csrFE) =\n", csrToMatrixRepresentation(CSRm)


#------------------------------------------------------------------
def csrPredFilter(CSRm, pred):
    # can be done in parallel (by rows)
    coo = CSRm.tocoo()
    triples = [[row,col,val] for row,col,val in zip(coo.row,coo.col,coo.data)
               if pred(val)]
    i, j, data = TRANS(triples)
    CSRm = scipy.sparse.coo_matrix((data,(i,j)),CSRm.shape).tocsr()
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> csrPredFilter"
    CSRm = csrPredFilter(csrProduct(csrFV, csrTranspose(csrEV)).T, GE(2)).T
    print "\nccsrPredFilter(csrFE) =\n", csrToMatrixRepresentation(CSRm)


#------------------------------------------------------------------
#--topology interface layer----------------------------------------
#------------------------------------------------------------------

#------------------------------------------------------------------
def csrCreateTotalChain(kn):
    csrMat = csrCreateFromCoo(cooCreateFromBrc(TRANS([kn*[0]])))
    return csrMat

if __name__ == "__main__" and True:
    print "\n>>> csrCreateTotalChain"
    csrMat = csrCreateTotalChain(csrGetNumberOfRows(csrFV))
    print "\ncsrCreateTotalChain(csrGetNumberOfRows(csrFV)) =\n", \
        csrToMatrixRepresentation(csrMat)

#------------------------------------------------------------------
def csrCreateUnitChain(kn,k):
    CSRout = lil_matrix((kn, 1))
    CSRout[k,0] = 1
    return CSRout.tocsr()

if __name__ == "__main__" and True:
    print "\n>>> csrCreateUnitChain"
    CSRm = csrCreateUnitChain(4,2)
    print "\ncsrCreateUnitChain(csrFV,2) =\n", \
        csrToMatrixRepresentation(CSRm)

#------------------------------------------------------------------
def csrExtractAllGenerators(CSRm):
    listOfListOfNumerals = [csrTranspose(CSRm)[k].tocoo().col.tolist()
                            for k in range(CSRm.shape[1])]
    return listOfListOfNumerals

if __name__ == "__main__" and True:
    print "\n>>> csrExtractAllGenerators"
    boundary_2_Op = csrMaxFilter(
                                 csrProduct(csrEV, csrTranspose(csrFV)))
    listOfListOfNumerals = csrExtractAllGenerators(boundary_2_Op)
    print "\ncsrExtractAllGenerators(boundary_2_Op) =\n", \
        listOfListOfNumerals

#------------------------------------------------------------------
def csrChainToCellList(CSRm):
    ListOfInt = CSRm.tocoo().row.tolist()
    return ListOfInt

if __name__ == "__main__" and True:
    print "\n>>> csrChainToCellList"
    total_2_chain = csrCreateTotalChain(csrGetNumberOfRows(csrFV))
    print "\ntotal_2_chain =\n", csrChainToCellList(total_2_chain)
    total_1_chain = csrCreateTotalChain(csrGetNumberOfRows(csrEV))
    print "\ntotal_1_chain =\n", csrChainToCellList(total_1_chain)

#------------------------------------------------------------------
#--topology query layer--------------------------------------------
#------------------------------------------------------------------

#------------------------------------------------------------------
def larCellAdjacencies(CSRm):
    CSRm = csrProduct(CSRm,csrTranspose(CSRm))
    return CSRm

if __name__ == "__main__" and True:
    print "\n>>> larCellAdjacencies"
    adj_2_cells = larCellAdjacencies(csrFV)
    print "\nadj_2_cells =\n", csrToMatrixRepresentation(adj_2_cells)
    adj_1_cells = larCellAdjacencies(csrEV)
    print "\nadj_1_cells =\n", csrToMatrixRepresentation(adj_1_cells)


#------------------------------------------------------------------
def larCellIncidences(CSRm1,CSRm2):
    return csrProduct(CSRm1, csrTranspose(CSRm2))

if __name__ == "__main__" and True:
    print "\n>>> larCellIncidences"
    print "\nlarCellIncidences =\n", csrToMatrixRepresentation(
                                                               larCellIncidences(csrFV,csrEV))


#------------------------------------------------------------------
# FV = d-chain;  EV = (d-1)-chain
def larBoundary(FV,EV):
    csrFV = csrCreateFromCoo(cooCreateFromBrc(FV))
    csrEV = csrCreateFromCoo(cooCreateFromBrc(EV))
    csrBoundary_2 = csrMaxFilter(larCellIncidences(csrEV,csrFV))
    return csrBoundary_2

if __name__ == "__main__" and True:
    print "\n>>> larBoundary"
    csrBoundary_2 = larBoundary(FV,EV)
    print "\ncsrBoundary.T =\n", csrToMatrixRepresentation(csrBoundary_2.T)
    print "\ncsrBoundary_2.T =\n", csrToBrc(csrBoundary_2.T)
    print "\ncsrcoBoundary_1.T =\n", csrToBrc(csrBoundary_2)

#------------------------------------------------------------------
def larBoundaryChain(csrBoundaryMat,brcCellList):
    n = csrGetNumberOfColumns(csrBoundaryMat)
    csrChain = sum([csrCreateUnitChain(n,k) for k in brcCellList])
    print "\nchain =", csrToMatrixRepresentation(csrChain)
    return csrBinFilter(csrProduct(csrBoundaryMat,csrChain))

if __name__ == "__main__" and True:
    print "\n>>> larBoundaryChain"
    csrBoundary_2 = larBoundary(FV,EV)
    chain_1 = larBoundaryChain(csrBoundary_2,[0,2])
    print "\nlarBoundaryChain =\n", csrChainToCellList(chain_1)
    chain_1 = larBoundaryChain(csrBoundary_2,[0])
    print "\nlarBoundaryChain =\n", csrChainToCellList(chain_1)
    chain_1 = larBoundaryChain(csrBoundary_2,
                               range(csrGetNumberOfColumns(csrBoundary_2)))
    print "\nlarBoundaryChain =\n", csrChainToCellList(chain_1)

#------------------------------------------------------------------
def larCoboundaryChain(csrCoBoundaryMat,brcCellList):
    m = csrGetNumberOfColumns(csrCoBoundaryMat)
    csrChain = sum([csrCreateUnitChain(m,k) for k in brcCellList])
    print "\nchain =", csrToMatrixRepresentation(csrChain)
    return csrBinFilter(csrProduct(csrCoBoundaryMat,csrChain))

if __name__ == "__main__" and True:
    print "\n>>> larCoboundaryChain"
    csrCoBoundary_1 = larBoundary(FV,EV).T
    chain_2 = larCoboundaryChain(csrCoBoundary_1,[0,2])
    print "\nlarCoBoundaryChain =\n", csrChainToCellList(chain_2)
    chain_2 = larCoboundaryChain(csrCoBoundary_1,[0])
    print "\nlarCoBoundaryChain =\n", csrChainToCellList(chain_2)
    chain_2 = larCoboundaryChain(csrCoBoundary_1,
                                 range(csrGetNumberOfColumns(csrCoBoundary_1)))
    print "\nlarCoBoundaryChain =\n", csrChainToCellList(chain_2)

#------------------------------------------------------------------
#--model geometry layer--------------------------------------------
#--larOp : model -> model------------------------------------------
#------------------------------------------------------------------
# model = (vertices, topology)
#------------------------------------------------------------------
# binary product of cell complexes

def larProduct(models):
    model1,model2 = models
    V, cells1 = model1
    W, cells2 = model2
    verts = collections.OrderedDict(); k = 0
    for v in V:
        for w in W:
            vertex = tuple(v+w)
            if not verts.has_key(vertex):
                verts[vertex] = k
                k += 1
    cells = [sorted([verts[tuple(V[v]+W[w])] for v in c1 for w in c2])
             for c1 in cells1 for c2 in cells2]
    model = AA(list)(verts.keys()), sorted(cells)
    return model

if __name__ == "__main__" and True:
    geom_0,topol_0 = [[0.],[1.],[2.],[3.],[4.]],[[0,1],[1,2],[3,4]]
    geom_1,topol_1 = [[0.],[1.],[2.]], [[0,1],[1,2]]
    mod_0 = (geom_0,topol_0)
    mod_1 = (geom_1,topol_1)
    
    squares = larProduct([mod_0,mod_1])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
    
    cubes = INSL(larProduct)([mod_0,mod_1,mod_0]) # ==
    cubes = larProduct([squares,mod_0])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
    
    line = [[0.],[1.]], [[0,1]]
    cube = INSL(larProduct)([line,line,line])
    geometry, topology = cube
    VIEW(MKPOL([geometry, LAR2PLASM(topology), None]))
    
    geom_2 = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
    topol_2 = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
    segments = mod_1
    triangles = (geom_2,topol_2)
    wedges = larProduct([triangles,segments])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(wedges)))

#------------------------------------------------------------------
# extrusion of simplicial complexes
# combinatorial algorithm

def cumsum(iterable):
    # cumulative addition: list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable:
        s = s + c
        yield s

def larExtrude(model,pattern):
    V,FV = model
    d = len(FV[0])
    offset = len(V)
    m = len(pattern)
    outcells = []
    for cell in FV:
        # create the indices of vertices in the cell "tube"
        tube = [v + k*offset for k in range(m+1) for v in cell]
        # take groups of d+1 elements, via shifting by one
        rangelimit = len(tube)-d
        cellTube = [tube[k:k+d+1] for k in range(rangelimit)]
        outcells += [scipy.reshape(cellTube,newshape=(m,d,d+1)).tolist()]
    outcells = AA(CAT)(TRANS(outcells))
    outcells = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    outVerts = VERTEXTRUDE((V,coords))
    newModel = outVerts, CAT(outcells)
    return newModel

if __name__ == "__main__" and True:
    V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
    FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8]]
    model = larExtrude((V,FV),2*[1,2,-3])
    VIEW(EXPLODE(1,1,1.2)(MKPOLS(model)))
    
    V0 = [[]]
    CV0 = [[0]]
    model = larExtrude((V0,CV0),4*[1,-1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,4*[1,-1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,4*[1,-1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))

#------------------------------------------------------------------
# extraction of facets of a cell complex

def setup(model,dim):
    V, cells = model
    csr = csrCreate(cells)
    csrAdjSquareMat = larCellAdjacencies(csr)
    csrAdjSquareMat = csrPredFilter(csrAdjSquareMat, GE(dim)) # ? HOWTODO ?
    facets = []
    pointsOfFirstCell = [V[k] for k in cells[0]]
    return V,cells,csr,csrAdjSquareMat,facets

def larFacets(model,dim=3):
    V,cells,csr,csrAdjSquareMat,facets = setup(model,dim)
    print "\ncsrAdjSquareMat =\n",csrToMatrixRepresentation(csrAdjSquareMat)
    # for each input cell i
    cellFacets = []
    for i in range(len(cells)):
        adjCells = csrAdjSquareMat[i].tocoo()
        cell1 = csr[i].tocoo().col
        pairs = zip(adjCells.col,adjCells.data)
        for j,v in pairs:
            if (i!=j):
                cell2 = csr[j].tocoo().col
                cell = list(set(cell1).intersection(cell2))
                cellFacets.append(sorted(cell))
    print "\ncellFacets =",cellFacets
    return V,cellFacets


if __name__ == "__main__" and True:
    # 1D & 2D & 3D CUBOIDAL example: OK
    geom_0,topol_0 = [[0.],[1.],[2.],[3.],[4.]],[[0,1],[1,2],[2,3]]
    geom_1,topol_1 = [[0.],[1.],[2.]], [[0,1],[1,2]]
    mod_0 = (geom_0,topol_0)
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_0,mod_1])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
    skel = larFacets(squares,2)
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
    cubes = INSL(larProduct)([mod_0,mod_1,mod_0]) # ==
    cubes = larProduct([squares,mod_0])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
    skel = larFacets(cubes)
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
#    skel = larFacets(skel,2)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))

#------------------------------------------------------------------


#def prepKey (args): return "["+", ".join(args)+"]"
#
#def fixedPrec(value):
#	if value>0: ret = '%0.1f'% (int(math.floor(value*1E1))/1E1)
#	else: ret = '%0.1f'% -(int(math.floor(-value*1E1))/1E1)
#	if ret == '-0.0': ret = '0.0'
#	return ret
#
#def code (vect): return prepKey(AA(COMP([str,fixedPrec]))(vect))
#
#def mappingOracle(cellPoints):
#    print "\ncellPoints =",cellPoints
#    homPoints = mat([AR([p,1.0]) for p in cellPoints])
#    print "\nhomPoints =\n",homPoints
#    singularValues = eval(code(scipy.linalg.svdvals(homPoints)))
#    print "\nsingularValues =",singularValues
#    rank = len([val for val in singularValues if val!=0.0])
#    print "\nrank =",rank
#    codimension = homPoints.shape[1]-rank
#    print "\ncodimension =",codimension
#    if codimension==0: return False
#    elif codimension==1: return True
#    else: raise Exception('codimension error')
#
#def setup(model,dim):
#    V, cells = model
#    csr = csrCreate(cells)
#    csrAdjSquareMat = larCellAdjacencies(csr)
#    csrAdjSquareMat = csrPredFilter(csrAdjSquareMat, GE(dim)) # ? HOWTODO ?
#    facets = []
#    pointsOfFirstCell = [V[k] for k in cells[0]]
#    isFacet = mappingOracle(pointsOfFirstCell)
#    return V,cells,csr,csrAdjSquareMat,facets,isFacet
#
## predicate to test if a boundary representation is closed
#def isCellClosed(cellFacets):
#    # works for any type of cell (to proof)
#    facets = AA(AA(lambda k: k+1))(cellFacets)
#    # return False
#    if facets != []:
#        theSum = sum(set(CAT(facets)))
#        return (3*theSum == sum(CAT(facets))) #or (2*theSum == sum(CAT(facets)))
#    # ? HOWTODO ?
#    else: return False
#
#def cellBoundaryFacets(points):
#    print "\npoints =",points
#    triangles = Delaunay(points).convex_hull
#    print "\n>>>triangles =",triangles
#    if len(points[0])==2: return triangles
#    centroid = centre(points)
#    print "\n>>>centroid =",centroid
#    points = [VECTDIFF([p,centroid]) for p in points]
#    print "\n>>>points =",points
#    tuples = [[points[k] for k in t] for t in triangles]
#    print "\n>>>tuples =",tuples
#    triangles = [tria if det(t) >= 0 else [tria[1]]+[tria[0]]+tria[2:]
#                 for t,tria in zip(tuples,triangles.tolist())]
#    print "\n>>>triangles =",triangles
#    tuples = [[points[k] for k in t] for t in triangles]
#    print "\n>>>tuples =",tuples
#    normals = [ VECTPROD([ VECTDIFF([p,triple[0]]) for p in triple[1:] ]) for triple in tuples ]
#    print "\n>>>normals =",normals
#    normals = [code(UNITVECT(normal)) for normal in normals]
#    print "\n>>>normals =",normals
#    
#    facets = collections.OrderedDict()
#    print "\n>>>facets =",facets
#    for triangle,normal in TRANS([triangles,normals]):
#        if facets.has_key(normal):
#            facet = facets[normal]
#            facets[normal] = facet.union(triangle)
#        else:
#            facets[normal] = set(triangle)
#    print "\n>>>facets =",facets
#    return AA(list)(facets.values())
#
#def centre(points):
#    n = float(len(points))
#    return AA(lambda coord: sum(coord)/n)(TRANS(points))
#
#def mapping(cellPoints):
#    # TODO
#    m = len(cellPoints[0])
#    vectors = mat([VECTDIFF([p,cellPoints[0]]) for p in cellPoints[1:m]])
#    vectors = vectors.I.tolist()
#    vectors = AA(eval)(AA(code)(vectors))
#    return vectors
#
#def mapping(points):
#    # TODO
#    m = len(points[0])
#    print "\nm =",m
#    print "\npoints =",points
#    matrix = mat([VECTDIFF([p,points[0]]) for p in points[1:m+1]])
#    print "\nmatrix =",matrix
#    vectors = matrix
#    if det(vectors) != 0.0:
#        vects = mat([VECTDIFF([p,points[0]]) for p in points[1:]])
#        print "\nvects =",vects
#        vectors = AA(eval)(AA(code)(matrix * vects.T))
#        print "\nvectors =",vectors
#    return vectors
#
#
#def larFacets(model,dim=3):
#    V,cells,csr,csrAdjSquareMat,facets,isFacet = setup(model,dim)
#    print "\ncsrAdjSquareMat =\n",csrToMatrixRepresentation(csrAdjSquareMat)
#    print "\nisFacet =",isFacet
#    # for each input cell i
#    for i in range(len(cells)):
#        cellFacets = []
#        adjCells = csrAdjSquareMat[i].tocoo()
#        cell1 = csr[i].tocoo().col
#        pairs = zip(adjCells.col,adjCells.data)
#        for j,v in pairs:
#            if (i!=j):
#                cell2 = csr[j].tocoo().col
#                cell = list(set(cell1).intersection(cell2))
#                cellFacets.append(sorted(cell))
#                cellFacets.append(sorted(list(set(cell1).difference(cell))))
#                cellFacets.append(sorted(list(set(cell2).difference(cell))))
#        cellFacets = [facet for facet in cellFacets
#                      if set(facet).issubset(cells[i]) and len(facet)>=dim ] # 3? ... !!!
#        print "\ncellFacets =",cellFacets
#        isClosed = isCellClosed(cellFacets)
#        print "\nisClosed =",isClosed
#        if not isClosed:
#            cellPoints = [V[k] for k in cells[i]]
#            if isFacet: cellPoints = mapping(cellPoints)
#            #facetAddrs = Delaunay(cellPoints).convex_hull
#            facetAddrs = cellBoundaryFacets(cellPoints)
#            print "\n<<<>>>facetAddrs =\n", facetAddrs
#            simplices = [sorted([cells[i][item] for item in pair])
#                         for pair in facetAddrs]
#            simplexLength = len(simplices[0])
#            for simplex in simplices:
#                subset = False
#                for cell in cellFacets:
#                    if set(simplex).issubset(cell): subset = True
#                if not subset: cellFacets += [simplex]
#        facets += cellFacets
#        facets = sorted(facets)
#        facets = [facet for k,facet in enumerate(facets[:-1])if facet!=facets[k+1]] + [facets[-1]]
#    return V,[face for face in facets]# if len(face) >= simplexLength]
#
#
#if __name__ == "__main__" and True:
#    print "\n>>> larFacetComputing"
#    
#    # 2D & 3D SIMPLICIAL example: OK
#    V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
#    FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8]]
#    VIEW(MKPOL([V, AA(AA(lambda k: k+1))(FV), None]))
#    model = larFacets((V,FV))
#    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
#    model = larExtrude((V,FV),1*[1])
#    VIEW(MKPOL([model[0], AA(AA(lambda k: k+1))(model[1]), None]))
#    model = larFacets(model)
#    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
##    model = larFacets(model,dim=2)
##    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
#    
#    # 2D & 3D SIMPLICIAL example: OK
#    V0 = [[]]
#    CV0 = [[0]]
#    model = larExtrude((V0,CV0),2*[1,-1,1])
#    model = larExtrude(model,2*[1,-1,1])
#    skel = larFacets(model)
#    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(skel)))
#    model = larExtrude(model,1*[1,-1,1])
#    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
#    skel = larFacets(model)
#    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(skel)))
##    skel = larFacets(skel,2)
##    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(skel)))
#    
#    # 1D & 2D & 3D CUBOIDAL example: OK
#    geom_0,topol_0 = [[0.],[1.],[2.],[3.],[4.]],[[0,1],[1,2],[3,4]]
#    geom_1,topol_1 = [[0.],[1.],[2.]], [[0,1],[1,2]]
#    mod_0 = (geom_0,topol_0)
#    mod_1 = (geom_1,topol_1)
#    squares = larProduct([mod_0,mod_1])
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
#    skel = larFacets(squares,2)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
#    cubes = INSL(larProduct)([mod_0,mod_1,mod_0]) # ==
#    cubes = larProduct([squares,mod_0])
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
#    skel = larFacets(cubes)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
##    skel = larFacets(skel,2)
##    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
#    
#    
#    # 1D & 2D & 3D CUBOIDAL example: OK
#    V = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]
#    FV = [[0,1,2,3,4,5,6,7]]
#    cube = (V,FV)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cube)))
#    skel = larFacets(cube)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
##    skel = larFacets(skel,2)
##    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
#    
#    
#    # 3D MIXED (general POLYTOPAL) example:  error to check !!
#    geom_2 = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
#    topol_2 = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
#    segments = mod_1
#    triangles = (geom_2,topol_2)
#    wedges = larProduct([triangles,segments])
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(wedges)))
#    skel = larFacets(wedges)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
##    skel = larFacets(skel,2)                           # => KO !!  (troppe 1-facce)
##    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel)))
#    
#    
#    
#    # 3D MIXED (general POLYTOPAL) example:  error to check !!
#    geom = [[0, 0], [1, 0], [0, 1], [1, 1]]
#    topol = [[0, 1, 2], [1, 2, 3]]
#    segment = [[0.],[1.]], [[0,1]]
#    triangle = (geom,topol)
#    wedge = larProduct([triangle,segment])
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(wedge)))
#    skel2D = larFacets(wedge)
#    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel2D)))
##    skel1D = larFacets(skel2D,2)                       # => KO !!  (troppe 1-facce)
##    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(skel1D)))
#
#

#------------------------------------------------------------------
#--application layer (demo)----------------------------------------
#------------------------------------------------------------------

if __name__ == "__main__" and True:
    
    # input of topology and geometry
    V = [[5.,29.],[17.,29.],[8.,25.],[11.,25.],[14.,25.],
         [ 0.,23.], [5.,23.],[17.,23.],[27.,23.],[0.,20.],
         [5.,20.],[8.,19.],[11.,19.],[11.,17.],[14.,17.],[0.,16.],
         [5.,16.],[14.,16.],[17.,16.],[23.,16.],[0.,10.],
         [14.,10.],[23.,10.],[27.,10.],[0.,6.],[5.,6.],[5.,3.],
         [20.,3.],[23.,3.],[20.,0.],[27.,0.]]
    
    EV = [[5,6,0,1,7],[7,18],[18,17],[17,21,20],
          [20,15,16,10,9,5], [12,11,2,3],[3,4,14,13,12],[7,8,23],[22,23],
          [22,19,18],[22,28,27],[26,27],[26,25,24,20],[23,30,29,27]]
    
    FV = [[0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,20,21],
          [7,8,18,19,22,23],[17,18,19,20,21,22,24,25,26,27,28],
          [22,23,27,28,29,30]]
    
    print V
    
    # characteristic matrices
    cooFV = cooCreateFromBrc(FV)
    cooEV = cooCreateFromBrc(EV)
    print "\ncoo(FV) =\n", cooFV
    print "\ncoo(EV) =\n", cooEV
    
    
    csrFV = csrCreateFromCoo(cooFV)
    csrEV = csrCreateFromCoo(cooEV)
    print "\ncsr(FV) =\n", csrFV
    print "\ncsr(EV) =\n", csrEV
    
    
    print "\nFV =\n", csrToMatrixRepresentation(csrFV)
    print "\nEV =\n", csrToMatrixRepresentation(csrEV)
    
    # transposition
    
    csrVF = csrTranspose(csrFV)
    print "\nFV.T =\n", csrToMatrixRepresentation(csrVF)
    
    # product
    csrEF = csrProduct(csrEV, csrVF)
    print "\nEF =\n", csrToMatrixRepresentation(csrEF)
    
    
    # product and transposition
    csrFE = csrTranspose(csrEF)
    print "\nFE =\n", csrToMatrixRepresentation(csrFE)
    
    
    # boundary and coboundary operators
    boundary = csrMaxFilter(csrEF)
    coboundary = csrTranspose(boundary)
    print "\ncoboundary =\n", csrToMatrixRepresentation(coboundary)
    
    
    # boundary 1-chains (edge numerals) of unit 2-chains (i.e. 2-cells)
    boundary_2_Op = boundary
    _2cells = csrExtractAllGenerators(boundary_2_Op)
    print "\n_2cells =\n", _2cells
    
    
    # boundary 0-chains (vertex numerals) of unit 1-chains (i.e. 1-cells)
    boundary_1_Op = csrTranspose(csrEV)
    _1cells = csrExtractAllGenerators(boundary_1_Op)
    print "\n_1cells =\n", _1cells
    
    
    # 2D polygon complex
    polygons2D = [STRUCT([ POLYLINE([V[v]
                                     for v in EV[edge]]) for edge in cell])
                  for cell in _2cells]
    VIEW(EXPLODE(1.2,1.2,1)(polygons2D))
    
    
    # 2D curved cell complex
    cells2D = [SOLIDIFY(STRUCT([ bezier([V[v] for v in EV[edge]])
                                for edge in cell]))
               for cell in _2cells]
    
    
    VIEW(EXPLODE(1.2,1.2,1)(cells2D))
    colors = [RED,GREEN,BLUE,YELLOW,CYAN,MAGENTA,WHITE,GRAY,BLACK]
    VIEW(STRUCT([COLOR(colors[k % 9])(cell) for k,cell in enumerate(cells2D)]))
    
    
    # csr column representation of the total (row-)chain
    total_2_chain = csrCreateTotalChain(csrGetNumberOfRows(csrFV))
    print "\ntotal_2_chain =\n", csrChainToCellList(total_2_chain)
    
    
    # boundary 1-chain computation
    boundary_1_chain = csrBinFilter( csrProduct(boundary, total_2_chain) )
    boundary_1_cells = csrChainToCellList( boundary_1_chain )
    print "\nboundary_1_cells =\n",boundary_1_cells
    
    
    # boundary 1-chain visualization
    boundary1D = AA(bezier)([[V[v] for v in EV[e]]
                             for e in boundary_1_cells ])
    VIEW(STRUCT(boundary1D))
    
    
    # computation of interior 1-cells
    interior_1_cells = [k for k in range(len(EV))
                        if k not in boundary_1_cells]
    print "\ninterior_1_cells =\n", interior_1_cells
    
    
    
    # visualization of interior 1-cells
    interior1D= AA(bezier)([[V[v] for v in EV[e]]
                            for e in interior_1_cells])
    VIEW(STRUCT(AA(COLOR(RED))(boundary1D) +
                cells2D + AA(COLOR(GREEN))(interior1D)))
    
    
    # computation of a 0-complex
    V_0 = [[0.],[5.],[10.],[15.],[20.]]
    VV_0 = AA(LIST)(range(len(V_0))) # [[0],[1],[2],[3],[4]]
    cells0D = AA(MK)([V_0[v[0]] for v in VV_0])
    floors2D = AA(PROD)(CART([cells2D,cells0D]))
    VIEW(EXPLODE(1.2,1.2,1.5)(floors2D))
    
    
    # computation of a 1-complex
    EV_1 = [[0,1],[1,2],[2,3],[3,4]]
    csrEV_1 = csrCreateFromCoo(cooCreateFromBrc(EV_1))
    csrVE_1 = csrTranspose(csrEV_1)
    boundary1 = csrVE_1   # by def: max(VE_1[K]) == 1.
    print "\nboundary1 =\n", csrToMatrixRepresentation(boundary1)
    
    
    # 1D cell complex
    total_1_chain = csrCreateTotalChain(csrGetNumberOfRows(csrEV_1))
    print "\ntotal_1_chain =\n", total_1_chain
    cells1D = [ POLYLINE([ V_0[v] for v in edge ]) for edge in EV_1 ]
    VIEW(EXPLODE(2,2,2)(cells1D))
    
    
    # Cartesian product complexes
    boundary2D=AA(PROD)(CART([boundary1D,cells1D]))
    interior2D=AA(PROD)(CART([interior1D,cells1D]))
    VIEW(EXPLODE(1.1,1.1,1.5)(interior2D + floors2D))
    solid3D = AA(PROD)(CART([cells2D,cells1D]))
    VIEW(EXPLODE(1.1,1.1,1.5)(solid3D))
    boundary3D = boundary2D + floors2D
    VIEW(EXPLODE(1.1,1.1,1.5)(boundary3D))



"""
    Data: Md, (p, q)
    Result: Mp,q
    if p + q = 0 then k ← 1;
    else if p ∗ q = 0
    then k ← p + q;
    else k ← min(p,q);
    for d > h ≥ k do
    Mh ← φ(Mh+1);
    end
    if p+q = 0
    then return M1tM1;
    else if q = 0 then return Mp;
    else if p = 0 then return Mqt;
    else return MpMqt
    """
