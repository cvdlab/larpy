from lar import *

# -----------------------------------------------------------------------------
# --- graphical test of the concept -------------------------------------------
# -----------------------------------------------------------------------------

V = [[-3,-5,3],[-1,-4,-1],[1,-5,-0.5],[4,-6,4],
     [-2,-2,1],[0,-1,0],[2,-2,0],[6,-0.5,-2],
     [-2,0,-3],[0,0,0],[2,0,0],[3,1,-3],
     [-4,4,-3],[0.5,3,1],[2,2,1],[4,3,-1.5],
     [0,6,-1],[1.5,5,0.5],[3,4,2],[4,6,-3]]

a = [[0,1,2,3],[4,5,6,7],[8,9,10,11]]
b = [[8,9,10,11],[12,13,14,15],[16,17,18,19]]

A,B = [[V[k] for k in list] for list in a], [[V[k] for k in list] for list in b]

scaffold_a = STRUCT(AA(POLYLINE)(CAT([A,B])))
scaffold_b = STRUCT(AA(POLYLINE)(CAT(AA(TRANS)([A,B]))))
scaffold = STRUCT([scaffold_a,scaffold_b])

VIEW(scaffold)

ca,cb,cta,ctb = AA(AA(BEZIER(S1)))([A,B,TRANS(A),TRANS(B)])
dom1D = INTERVALS(1)(20)
dom2D = PROD([dom1D,dom1D])

ca0,ca1,ca2,cb0,cb1,cb2 = CONS(AA(MAP)(CAT([ca,cb])))(dom1D)
cgrid_a = STRUCT([ca0,ca1,ca2,cb0,cb1,cb2])

cta0,cta1,cta2,cta3,ctb0,ctb1,ctb2,ctb3 = CONS(AA(MAP)(CAT([cta,ctb])))(dom1D)
cgrid_b = STRUCT([cta0,cta1,cta2,cta3,ctb0,ctb1,ctb2,ctb3])

VIEW(STRUCT([scaffold,cgrid_a,cgrid_b]))

sa = COLOR(RED)(MAP(BEZIER(S2)(ca))(dom2D))
sb = COLOR(GREEN)(MAP(BEZIER(S2)(cb))(dom2D))

VIEW(STRUCT([scaffold,cgrid_a,cgrid_b,sa,sb]))

# -----------------------------------------------------------------------------
# topology computation

FV = AA(CAT)([a,b])
EV = [a[0],a[-1],b[-1],TRANS(a)[0], TRANS(a)[-1],TRANS(b)[0], TRANS(b)[-1]]

FV = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
      [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
EV = [[0, 1, 2, 3], [8, 9, 10, 11], [16, 17, 18, 19],
      [0, 4, 8], [3, 7, 11], [8, 12, 16], [11, 15, 19]]

csrFV = csrCreate(FV)
csrEV = csrCreate(EV)
csrBoundary_2 = larBoundary(FV,EV)
print "\ncsrBoundary_2.T =\n", csrToMatrixRepresentation(csrBoundary_2.T)
chain_1 = larBoundaryChain(csrBoundary_2, range(csrGetNumberOfColumns(csrBoundary_2)))
print "\nlarBoundaryChain =\n", csrChainToCellList(chain_1)
_1cells = csrExtractAllGenerators(chain_1)[0]

boundaryEV = [[V[v] for v in EV[e]] for e in _1cells]
boundaryMaps = AA(BEZIER(S1))(boundaryEV)
boundary = STRUCT(CONS(AA(MAP)(boundaryMaps))(dom1D))
VIEW(boundary)

