from lar import *


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
    cellFacets = [cell for k,cell in enumerate(sorted(cellFacets)) if k%2==0]
    print "\ncellFacets =",cellFacets
    return V,sorted(cellFacets)


V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
EV = [[0,1],[1,2],[0,3],[1,4],[2,5],[3,4],[4,5],[3,6],[4,7],[5,8],[6,7],[7,8]]
FV = [[0,1,3,4],[1,2,4,5],[3,4,6,7],[4,5,7,8],[0,1,2,3,5,6,7,8]]

VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EV))))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV[:-1]))))

csrEV = csrCreate(EV)
csrFV = csrCreate(FV)

print "\ncsrEV =\n", csrToMatrixRepresentation(csrEV)
print "\ncsrFV =\n", csrToMatrixRepresentation(csrFV)

csrEE = larCellAdjacencies(csrEV)
csrFF = larCellAdjacencies(csrFV)

print "\ncsrEE =\n", csrToMatrixRepresentation(csrEE)
print "\ncsrFF =\n", csrToMatrixRepresentation(csrFF)

facets = larFacets((V,FV),dim=2)
#EV = facets[1]
EV = [[1, 0, 3], [1, 2, 5], [1, 4], [3, 4], [3, 6, 7], [4, 5], [4, 7], [5, 8, 7]]
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(facets)))

skel1D = AA(bezier)([[V[v] for v in EV[e]] for e,ev in enumerate(EV) ])
VIEW(STRUCT(skel1D))
skel1D = AA(POLYLINE)([[V[v] for v in EV[e]] for e,ev in enumerate(EV) ])
VIEW(STRUCT(skel1D))

csrVV = csrProduct(csrTranspose(csrFV),csrFV)
print "\ncsrVV =\n", csrToMatrixRepresentation(csrVV)

csrVV = csrProduct(csrTranspose(csrFV[:-1]),csrFV[:-1])
print "\ncsrVV' =\n", csrToMatrixRepresentation(csrVV)

csrVV = csrProduct(csrTranspose(csrEV),csrEV)
print "\ncsrVV'' =\n", csrToMatrixRepresentation(csrVV)

#print "\nfacets' =\n", facets[1]
#print "\nboundary =\n", csrToMatrixRepresentation(larBoundary(FV,facets[1]))
#facets = [[1, 4], [3, 4], [4, 5], [4, 7], [8,8]]