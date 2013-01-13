from lar import *


def larFacets(model,dim=3):
    V,cells,csr,csrAdjSquareMat,facets = setup(model,dim)
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

# VV computation via ALL the 2-cells (including the environment one)
csrVV = csrProduct(csrTranspose(csrFV),csrFV)
print "\ncsrVV =\n", csrToMatrixRepresentation(csrVV)
## VV computation via the interior 2-cells 
#csrVV = csrProduct(csrTranspose(csrFV[:-1]),csrFV[:-1])
#print "\ncsrVV' =\n", csrToMatrixRepresentation(csrVV)
## VV computation via the 1-cells
#csrVV = csrProduct(csrTranspose(csrEV),csrEV)
#print "\ncsrVV'' =\n", csrToMatrixRepresentation(csrVV)

# computation of "corner" 0-cells, i.e. incident to only one d-cell (d=maxdim)
csrVV = csrProduct(csrTranspose(csrFV),csrFV)
cornerVertices = [k for k in range(csrGetNumberOfRows(csrVV)) if csrVV[k,k]==2]
print "\ncornerVertices =\n", cornerVertices


geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
mod_1 = (geom_1,topol_1)
squares = larProduct([mod_1,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
cubes = INSL(larProduct)([mod_1,mod_1,mod_1]) # ==
cubes = larProduct([squares,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))


def test (vmin, vmax):
    # look whether v is on the boundary of the (multidimensional) interval [vmin,vmax]
	def test0 (v):
		return OR(AA(EQ)(CAT(AA(TRANS)([[vmin,v],[vmax,v]]))))
	return test0


V = VERTS([range(4),range(4),range(4)])
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
FV = list(cubes[1]+[outverts])
csrFV = csrCreate(FV)
csrVF = csrTranspose(csrFV)
csrFF = csrProduct(csrFV,csrVF)
print "\ncsrFF =\n", csrToMatrixRepresentation(csrFF)
facets = larFacets((V,FV),dim=3)
VIEW(EXPLODE(2.5,2.5,2.5)(MKPOLS(facets)))

# computation of "corner" 0-cells, i.e. incident to one/two/four d-cells (d=maxdim)
csrVV = csrProduct(csrTranspose(csrFV[:-1]),csrFV[:-1])
print "\ncsrVV =\n", csrToMatrixRepresentation(csrVV)
vcorners = set([k for k,vert in enumerate(csrVV) if vert[0,k]==1])
print "\ncornerVertices =\n", vcorners
vedges = set([k for k,vert in enumerate(csrVV) if vert[0,k]==2])
print "\nedgeVertices =\n", vedges
vfaces = set([k for k,vert in enumerate(csrVV) if vert[0,k]==4])
print "\nfaceVertices =\n", vfaces

def vertexClassification(vclasses):
    def vertexClassification0(facet):
        return [list(vclass.intersection(facet)) for vclass in vclasses]
    return vertexClassification0

classify3D = vertexClassification([vcorners,vedges,vfaces])
vertexSubsets = [classify3D(facet)  if len(facet)>4 else [facet]  for facet in facets[1]]
print "\nvertexSubsets =\n", vertexSubsets

vertexCells = [set() for v in V]
for k,cell in enumerate(FV[:-1]):
    for v in cell:
        vertexCells[v] = vertexCells[v].union([k])
vertexCellSubsets = [ [[list(vertexCells[v]) for v in sub] for sub in subset] if len(subset)>1 else [subset]
                     for subset in vertexSubsets ]
print "\nvertexCellSubsets =\n", vertexCellSubsets


# combinazion of facet vertices of various groups
from itertools import combinations

def mapVerts2Cells (vertexSubsets,vertexCellSubsets):
    pairs = zip(vertexSubsets,vertexCellSubsets)
    fun = lambda args: CAT(AA(TRANS)(TRANS(args)))
    return dict(CAT(AA(fun)(pairs)))

def combFacetVertices (cverts):
    def combFacetVertices0 (vertGroups):
        print "\nvertGroups =",vertGroups
        if len(vertGroups) == 1: return vertGroups
        else:
            K_1 = [vertGroups[0] + list(pair) for pair in combinations(vertGroups[1],2) ]
            print "\nK_1 =", K_1
            if len(K_1) == 3: K_2 = [sub + [vertGroups[2][k]] for k,sub in enumerate(K_1)]
            elif len(K_1) == 1: K_2 = [K_1[0] + vertGroups[2][:2], K_1[0] + vertGroups[2][2:]]
            else: raise Exception('facet error')
        return K_2
    return combFacetVertices0


cverts = mapVerts2Cells(vertexSubsets,vertexCellSubsets)
print "\ncverts =",cverts
all = CAT([combFacetVertices(cverts)(vertGroups) for vertGroups in vertexSubsets])
print "\nall =", all
VIEW(EXPLODE(2.5,2.5,2.5)(MKPOLS((V,all))))








