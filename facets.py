# algoritmo estrazione faccette
# =============================

from lar import *

# input:  matrice caratteristica $M_d$;
# output:  matrice caratteristica $M_{d-1}$.


def larFacets(model,dim=3):
    """
    Estraction of (d-1)-cellFacets from model := (V,d-cells)
    Return (V, (d-1)-cellFacets)
    """
    V,cells,csr,csrAdjSquareMat,facets = setup(model,dim)
    cellFacets = []
    # for each input cell i
    for i in range(len(cells)):
        adjCells = csrAdjSquareMat[i].tocoo()
        cell1 = csr[i].tocoo().col
        pairs = zip(adjCells.col,adjCells.data)
        for j,v in pairs:
            if (i!=j):
                cell2 = csr[j].tocoo().col
                cell = list(set(cell1).intersection(cell2))
                cellFacets.append(sorted(cell))
    # sort and remove duplicates
    cellFacets = sorted(cellFacets)
    cellFacets = [facet for k,facet in enumerate(cellFacets[:-1]) if facet != cellFacets[k+1]] + [cellFacets[-1]]
    return V,cellFacets

def larSkeletons (model,dim=3):
    """
    Estraction of all skeletons from model := (V,d-cells)
    Return (V, [d-cells, (d-1)-cells, ..., 1-cells]) where p-cells is a list_of_lists_of_integers
    """
    faces = []
    faces.append(model[1])
    for p in range(dim,0,-1):
        model = larFacets(model,dim=p)
        faces.append(model[1])
    return model[0], REVERSE(faces)


"""
Predicate tu test wheather 
"""
def gridTest(bounds, relativeCoords=True):
    if relativeCoords: 
        bounds = [ [0] + PROGRESSIVESUM([0]+bound) for bound in list(abs(scipy.array(bounds))) ]
        print "\nbounds =",bounds
    def multiBoxBoundaryTest0(point):
        return OR(AA(EQ)(CAT(AA(DISTR)(TRANS([bounds,point])))))
    return multiBoxBoundaryTest0

multiBoxBoundaryTest([[1,1,1],[1,1,1]])([1,1])
  
  

def test (bounds):
    """
    Look whether v is on the boundary of a unconnected (multi-dim) interval [vmin_0,vmax_0, ... ,vmin_n,vmax_n]
    Return a Boolean value
    """
    def test0 (v):
        return OR(AA(EQ)(CAT(AA(TRANS)(DISTR([bounds,v])))))
    return test0


if __name__=="__main__":

    # test data (3D cuboidal complex)
    #geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
    interval = range(3+1)
    geom_1,topol_1 = AA(LIST)(interval), [[k,k+1] for k,v in enumerate(interval[:-1])]
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_1,mod_1])
    cubes = larProduct([squares,mod_1])
    #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
    # test data (3D cuboidal complex)
    V = VERTS([interval,interval,interval])
    outverts = [ k for k,v in enumerate(V) if test([V[0],V[-1]])(V[k]) ]
    F3V = list(cubes[1]+[outverts])


    model = (V,F3V)
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    V = model[0]
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F3V[:-1])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V)) ))


    # Testing sparse-matrix representation
    csrF3V = csrCreate(F3V, shape=(len(F3V),len(V)))
    csrF2V = csrCreate(F2V, shape=(len(F2V),len(V)))
    csrF1V = csrCreate(F1V, shape=(len(F1V),len(V)))
    csrF0V = csrCreate(F0V, shape=(len(F0V),len(V)))
    print "\nrepr(csrF3V) =",repr(csrF3V)
    print "\nrepr(csrF2V) =",repr(csrF2V)
    print "\nrepr(csrF1V) =",repr(csrF1V)
    print "\nrepr(csrF0V) =",repr(csrF0V)




    edges_2v = []
    edges_3v = []
    verts = set(CAT(faces[0])) # boundaries of edges
    for edge in faces[1]:
        if len(edge) == 2: edges_2v.append(edge)
        elif len(edge) == 3:
            if edge[1] not in verts: edges_3v.append(edge)
            elif edge[0] not in verts: edges_3v.append([edge[1],edge[0],edge[2]])
            else: edges_3v.append([edge[1],edge[2],edge[0]])
    edges = CAT([edges_2v,edges_3v])
    VIEW(STRUCT(AA(bezier)([[V[v] for v in edge] for edge in edges])))

"""


EV = edges
FV = faces[2]
csrFV = csrCreate(FV)
csrEV = csrCreate(EV,shape=(84,64))
csrVF = csrTranspose(csrFV)
csrEF = csrProduct(csrEV, csrVF)
boundary = csrMaxFilter(csrEF)
# 2D curved cell complex
_2cells = csrExtractAllGenerators(boundary)


def cellMapping(cell):
    return BEZIER(S2)([ BEZIER(S1)([V[v] for v in EV[edge]]) for edge in cell])


dom2D = PROD([INTERVALS(1)(10),INTERVALS(1)(10)])
cells2D = [MAP(cellMapping(cell))(dom2D) for cell in _2cells]
#cells2D = [SOLIDIFY(STRUCT([ bezier([V[v] for v in EV[edge]]) for edge in cell])) for cell in _2cells]
VIEW(EXPLODE(2,2,2)(cells2D))


##############################


V0 = [[]]
CV0 = [[0]]
model = larExtrude((V0,CV0),2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
model = larExtrude(model,2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))

model,faces = larSkeletons(model,dim=2)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((model[0],faces[2]))))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((model[0],faces[1]))))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((model[0],faces[0]))))


##############################


V0 = [[]]
CV0 = [[0]]
model = larExtrude((V0,CV0),2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
model = larExtrude(model,2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
#model = larExtrude(model,2*[1,1,1])
#VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))


V = model[0]
outverts = [ k for k,v in enumerate(V) if test([V[0],V[-1]])(V[k]) ]
EV = list(model[1]+[outverts])

# 2. Calcolo degli spigoli delle faccette

edges = larFacets((V,EV),dim=2)
VIEW(EXPLODE(2,2,2)(MKPOLS(edges)))

# 3. Calcolo dei vertici di bordo degli spigoli delle faccette

vertices = larFacets((V,edges[1]),dim=1)
VIEW(EXPLODE(2,2,2)(MKPOLS(vertices)))

# 4. ordering of facet edges

edges_2v = []
edges_3v = []
verts = set(CAT(vertices[1])) # boundaries of edges
for edge in edges[1]:
    if len(edge) == 2: edges_2v.append(edge)
    elif len(edge) == 3:
        if edge[1] not in verts: edges_3v.append(edge)
        elif edge[0] not in verts: edges_3v.append([edge[1],edge[0],edge[2]])
        else: edges_3v.append([edge[1],edge[2],edge[0]])
edges = CAT([edges_2v,edges_3v])
VIEW(STRUCT(AA(bezier)([[V[v] for v in edge] for edge in edges])))



##############################


V0 = [[]]
CV0 = [[0]]
model = larExtrude((V0,CV0),2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
model = larExtrude(model,2*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
model = larExtrude(model,1*[1,1,1])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))


V = model[0]
outverts = [ k for k,v in enumerate(V) if test([V[0],V[-1]])(V[k]) ]
FV = list(model[1]+[outverts])


# 1. Calcolo faccette come sottoinsiemi $S_{i,j}$ di vertici comuni alle coppie di celle $S_i, S_j \in \Lambda_d(X)$ ($i,j\in \N$). $\sharp S_{i,j} \geq \dim(S_i) = \dim(S_j) = d$

facets = larFacets((V,FV),dim=3)
VIEW(EXPLODE(2,2,2)(MKPOLS(facets)))

# 2. Calcolo degli spigoli delle faccette

edges = larFacets(facets,dim=2)
VIEW(EXPLODE(2,2,2)(MKPOLS(edges)))

# 3. Calcolo dei vertici di bordo degli spigoli delle faccette

vertices = larFacets((V,edges[1]),dim=1)
VIEW(EXPLODE(2,2,2)(MKPOLS(vertices)))

# 4. ordering of facet edges

edges_2v = []
edges_3v = []
verts = set(CAT(vertices[1])) # boundaries of edges
for edge in edges[1]:
    if len(edge) == 2: edges_2v.append(edge)
    elif len(edge) == 3:
        if edge[1] not in verts: edges_3v.append(edge)
        elif edge[0] not in verts: edges_3v.append([edge[1],edge[0],edge[2]])
        else: edges_3v.append([edge[1],edge[2],edge[0]])
edges = CAT([edges_2v,edges_3v])
VIEW(STRUCT(AA(bezier)([[V[v] for v in edge] for edge in edges])))
"""
