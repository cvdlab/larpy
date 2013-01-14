# algoritmo estrazione faccette
# =============================

from lar import *

# input:  matrice caratteristica $M_d$;
# output:  matrice caratteristica $M_{d-1}$.

# test data (3D cuboidal complex)
geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
geom_1,topol_1 = AA(LIST)(range(10)), [[k,k+1] for k,v in enumerate(range(10)[:-1])]
mod_1 = (geom_1,topol_1)
squares = larProduct([mod_1,mod_1])
cubes = larProduct([squares,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))

def larFacets(model,dim=3):
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

# look whether v is on the boundary of the (multidimensional) interval [vmin,vmax]
def test (vmin, vmax):
	def test0 (v):
		return OR(AA(EQ)(CAT(AA(TRANS)([[vmin,v],[vmax,v]]))))
	return test0


# test data (3D cuboidal complex)
V = VERTS([range(4),range(4),range(4)])
V = VERTS([range(10),range(10),range(10)])
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
FV = list(cubes[1]+[outverts])


# decomposizione iniziale
# -----------------------

# 1. Calcolo faccette come sottoinsiemi $S_{i,j}$ di vertici comuni alle coppie di celle $S_i, S_j \in \Lambda_d(X)$ ($i,j\in \N$). $\sharp S_{i,j} \geq \dim(S_i) = \dim(S_j) = d$

facets = larFacets((V,FV),dim=3)
VIEW(EXPLODE(2,2,2)(MKPOLS(facets)))

# 2. Calcolo degli spigoli delle faccette 

edges = larFacets((V,facets[1]),dim=2)
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


##################

# test data (2D cuboidal complex)
geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
geom_1,topol_1 = AA(LIST)(range(10)), [[k,k+1] for k,v in enumerate(range(10)[:-1])]
mod_1 = (geom_1,topol_1)
squares = larProduct([mod_1,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))

V = squares[0]
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
EV = list(squares[1]+[outverts])

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
#model = larExtrude(model,2*[1,1,1])
#VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))


V = model[0]
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
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
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
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

