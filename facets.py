# algoritmo estrazione faccette
# =============================

from lar import *

# input:  matrice caratteristica $M_d$;
# output:  matrice caratteristica $M_{d-1}$.

# test data (3D cuboidal complex)
geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
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
    print "\ncellFacets =",cellFacets
    return V,cellFacets

# look whether v is on the boundary of the (multidimensional) interval [vmin,vmax]
def test (vmin, vmax):
	def test0 (v):
		return OR(AA(EQ)(CAT(AA(TRANS)([[vmin,v],[vmax,v]]))))
	return test0


# test data (3D cuboidal complex)
V = VERTS([range(4),range(4),range(4)])
outverts = [ k for k,v in enumerate(V) if test(V[0],V[-1])(V[k]) ]
FV = list(cubes[1]+[outverts])


# decomposizione iniziale
# -----------------------

# 1. Calcolo faccette come sottoinsiemi $S_{i,j}$ di vertici comuni alle coppie di celle $S_i, S_j \in \Lambda_d(X)$ ($i,j\in \N$). $\sharp S_{i,j} \geq \dim(S_i) = \dim(S_j) = d$

facets = larFacets((V,FV),dim=3)
VIEW(EXPLODE(2,2,2)(MKPOLS(facets)))

# 2. Calcolo degli spigoli delle faccette e delle faccette di bordo

edges = larFacets((V,facets[1]),dim=2)
VIEW(EXPLODE(2,2,2)(MKPOLS(edges)))





boundaryEdges = larBoundary(facets[1],edges[1])
_2cells = csrExtractAllGenerators(boundaryEdges)
polylines = [ POLYLINE([V[v] for v in edges[1][edge]]) for edge in range(csrGetNumberOfColumns(boundaryEdges))]
VIEW(EXPLODE(2,2,2)(polylines))
polygons2D = [STRUCT([ POLYLINE([V[v] for v in edges[1][edge]]) for edge in cell]) for cell in _2cells]
VIEW(EXPLODE(3,3,3)(polygons2D))


# 2. Ordinamento delle faccette di bordo in funzione crescente del numero di vertici.

vertexOrder,FV = TRANS([pair for pair in sorted([(len(facet),facet) for facet in facets[1]])])
print "\nsortedFacets =",FV
print "\nvertexOrder =",vertexOrder

# 3. Calcolo di un indice "invertito":  v -> incident facets

cellsOnVert = [[] for k in range(len(V))]
[cellsOnVert[v].append(f) for f,facet in enumerate(FV) for v in facet]
print "\ncellsOnVert =",cellsOnVert


# 2. Partizione delle faccette di cardinalita' $> d-1$ in sottoinsiemi "piatti", ovvero con dimensione del guscio affine pari a $d-1$.


# 3. Sia $M_{d-1}$ la matrice caratteristica dei sottoinsiemi $S_{i,j}$, ordinata per righe rispetto al numero (crescente) di vertici incidenti. Si vuole sostituire ad ogni riga con somma maggiore di $d$ un insieme di righe di somma $d$ corrispondenti ad un insieme connesso di faccette piatte con gli stessi vertici, ovvero con la stessa somma sulle colonne.


# 4.  Questa decomposizione delle righe e' unica? probabilmente no. Dovendo scegliere tra piu' decomposizioni ammissibili, sceglieremmo quella (unica?) dove la catena cobordo del bordo abbia---a tratti---lo stesso guscio affine delle faccette esterne adiacenti al bordo.


# decomposizione fine
# -------------------

# 5. Consideriamo la prima faccetta di cardinalita' maggiore di $d$, sia $S_k$, di indice $k$ ed estraiamo la sottomatrice massimale $N_k$ delle righe di $M_{d-1}$ con indice di riga minore di $k$  e somma per colonne pari a $2M_k$, ovvero delle faccette adiacenti di cardinalita' eguale a $d$.

# 6. Partizioniamo $M_k$ (utilizzando $N_k$ e/o $N_k^T.N_k$ ed $N_k.N_k^T$) in piu' righe $R_k$, che sostituiamo ad $M_k$.

# 7. Incrementiamo $k$ e ritorniamo al punto 5, fino a terminare.
