from lar import *

"""
    Examples of skeleton and facets extraction via topological methods.
    Only grid (cuboidal and simplicial) examples are provided.
    
    TODO: further (general) examples needed: e.g. Voronoi complexes.
"""


if __name__=="__main__":
    
    # -----------------------------------------------------------------
    # 1. rounded cuboidal 3-grid
    
    
    # test data (3D cuboidal complex)
    #geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
    interval = range(3+1)
    geom_1,topol_1 = AA(LIST)(interval), [[k,k+1] for k,v in enumerate(interval[:-1])]
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_1,mod_1])
    cubes = larProduct([squares,mod_1])
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
    # test data (3D cuboidal complex)
    V = VERTS([interval,interval,interval])
    outverts = [[ k for k,v in enumerate(V) if outerVertexTest([V[0],V[-1]])(V[k]) ]]
    F3V = cubes[1]+ outverts
    # OK
    
    model = (V,F3V)
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))
    V,faces = larSkeletons(model,dim=3,grid=True)
    F0V, F1V, F2V, F3V = faces
    print "AA(LEN)([F0V, F1V, F2V, F3V]) =", AA(LEN)([F0V, F1V, F2V, F3V])
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F3V[:-1])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V)) ))
    
    edges_2v = []
    edges_3v = []
    # ordering of vertices of quadratic curves
    verts = set(CAT(faces[0])) # boundaries of edges
    for edge in faces[1]:
        if len(edge) == 2: edges_2v.append(edge)
        elif len(edge) == 3:
            if edge[1] not in verts: edges_3v.append(edge)
            elif edge[0] not in verts: edges_3v.append([edge[1],edge[0],edge[2]])
            else: edges_3v.append([edge[1],edge[2],edge[0]])
    edges = CAT([edges_2v,edges_3v])
    VIEW(STRUCT(AA(bezier)([[V[v] for v in edge] for edge in edges])))
    
    # -----------------------------------------------------------------
    # 2. rounded cuboidal 2-grid
    
    
    # test data (3D cuboidal complex)
    #geom_1,topol_1 = [[0.],[1.],[2.],[3.]], [[0,1],[1,2],[2,3]]
    interval = range(10+1)
    geom_1,topol_1 = AA(LIST)(interval), [[k,k+1] for k,v in enumerate(interval[:-1])]
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_1,mod_1])

    V = VERTS([interval,interval])
    outverts = [[ k for k,v in enumerate(V) if outerVertexTest([V[0],V[-1]])(V[k]) ]]
    F2V = squares[1]+ outverts
    
    model = (V,F2V)
    V,faces = larSkeletons(model,dim=2,grid=True)
    F0V, F1V, F2V= faces
    print "AA(LEN)([F0V, F1V, F2V]) =", AA(LEN)([F0V, F1V, F2V])
    V = model[0]
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V[:-1])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V)) ))
    
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
    


    # -----------------------------------------------------------------
    # 2. cuboidal 3-grid
    
    interval = range(3+1)
    geom_1,topol_1 = AA(LIST)(interval), [[k,k+1] for k,v in enumerate(interval[:-1])]
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_1,mod_1])
    cubes = larProduct([squares,mod_1])
    V = VERTS([interval,interval,interval])
    outverts = boundarGrid(cubes,[0,0,0],[3,3,3])
    F3V = cubes[1]+ outverts
    
    
    model = (V,F3V)
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    print "AA(LEN)([F0V, F1V, F2V, F3V]) =", AA(LEN)([F0V, F1V, F2V, F3V])
    V = model[0]
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F3V[:-6])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))
    
    # Testing sparse-matrix representation
    csrF3V = csrCreate(F3V, shape=(len(F3V),len(V)))
    csrF2V = csrCreate(F2V, shape=(len(F2V),len(V)))
    csrF1V = csrCreate(F1V, shape=(len(F1V),len(V)))
    csrF0V = csrCreate(F0V, shape=(len(F0V),len(V)))
    print "\nrepr(csrF3V) =",repr(csrF3V)
    print "\nrepr(csrF2V) =",repr(csrF2V)
    print "\nrepr(csrF1V) =",repr(csrF1V)
    print "\nrepr(csrF0V) =",repr(csrF0V)


    # -----------------------------------------------------------------
    # 4. cuboidal 2-grid
    
    interval = range(3+1)
    geom_1,topol_1 = AA(LIST)(interval), [[k,k+1] for k,v in enumerate(interval[:-1])]
    mod_1 = (geom_1,topol_1)
    squares = larProduct([mod_1,mod_1])
    V = VERTS([interval,interval])
    outverts = boundarGrid(squares,[0,0],[3,3])
    F2V = squares[1]+ outverts
    
    
    model = (V,F2V)
    V,faces = larSkeletons(model,dim=2)
    F0V, F1V, F2V = faces
    print "AA(LEN)([F0V, F1V, F2V]) =", AA(LEN)([F0V, F1V, F2V])
    V = model[0]
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))
    

    # -----------------------------------------------------------------
    # 5. simplicial 3-grid

    V0 = [[]]
    CV0 = [[0]]
    model = larExtrude((V0,CV0),1*[1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,1*[1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,1*[1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    outverts = boundarGrid(model,[0,0,0],[3,3,3])
    F3V = model[1]+ outverts

    V = model[0]
    model = (V,F3V)
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    print "AA(LEN)([F0V, F1V, F2V, F3V]) =", AA(LEN)([F0V, F1V, F2V, F3V])
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F3V[:-6])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V)) ))
    VIEW(STRUCT( MKPOLS((V,F1V)) ))

    # -----------------------------------------------------------------
    # 6. simplicial 2-grid

    V0 = [[]]
    CV0 = [[0]]
    model = larExtrude((V0,CV0),1*[1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,1*[1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    outverts = boundarGrid(model,[0,0],[3,3])
    F2V = model[1]+ outverts

    V = model[0]
    model = (model[0],F2V)
    V,faces = larSkeletons(model,dim=2)
    F0V, F1V, F2V = faces
    print "AA(LEN)([F0V, F1V, F2V]) =", AA(LEN)([F0V, F1V, F2V])
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(2,2,2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))
    VIEW(STRUCT( MKPOLS((V,F1V)) ))


