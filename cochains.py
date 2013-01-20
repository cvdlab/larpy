from lar import *

# -----------------------------------------------------------------------------
# --- graphical test of the concept -------------------------------------------
# -----------------------------------------------------------------------------

a0 = [[-3,-5,3],[-1,-4,-1],[1,-5,-0.5],[4,-6,4]]
a1 = [[-2,-2,1],[0,-1,0],[2,-2,0],[6,-0.5,-2]]
a2 = [[-2,0,-3],[0,0,0],[2,0,0],[3,1,-3]]

a = [a0,a1,a2]

b0 = [[-2,0,-3],[0,0,0],[2,0,0],[3,1,-3]]
b1 = [[-4,4,-3],[0.5,3,1],[2,2,1],[4,3,-1.5]]
b2 = [[0,6,-1],[1.5,5,0.5],[3,4,2],[4,6,-3]]

b = [b0,b1,b2]

scaffold_a = STRUCT(AA(POLYLINE)(CAT([a,b])))
scaffold_b = STRUCT(AA(POLYLINE)(CAT(AA(TRANS)([a,b]))))
scaffold = STRUCT([scaffold_a,scaffold_b])

VIEW(scaffold)

ca,cb,cta,ctb = AA(AA(BEZIER(S1)))([a,b,TRANS(a),TRANS(b)])
dom1D = INTERVALS(1)(20)

ca0,ca1,ca2,cb0,cb1,cb2 = CONS(AA(MAP)(CAT([ca,cb])))(dom1D)
cgrid_a = STRUCT([ca0,ca1,ca2,cb0,cb1,cb2])

cta0,cta1,cta2,cta3,ctb0,ctb1,ctb2,ctb3 = CONS(AA(MAP)(CAT([cta,ctb])))(dom1D)
cgrid_b = STRUCT([cta0,cta1,cta2,cta3,ctb0,ctb1,ctb2,ctb3])

VIEW(STRUCT([scaffold,cgrid_a,cgrid_b]))

sa = COLOR(RED)(MAP(BEZIER(S2)(ca))(dom2D))
sb = COLOR(GREEN)(MAP(BEZIER(S2)(cb))(dom2D))
dom2D = PROD([dom1D,dom1D])

VIEW(STRUCT([scaffold,cgrid_a,cgrid_b,sa,sb]))




