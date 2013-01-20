from lar import *

# -----------------------------------------------------------------------------
# --- graphical test of the concept -------------------------------------------
# -----------------------------------------------------------------------------

a0 = [[-3,-5,3],[-1,-4,-1],[1,-5,-0.5],[4,-6,4]]
a1 = [[-2,-2,1],[0,-1,0],[2,-2,0],[6,-0.5,-2]]
a2 = [[-2,0,-3],[0,0,0],[2,0,0],[3,1,-3]]

a = [a0,a1,a2]
ta = TRANS(a)

b0 = [[-2,0,-3],[0,0,0],[2,0,0],[3,1,-3]]
b1 = [[-4,4,-3],[0.5,3,1],[2,2,1],[4,3,-1.5]]
b2 = [[0,6,-1],[1.5,5,0.5],[3,4,2],[4,6,-3]]

b = [b0,b1,b2]
tb = TRANS(b)

scaffold_a = STRUCT(AA(POLYLINE)(CAT([a,b])))
scaffold_b = STRUCT(AA(POLYLINE)(CAT([ta,tb])))
scaffold = STRUCT([scaffold_a,scaffold_b])

VIEW(scaffold)

ca,cb = AA(AA(BEZIER(S1)))([a,b])
cta,ctb = AA(AA(BEZIER(S1)))([ta,tb])
dom1D = INTERVALS(1)(20)


ca0 = MAP(ca[0])(dom1D)
ca1 = MAP(ca[1])(dom1D)
ca2 = MAP(ca[2])(dom1D)
cb0 = MAP(cb[0])(dom1D)
cb1 = MAP(cb[1])(dom1D)
cb2 = MAP(cb[2])(dom1D)

cgrid_a = STRUCT([ca0,ca1,ca2,cb0,cb1,cb2])

cta0 = MAP(cta[0])(dom1D)
cta1 = MAP(cta[1])(dom1D)
cta2 = MAP(cta[2])(dom1D)
cta3 = MAP(cta[3])(dom1D)
ctb0 = MAP(ctb[0])(dom1D)
ctb1 = MAP(ctb[1])(dom1D)
ctb2 = MAP(ctb[2])(dom1D)
ctb3 = MAP(ctb[3])(dom1D)

cgrid_b = STRUCT([cta0,cta1,cta2,cta3,ctb0,ctb1,ctb2,ctb3])

VIEW(STRUCT([scaffold,cgrid_a,cgrid_b]))





