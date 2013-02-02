from lar import *


#-------------------------------------------------------------------
# dimension-independent conversion from discrete image to LAR

def larFromImageBlocks(blocks):
    dim = len(blocks[0][0])
    print "\ndim =\n",dim
    blocks = AA(CAT)(blocks)
    print "\nblocks =\n",blocks
    bounds = TRANS(blocks)
    print "\nbounds =\n",bounds
    cells = [[] for k in range(len(blocks))]
    print "\ncells =\n",cells
    minMaxCoords = zip(AA(min)(bounds[:dim]), AA(max)(bounds[dim:]))
    print "\nminMaxCoords =\n",minMaxCoords
    gridPoints = CART(AA(FROMTO)(AA(list)(minMaxCoords)))
    print "\ngridPoints =\n",gridPoints
    counter,V,i = [],[],0
    for point in gridPoints:
        for k,cell in enumerate(blocks):
            print k,cell
            pmin,pmax = cell[:dim],cell[dim:]
            print "\npmin,pmax =", (pmin,pmax)
            classify = AND(CAT([AA(ISGE)(TRANS([pmin,point])),AA(ISLE)(TRANS([pmax,point]))]))
            if classify: counter.append(k)
        if len(counter) >= dim+1:
            [cells[k].append(i) for k in counter]
            V += [point]
            i += 1
        counter = []
    return V,cells

#-------------------------------------------------------------------
# 2D image example

if __name__=="__main__":
    blocks = [ [[0,0],[5,10]], [[5,0],[9,3]], [[9,0],[13,3]], [[5,3],[8,10]],  [[8,3],[13,10]], [[0,10],[9,12]], [[9,10],[13,12]], [[0,0],[0,12]], [[0,0],[13,0]], [[13,0],[13,12]], [[0,12],[13,12]] ]
    
    model = larFromImageBlocks(blocks)
    V,cells = model
    print "\nV =\n",V
    print "\ncells =\n",cells
    V,faces = larSkeletons(model,dim=2)
    F0V, F1V, F2V = faces
    print "AA(LEN)([F0V, F1V, F2V]) =", AA(LEN)([F0V, F1V, F2V])
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V[:-4])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V[:-4])) ))

#-------------------------------------------------------------------
# 3D image example

if __name__=="__main__":
    blocks = [ [[0,0,0],[5,10,3]], [[5,0,0],[9,3,3]], [[9,0,0],[13,3,3]], [[5,3,0],[8,10,3]],  [[8,3,0],[13,10,3]], [[0,10,0],[9,12,3]], [[9,10,0],[13,12,3]], [[0,0,0],[0,12,3]], [[0,0,0],[13,0,3]], [[13,0,0],[13,12,3]], [[0,12,0],[13,12,3]], [[0,0,0],[13,12,0]], [[0,0,3],[13,12,3]] ]
    
    model = larFromImageBlocks(blocks)
    V,cells = model
    print "\nV =\n",V
    print "\ncells =\n",cells
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    print "AA(LEN)([F0V, F1V, F2V, F3V]) =", AA(LEN)([F0V, F1V, F2V, F3V])
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V[:-6])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))

