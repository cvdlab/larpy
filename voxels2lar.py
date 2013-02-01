from lar import *

def larImportImageBlocks(blocks):
    cells = [[] for k in range(len(blocks))]
    print "\ncells =\n",cells
    bounds = mat(blocks).T
    boundingBox = (bounds[0].min(), bounds[1].min(), bounds[2].max(), bounds[3].max())
    print "\nboundingBox =\n",boundingBox
    counter,V,i = [],[],0
    for x in range(boundingBox[0],boundingBox[2]+1):
        for y in range(boundingBox[1],boundingBox[3]+1):
            for k,cell in enumerate(blocks):
                xmin,ymin,xmax,ymax = cell
                if (xmin <= x <= xmax) and (ymin <= y <= ymax): counter.append(k)
            if len(counter) >= 3:
                [cells[k].append(i) for k in counter]
                V += [[x,y]]
                i += 1
            counter = []
    return V,cells


if __name__=="__main__":
    blocks = [ [0,0,5,10], [5,0,9,3], [9,0,13,3], [5,3,8,10],  [8,3,13,10], [0,10,9,12], [9,10,13,12], [0,0,0,12], [0,0,13,0], [13,0,13,12], [0,12,13,12] ]

    model = larImportImageBlocks(blocks)
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


def larImportImageBlocks(blocks):
    dim = len(blocks[0])/2
    cells = [[] for k in range(len(blocks))]
    print "\ncells =\n",cells
    bounds = TRANS(blocks)
    #boundingBox = (min(bounds[0]), bounds[1].min(), bounds[2].max(), bounds[3].max())
    boundingBox = CAT([AA(min)(bounds[:dim]), AA(max)(bounds[dim:])])
    print "\nboundingBox =\n",boundingBox
    counter,V,i = [],[],0
    for x in range(boundingBox[0],boundingBox[3]+1):
        for y in range(boundingBox[1],boundingBox[4]+1):
            for z in range(boundingBox[2],boundingBox[5]+1):
                for k,cell in enumerate(blocks):
                    xmin,ymin,zmin,xmax,ymax,zmax = cell
                    if (xmin <= x <= xmax) and (ymin <= y <= ymax) and (zmin <= z <= zmax): counter.append(k)
                if len(counter) >= 4:
                    [cells[k].append(i) for k in counter]
                    V += [[x,y,z]]
                    i += 1
                counter = []
    return V,cells


if __name__=="__main__":
    blocks = [ [0,0,0,5,10,3], [5,0,0,9,3,3], [9,0,0,13,3,3], [5,3,0,8,10,3],  [8,3,0,13,10,3], [0,10,0,9,12,3], [9,10,0,13,12,3], [0,0,0,0,12,3], [0,0,0,13,0,3], [13,0,0,13,12,3], [0,12,0,13,12,3], [0,0,0,13,12,0], [0,0,3,13,12,3] ]
    
    model = larImportImageBlocks(blocks)
    print "\nmodel =\n",model
    V,faces = larSkeletons(model,dim=3)
    F0V, F1V, F2V, F3V = faces
    print "AA(LEN)([F0V, F1V, F2V, F3V]) =", AA(LEN)([F0V, F1V, F2V, F3V])
    V = model[0]
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F3V[:-6])) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F2V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F1V)) ))
    VIEW(EXPLODE(1.2,1.2,1.2)( MKPOLS((V,F0V+F1V+F2V+F3V[:-6])) ))
