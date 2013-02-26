from lar import *
import collections
from time import time


solid = [[28, 14, 8, 18], [36, 18, 4, 14], [27, 32, 9, 11], [36, 32, 3, 7], [18, 35, 9, 9], [10, 12, 18, 9], [28, 12, 2, 2], [7, 7, 3, 12], [5, 9, 2, 7], [6, 16, 1, 2], [10, 10, 10, 2], [15, 21, 13, 2], [18, 23, 10, 1], [27, 43, 8, 2], [10, 8, 6, 2], [16, 9, 1, 1], [9, 19, 1, 1], [36, 15, 2, 3], [13, 21, 2, 1], [20, 33, 7, 2], [24, 24, 4, 1], [21, 44, 6, 2], [28, 45, 1, 1], [36, 39, 2, 2], [33, 13, 1, 1], [16, 39, 2, 1], [17, 38, 1, 1], [24, 11, 1, 1], [31, 13, 2, 1], [19, 44, 2, 1], [38, 17, 1, 1], [36, 41, 1, 1], [12, 7, 2, 1], [26, 32, 1, 1], [23, 11, 1, 1], [10, 6, 2, 2], [12, 21, 1, 1], [6, 8, 1, 1], [40, 24, 1, 6], [39, 32, 1, 3], [22, 11, 1, 1], [17, 40, 1, 3], [19, 34, 1, 1], [27, 31, 1, 1], [20, 11, 2, 1], [30, 13, 1, 1], [40, 21, 1, 3], [17, 36, 1, 2], [27, 30, 1, 1], [4, 12, 1, 1], [27, 29, 1, 1], [27, 45, 1, 1], [40, 20, 1, 1]]

empty = [[12, 0, 38, 7], [20, 7, 16, 4], [36, 7, 14, 7], [38, 14, 1, 3], [39, 14, 11, 4], [41, 18, 9, 32], [0, 21, 11, 28], [11, 23, 5, 26], [10, 49, 1, 1], [11, 49, 30, 1], [0, 0, 4, 19], [0, 19, 7, 2], [2, 49, 8, 1], [29, 45, 6, 1], [35, 44, 6, 2], [16, 46, 19, 3], [35, 46, 6, 3], [14, 7, 6, 1], [4, 0, 8, 6], [0, 49, 2, 1], [38, 41, 1, 1], [39, 39, 2, 3], [37, 42, 2, 2], [39, 42, 2, 2], [39, 35, 2, 4], [5, 6, 5, 1], [4, 7, 2, 2], [6, 7, 1, 1], [16, 23, 1, 13], [17, 25, 9, 8], [17, 33, 3, 1], [11, 22, 4, 1], [26, 27, 1, 5], [27, 25, 1, 4], [29, 11, 7, 1], [23, 24, 1, 1], [9, 20, 1, 1], [25, 11, 4, 1], [37, 41, 1, 1], [7, 19, 2, 2], [38, 39, 1, 2], [17, 8, 3, 2], [4, 16, 2, 3], [34, 12, 2, 2], [16, 40, 1, 6], [26, 26, 1, 1], [17, 35, 1, 1], [32, 12, 2, 1], [18, 24, 5, 1], [17, 44, 2, 2], [36, 43, 1, 1], [30, 12, 2, 1], [18, 34, 1, 1], [16, 36, 1, 3], [37, 14, 1, 1], [19, 45, 2, 1], [36, 14, 1, 1], [26, 25, 1, 1], [17, 24, 1, 1], [6, 18, 1, 1], [16, 8, 1, 1], [17, 34, 1, 1], [4, 14, 1, 2], [36, 42, 1, 1], [35, 43, 1, 1], [4, 11, 1, 1], [40, 19, 1, 1], [11, 21, 1, 1], [40, 31, 1, 4], [4, 9, 1, 2], [40, 30, 1, 1], [4, 13, 1, 1], [4, 6, 1, 1], [40, 18, 1, 1], [17, 23, 1, 1], [17, 43, 1, 1]]

"""
solid = [[0,0,4,2],[0,2,2,2],[2,2,2,2]]
empty = []
    """

nsolid = len(solid)
nempty = len(empty)

blocks = solid + empty

def writeBlock(cooStore):
    def writeBlock0(block):
        cooStore.append([0,0,2])
        x,y,dx,dy = block
        i = x
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        i = x+dx
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        j = y
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        j = y+dy
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        cooStore.append([i,j,2])
        return cooStore
    return writeBlock0

start = time()

cooStore = []
for block in blocks:
    writeBlock(cooStore)(block)

print cooStore

lilStore = csrCreateFromCoo(cooStore).tolil()

def readBlock(lilStore):
    def readBlock0(block):
        x,y,dx,dy = block
        outBlock = [[(i,y) for i in range(x,x+dx) if lilStore[i,y] > 2]]
        outBlock.append([(x+dx,j) for j in range(y,y+dy) if lilStore[x+dx,j] > 2])
        outBlock.append([(i,y+dy) for i in range(x+dx,x,-1) if lilStore[i,y+dy] > 2])
        outBlock.append([(x,j) for j in range(y+dy,y,-1) if lilStore[x,j] > 2])
        return CAT(outBlock)
    return readBlock0


updatedBlock = [readBlock(lilStore)(block) for block in blocks]

verts = collections.OrderedDict(); k = 0

index = 0
for block in updatedBlock:
    for vert in block:
        if vert not in verts:
            verts[vert] = k
            k += 1


V = AA(list)(AA(AA(float))(verts.keys()))
FV = [[verts[vert] for vert in block] for block in updatedBlock]

end = time()

print "\ntotal time (sec) =", end-start
print "\nV = ", V
print "\nFV = ", FV, "\n"

solids = AA(COLOR(BLUE))(MKPOLS((V,FV[:nsolid])))
voids = AA(COLOR(RED))(MKPOLS((V,FV[nsolid:])))

VIEW(EXPLODE(1.2,1.2,1.2)(solids + voids))



