from lar import *
import collections
from time import time


solid = [[6, 8, 8, 8], [14, 12, 16, 4], [10, 16, 4, 5], [14, 16, 24, 6], [28, 22, 11, 14], [18, 36, 20, 5], [38, 36, 1, 3], [17, 41, 19, 2], [18, 43, 17, 1], [30, 13, 4, 3], [21, 44, 12, 2], [7, 16, 3, 3], [10, 6, 2, 2], [8, 7, 2, 1], [18, 22, 10, 2], [7, 7, 1, 1], [9, 19, 1, 1], [23, 11, 2, 1], [4, 12, 2, 1], [14, 10, 6, 2], [39, 20, 2, 10], [34, 14, 1, 1], [35, 15, 3, 1], [19, 34, 9, 2], [26, 32, 2, 2], [18, 35, 1, 1], [14, 8, 2, 2], [26, 24, 2, 1], [21, 33, 5, 1], [5, 15, 1, 1], [22, 11, 1, 1], [34, 15, 1, 1], [21, 11, 1, 1], [20, 11, 1, 1], [6, 17, 1, 1], [25, 24, 1, 1], [15, 22, 3, 1], [5, 14, 1, 1], [24, 24, 1, 1], [12, 7, 2, 1], [35, 14, 1, 1], [17, 40, 1, 1], [39, 33, 1, 2], [19, 44, 2, 1], [17, 38, 1, 2], [5, 10, 1, 2], [6, 16, 1, 1], [39, 31, 1, 2], [38, 17, 1, 5], [36, 41, 1, 1], [13, 21, 1, 1], [5, 9, 1, 1], [17, 36, 1, 2], [5, 13, 1, 1], [39, 19, 1, 1], [12, 21, 1, 1], [39, 30, 1, 1], [16, 39, 1, 1], [39, 18, 1, 1], [27, 29, 1, 3], [20, 33, 1, 1], [16, 9, 1, 1]]

empty = [[0, 0, 26, 6], [26, 0, 24, 12], [0, 16, 5, 29], [5, 23, 11, 22], [0, 45, 5, 5], [5, 45, 14, 5], [41, 12, 9, 31], [36, 43, 14, 7], [19, 46, 12, 4], [31, 45, 5, 5], [0, 10, 1, 1], [1, 6, 4, 5], [0, 11, 1, 5], [1, 11, 3, 5], [9, 6, 1, 1], [34, 12, 7, 2], [40, 33, 1, 10], [17, 6, 9, 4], [20, 10, 6, 1], [5, 19, 4, 1], [5, 20, 5, 3], [16, 24, 8, 1], [16, 25, 3, 9], [19, 25, 8, 7], [7, 6, 2, 1], [12, 6, 5, 1], [34, 44, 2, 1], [4, 11, 1, 1], [4, 13, 1, 3], [31, 12, 3, 1], [10, 21, 1, 1], [11, 22, 4, 1], [11, 21, 1, 1], [37, 14, 4, 1], [36, 14, 1, 1], [5, 6, 2, 3], [10, 22, 1, 1], [37, 41, 3, 2], [39, 15, 1, 2], [40, 17, 1, 3], [0, 9, 1, 1], [21, 32, 5, 1], [27, 26, 1, 3], [35, 43, 1, 1], [16, 34, 2, 2], [18, 34, 1, 1], [19, 45, 2, 1], [18, 44, 1, 1], [16, 43, 2, 2], [5, 18, 2, 1], [30, 12, 1, 1], [40, 31, 1, 2], [16, 36, 1, 3], [39, 17, 1, 1], [15, 7, 2, 1], [25, 11, 1, 1], [33, 44, 1, 1], [0, 6, 1, 3], [39, 36, 1, 5], [14, 7, 1, 1], [40, 30, 1, 1], [16, 41, 1, 2], [36, 42, 1, 1], [40, 16, 1, 1], [40, 15, 1, 1], [38, 40, 1, 1], [38, 16, 1, 1], [38, 39, 1, 1], [20, 32, 1, 1], [19, 33, 1, 1], [17, 23, 1, 1], [38, 15, 1, 1], [5, 17, 1, 1], [5, 16, 1, 1], [27, 25, 1, 1], [39, 35, 1, 1], [19, 32, 1, 1], [16, 40, 1, 1], [16, 8, 1, 1], [16, 23, 1, 1]]

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



