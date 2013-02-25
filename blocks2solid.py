from lar import *

solid = [[10, 13, 21, 8], [31, 15, 7, 6], [28, 21, 3, 24], [31, 21, 8, 18], [7, 7, 6, 4], [13, 8, 3, 3], [5, 11, 8, 2], [13, 11, 12, 2], [11, 6, 1, 1], [31, 39, 7, 2], [31, 41, 2, 4], [33, 41, 3, 2], [12, 21, 16, 1], [6, 13, 1, 5], [7, 13, 3, 6], [21, 33, 7, 13], [19, 34, 2, 2], [17, 36, 4, 7], [31, 13, 3, 2], [39, 20, 2, 10], [33, 43, 2, 1], [39, 34, 1, 1], [22, 22, 6, 1], [18, 23, 10, 1], [21, 22, 1, 1], [4, 12, 1, 1], [18, 22, 3, 1], [15, 22, 3, 1], [9, 19, 1, 1], [5, 9, 2, 2], [26, 12, 4, 1], [24, 24, 4, 1], [25, 12, 1, 1], [36, 41, 1, 1], [26, 32, 2, 1], [5, 14, 1, 2], [35, 14, 1, 1], [13, 7, 1, 1], [38, 18, 2, 2], [27, 29, 1, 3], [17, 10, 3, 1], [5, 13, 1, 1], [28, 45, 1, 1], [6, 8, 1, 1], [19, 43, 2, 2], [20, 33, 1, 1], [18, 35, 1, 1], [18, 43, 1, 1], [16, 39, 1, 1], [38, 17, 1, 1], [38, 20, 1, 1], [10, 6, 1, 1], [34, 14, 1, 1], [16, 10, 1, 1], [16, 9, 1, 1], [39, 33, 1, 1], [39, 30, 1, 3]]

empty = [[0, 22, 15, 12], [15, 24, 4, 10], [0, 34, 15, 16], [15, 34, 1, 16], [41, 0, 9, 33], [40, 33, 10, 17], [40, 6, 1, 14], [12, 0, 28, 7], [20, 7, 20, 4], [0, 0, 12, 6], [17, 7, 3, 3], [0, 13, 5, 9], [19, 24, 5, 1], [19, 25, 3, 8], [22, 25, 5, 7], [30, 11, 10, 2], [0, 6, 5, 6], [0, 12, 4, 1], [35, 13, 5, 1], [11, 21, 1, 1], [16, 46, 24, 4], [30, 45, 10, 1], [40, 5, 1, 1], [5, 20, 5, 2], [7, 19, 2, 1], [5, 18, 2, 2], [16, 8, 1, 1], [25, 32, 1, 1], [40, 0, 1, 5], [27, 11, 3, 1], [9, 6, 1, 1], [36, 14, 4, 1], [16, 34, 2, 5], [14, 7, 3, 1], [35, 43, 5, 2], [24, 32, 1, 1], [33, 44, 2, 1], [16, 44, 3, 2], [7, 6, 2, 1], [22, 32, 2, 1], [36, 42, 4, 1], [5, 6, 2, 2], [5, 8, 1, 1], [5, 16, 1, 2], [15, 23, 3, 1], [38, 15, 2, 2], [18, 34, 1, 1], [25, 11, 2, 1], [38, 39, 2, 3], [37, 41, 1, 1], [20, 45, 1, 1], [16, 42, 1, 2], [19, 45, 1, 1], [19, 33, 1, 1], [16, 41, 1, 1], [10, 21, 1, 1], [17, 43, 1, 1], [34, 13, 1, 1], [16, 40, 1, 1], [29, 45, 1, 1], [27, 25, 1, 4], [40, 31, 1, 2], [40, 30, 1, 1], [39, 17, 1, 1], [39, 37, 1, 2], [39, 36, 1, 1], [39, 35, 1, 1]]

"""
solid = [[0,0,4,2],[0,2,2,2],[2,2,2,2]]
empty = []
"""

blocks = solid + empty

print blocks

def writeBlock(cooStore):
    def writeBlock0(block):
        x,y,dx,dy = block
        i = x
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        i = x+dx
        for j in range(y,y+dy+1): cooStore.append([i,j,1])
        j = y
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        j = y+dy
        for i in range(x,x+dx+1): cooStore.append([i,j,1])
        return cooStore
    return writeBlock0

cooStore = []
for block in blocks:
    writeBlock(cooStore)(block)

print cooStore

lilStore = csrCreateFromCoo(cooStore).tolil()

xmax,ymax = lilStore.shape
lilStore[0,0] += 1
lilStore[xmax-1,0] += 1
lilStore[0,ymax-1] += 1
lilStore[xmax-1,ymax-1] += 1


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


print "\n\n", updatedBlock

vertIndex = dict()

index = 0
for block in updatedBlock:
    for vert in block:
        if vert not in vertIndex:
            vertIndex[vert] = index
            index += 1

print "\n", vertIndex


print "\n FV = ", [[vertIndex[vert] for vert in block] for block in updatedBlock]
print "\n V = ", AA(list)(vertIndex.keys())




