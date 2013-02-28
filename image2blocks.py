from lar import *
from scipy import ndimage
from random import random
from copy import copy
import numpy as np
import pylab
import pymorph
import mahotas
from time import time

#-------------------------------------------------------------------
# data input of 2D image

solids = [[283,388,12,7],[288,395,10,16],[291,411,26,8],[309,401,12,10],[317,411,4,3],[298,405,3,6],[301,409,8,2],[296,419,10,3],[306,419,5,2],[286,395,2,6],[295,391,2,4],[298,399,2,6],[285,386,7,2],[284,395,2,2],[289,411,2,4],[317,414,3,2],[305,408,4,1],[321,402,1,8],[287,401,1,5],[298,396,1,3],[297,393,1,2],[282,391,1,2],[292,387,2,1],[294,419,2,2],[311,419,4,1],[317,416,2,1],[284,387,1,1],[288,385,1,1],[285,397,1,1],[308,407,1,1],[293,419,1,1],[317,417,1,1],[319,399,1,1],[315,397,1,1],[311,399,1,1],[310,400,2,1],[319,400,2,1],[312,398,7,3] ]
voids = [[299,381,27,16],[276,422,50,9],[276,398,10,24],[301,397,8,10],[322,397,4,25],[276,381,6,17],[282,381,17,4],[295,385,4,5],[286,415,4,7],[317,418,5,4],[311,420,6,2],[301,407,4,2],[286,411,3,4],[286,406,2,5],[290,419,3,3],[282,395,2,3],[282,385,2,3],[292,385,3,2],[297,390,2,3],[309,397,3,2],[319,397,3,2],[319,416,3,2],[320,414,2,2],[305,407,3,1],[299,397,2,2],[300,399,1,6],[298,393,1,3],[284,385,4,1],[306,421,5,1],[289,385,3,1],[312,397,3,1],[316,397,3,1],[282,388,1,3],[309,399,1,2],[321,399,1,3],[321,410,1,4],[293,421,3,1],[286,401,1,5],[282,393,1,2],[290,415,1,4],[284,386,1,1],[294,387,1,1],[295,390,2,1],[284,397,1,1],[293,420,1,1],[318,417,1,1],[315,419,2,1],[320,399,1,1],[310,399,1,1] ]
rects = solids + voids
#-------------------------------------------------------------------
# input testing

blocks = [[[rect[0],rect[1]],[rect[0]+rect[2], rect[1]+rect[3]]] for rect in rects]
boundary = [ [[276,381],[326,381]], [[276,381],[276,431]], [[276,431],[326,431]], [[326,326],[326,431]] ]

Image = STRUCT([COLOR(BLUE)(STRUCT([T([1,2])(rect[:2])(CUBOID(rect[2:])) for rect in solids])),
             COLOR(RED)(STRUCT([T([1,2])(rect[:2])(CUBOID(rect[2:])) for rect in voids]))])
IMAGE = S(1)(-1)(R([1,2])(PI/2)(T([1,2])([-276.0, -381.0])(Image)))

#-------------------------------------------------------------------

"""
   Decomposition of a 2D image into rectangular blocks.

"""
def BMPimage2blocks(filename):
    
    def crossGrow(image,cross):
        xm,x,xM,ym,y,yM = cross
        pointValue = image[x,y]
        if xm > 0 and image[xm-1,y] == pointValue: xm -= 1
        if xM < 49 and image[xM+1,y] == pointValue: xM += 1
        if ym > 0 and image[x,ym-1] == pointValue: ym -= 1
        if yM < 49 and image[x,yM+1] == pointValue: yM += 1
        return xm,x,xM,ym,y,yM

    def xbar(image,span,xspan):
        xm,x0,xM,y0 = span
        xmin,xMAX = xspan
        pointValue = image[x0,y0]
        if xm > xmin and image[xm-1,y0] == pointValue: xm -= 1
        if xM < xMAX and image[xM+1,y0] == pointValue: xM += 1
        return xm,x0,xM,y0

    def ybar(image,span,yspan):
        x0,ym,y0,yM = span
        ymin,yMAX = yspan
        pointValue = image[x0,y0]
        if ym > ymin and image[x0,ym-1] == pointValue: ym -= 1
        if yM < yMAX and image[x0,yM+1] == pointValue: yM += 1
        return x0,ym,y0,yM

    def delta(p,q): # area of rectangle (q,p)
        vect = AA(abs)(VECTDIFF([p,q]))
        return vect[0] * vect[1]


    # Input of image into persistent storage structures

    image = mahotas.imread(filename)
    pylab.imshow(image)
    pylab.show()

    imWidth,imHeight = image.shape
    Inside = copy(image)
    Outside = copy(image)
    insideBlocks = []
    outsideBlocks = []

    # Initialization of the sparse matrix used as query point storage

    Points = lil_matrix(csrCreate([range(imWidth) for k in range(imHeight)]), dtype=int)

    def pointSelection(Points):
        indices = np.nonzero(Points)
        nonEmptyRows = list(set(indices[0]))
        numberOfNonEmptyRows = len(nonEmptyRows)
        rowOrdinal = int(numberOfNonEmptyRows * random())
        rowIndex = nonEmptyRows[rowOrdinal]
        numberOfNonZeroElementsInRow = Points[rowIndex,:].nnz
        columnOrdinal = int(numberOfNonZeroElementsInRow * random())
        indices = np.nonzero(Points[rowIndex,:])[1]
        columnIndex = indices[columnOrdinal]
        return rowIndex,columnIndex

    # Generation of a random point

    while True:
    #for k in range(1):
        
        mask0 = copy(image)
        mask1 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
        mask2 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
        mask3 = scipy.zeros((imWidth,imHeight),dtype=image.dtype)
        mask = scipy.zeros((imWidth,imHeight),dtype=int)
        
        #p = int(imWidth*random()),int(imWidth*random())
        if Points.nnz != 0:
            p = pointSelection(Points)
            #p = 39,48
        else: break
        #print "\np =",p
        if Points[p] != 0:
            # reverse p value, in order to show
            image[p] = not(image[p])
            #pylab.imshow(image)
            #pylab.show()
            # reverse back p value to its original value
            image[p] = not(image[p])

            # field value of p

            outsidePoint = image[p]
            if outsidePoint: theImage = Outside
            else: theImage = Inside
            #print "\noutsidePoint =",outsidePoint

            # Computation of the sub-image of a point

            x0,y0 = p
            cross = x0,x0,x0,y0,y0,y0
            while cross != crossGrow(theImage,cross,):
                cross = crossGrow(theImage,cross,)
            xm,x0,xM,ym,y0,yM = cross
            xspan = xm,xM
            yspan = ym,yM

            # Draw the computed subimage of the point (by value inversion)

            if outsidePoint:
                for x in range(xm,xM+1):
                    for y in range(ym,yM+1):
                        mask0[x,y] = not(Outside[x,y])
                for x in range(xm,xM+1): mask0[x,y0] = Outside[x,y0]
                for y in range(ym,yM+1): mask0[x0,y] = Outside[x0,y]
            else:
                for x in range(xm,xM+1):
                    for y in range(ym,yM+1):
                        mask0[x,y] = not(Inside[x,y])
                for x in range(xm,xM+1): mask0[x,y0] = Inside[x,y0]
                for y in range(ym,yM+1): mask0[x0,y] = Inside[x0,y]

            #pylab.imshow(mask0)
            #pylab.show()

            # Computation of the 1D bar-charts of the point

            pointValue = image[p]
            x0,y0 = p

            Y0,Y1 = yspan
            for y in range(Y0,Y1+1):
                span = x0,x0,x0,y
                while span != xbar(theImage,span,xspan):
                    span = xbar(theImage,span,xspan)
                xm,x,xM,y = span
                # drawing of a span
                for x in range(xm,xM+1): mask1[x,y] = not(mask1[x,y])
            #pylab.imshow(mask1)
            #pylab.show()

            X0,X1 = xspan
            for x in range(X0,X1+1):
                span = x,y0,y0,y0
                while span != ybar(theImage,span,yspan):
                    span = ybar(theImage,span,yspan)
                x,ym,y,yM = span
                # drawing of a span
                for y in range(ym,yM+1): mask2[x,y] = not(mask2[x,y])
            #pylab.imshow(mask2)
            #pylab.show()
        
            # Check against all (X0,Y0, X1,Y1) not generated by a pixel of the original image

            if outsidePoint and X1-X0+1 == imHeight and Y1-Y0+1 == imWidth:
                break

            # Visible set: intersection of 1D bar-charts of the point

            for x in range(X0,X1+1):
                for y in range(Y0,Y1+1):
                    if Inside[p]:
                        mask3[x,y] = not(not(mask1[x,y]) and not(mask2[x,y]))
                    elif Outside[p]:
                        mask3[x,y] = not(mask1[x,y] and mask2[x,y])
            #pylab.imshow(mask3)
            #pylab.show()

            # Computation of weights of visible set

            x0,y0 = p
            spans = []
            for x in range(X0,X1+1):
                span = x,y0,y0,y0
                while span != ybar(theImage,span,yspan):
                    span = ybar(theImage,span,yspan)
                spans += [span]

            for span in spans:
                x,ym,y,yM = span
                for y in range(ym,yM+1):
                    if not(mask3[x,y]):
                        mask[x,y] = delta(p,(x,y))
            #print "\nmask =\n", mask[X0:(X1+1),Y0:(Y1+1)]


            # Sorting pixels of visible subregions

            def sortVisible(p,xspan,yspan):
                subregions = []
                (x0,y0),(X0,X1),(Y0,Y1) = p,xspan,yspan
                subregions = []
                subregion = []
                for x in range(X0,x0+1):
                    for y in range(Y0,y0+1):
                        subregion += [[mask[x,y],x,y]]
                subregions += [sorted(subregion,reverse=True)[0]]
                subregion = []
                for x in range(x0,X1+1):
                    for y in range(Y0,y0+1):
                        subregion += [[mask[x,y],x,y]]
                subregions += [sorted(subregion,reverse=True)[0]]
                subregion = []
                for x in range(X0,x0+1):
                    for y in range(y0,Y1+1):
                        subregion += [[mask[x,y],x,y]]
                subregions += [sorted(subregion,reverse=True)[0]]
                subregion = []
                for x in range(x0,X1+1):
                    for y in range(y0,Y1+1):
                        subregion += [[mask[x,y],x,y]]
                subregions += [sorted(subregion,reverse=True)[0]]
                return subregions

            subregions = sortVisible(p,xspan,yspan)

            # block generation (with corrections to eliminate the pixel area)
            # [check border conditions:  possible BUG]

            x00,y00 = subregions[0][1:]
            x10,y10 = subregions[1][1:]
            x10 = x10+1
            x01,y01 = subregions[2][1:]
            y01 = y01+1
            x11,y11 = subregions[3][1:]
            x11,y11 = x11+1,y11+1
            
            # modify position for pixel origin, to maximize cross covering
            def pixelOrigin(x0,y0, x00,y00, x10,y10, x01,y01, x11,y11):
                def area (x0,y0):
                    return abs((x00-x0)*(y00-y0) + (x10-x0)*(y10-y0) + (x01-x0)*(y01-y0) + (x11-x0)*(y11-y0))
                areas  = [ (area(x0,y0), (x0,y0)) ]
                areas += [ (area(x0+1,y0), (x0+1,y0)) ]
                areas += [ (area(x0,y0+1), (x0,y0+1)) ]
                areas += [ (area(x0+1,y0+1), (x0+1,y0+1)) ]
                best = sorted(areas)[-1]
                x0,y0 = best[1]
                return x0,y0
            
            x0,y0 = pixelOrigin(x0,y0, x00,y00, x10,y10, x01,y01, x11,y11)
            
            dx00,dy00 = AA(abs)(VECTDIFF([(x00,y00),(x0,y0)]))
            dx10,dy10 = AA(abs)(VECTDIFF([(x10,y10),(x0,y0)]))
            dx01,dy01 = AA(abs)(VECTDIFF([(x01,y01),(x0,y0)]))
            dx11,dy11 = AA(abs)(VECTDIFF([(x11,y11),(x0,y0)]))

            block00 = [x00,y00] + [dx00,dy00]
            block10 = [x0,y10] + [dx10,dy10]
            block01 = [x01,y0] + [dx01,dy01]
            block11 = [x0,y0] + [dx11,dy11]

            #print "\nblock00,block01,block10,block11 =", (block00,block10,block01,block11)
            
            # 

            # possible joins of the four blocks (Block = [x,y,dx,dy])

            """
                if 0010 and 0111 and 0001 and 1011: A
                elif 0010 and 0111: B
                elif 0010: D
                elif 0111: E
                elif 0001 and 1011: C
                elif 0001: F
                elif 1011: G
                else H
            """
            pred_0010 = block00[3] == block10[3]
            pred_0111 = block01[3] == block11[3]
            pred_0001 = block00[2] == block01[2]
            pred_1011 = block10[2] == block11[2]

            block_A = copy(block00)
            block_A[2] = block00[2] + block10[2]
            block_A[3] = block00[3] + block01[3]
            block_a = copy(block00)
            block_a[2] = block00[2] + block10[2]
            block_b = copy(block01)
            block_b[2] = block01[2] + block11[2]
            block_c = copy(block00)
            block_c[3] = block00[3] + block01[3]
            block_d = copy(block10)
            block_d[3] = block10[3] + block11[3]

            blocks = []
            if pred_0010 and pred_0111 and pred_0001 and pred_1011:
                blocks = [block_A]
            else:
                if pred_0010 and pred_0111:
                    blocks = [block_a,block_b]
                elif pred_0001 and pred_1011:
                    blocks = [block_c,block_d]

                else:
                    if blocks == [] and pred_0010:
                        blocks = [block_a,block01,block11]
                    elif blocks == [] and pred_0111:
                        [block_b,block00,block10]

                    elif blocks == [] and pred_0001:
                        blocks = [block_c,block10,block11]
                    elif blocks == [] and pred_1011:
                        blocks = [block_d,block00,block01]

            if blocks == []: blocks = [block00,block10,block01,block11]
            #print "\nblocks =", blocks

            # filter out empty boxes
        
            blocks = [ block for block in blocks if not( block[2]==0 or block[3]==0) ]
        
            # store current blocks in persistent data structures

            def fill(image,block):
                x0,y0,dx,dy = block
                for x in range(x0,x0+dx):
                    for y in range(y0,y0+dy):
                        image[x,y] = not(pointValue)
                        Points[x,y] = 0
                return image

            if outsidePoint:
                outsideBlocks += blocks
                for block in blocks:
                    Outside = fill(Outside,block)
            else:
                insideBlocks += blocks
                for block in blocks:
                    Inside = fill(Inside,block)

    return insideBlocks,outsideBlocks
#-------------------------------------------------------------------

if __name__ == "__main__":

    # print computing time

    start = time()
    insideBlocks,outsideBlocks = BMPimage2blocks('test1.bmp')
    end = time()
    print "\ntime =", end - start
    print "len(insideBlocks+outsideBlocks) =",len(insideBlocks+outsideBlocks), "\n"

    # display current blocks in pyPlasm

    def block(data):
        x,y,dx,dy = AA(float)(data)
        rectangle = CUBOID([dx,dy])
        return T([1,2])([x,y])(rectangle)

    if insideBlocks != []:
        solids = AA(COLOR(BLACK))(AA(#ID
                                     SKELETON(1)
                                    )(AA(block)(insideBlocks)))
    else: solids = []
    if outsideBlocks != []:
        voids = AA(COLOR(WHITE))(AA(#ID
                                    SKELETON(1)
                                    )(AA(block)(outsideBlocks)))
    else: voids = []

    View(#STRUCT([T(3)(-0.1)(IMAGE),
         STRUCT(voids + solids))#]))
    print "\nsolid =",insideBlocks
    print "\nempty =",outsideBlocks
