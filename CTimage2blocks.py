import dicom
import pylab
from lar import *
from copy import copy
from random import random
from scipy import ndimage,sparse,amax,amin
import numpy as np

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

def quantize (min,max,nparts):
    scalingFact = float(nparts)/(max-min)
    print "\n scalingFact =", scalingFact
    def quantize0 (vect):
        vect = AA(int)((vect-min)*scalingFact + 0.5)
        print "\n vect =", vect
        vect = np.array(vect)
        print "\n vect =", vect
        vect = vect/scalingFact + min
        return vect
    return quantize0

# Input of CT image into persistent storage structures

ct = dicom.read_file("CT_small.dcm")
image = copy(ct.pixel_array)
imWidth,imHeight = image.shape
pylab.imshow(image, cmap=pylab.cm.bone)
pylab.show()

Visited = copy(ct.pixel_array)
for i in range(128):
    for j in range(128):
        Visited[i,j] = 0

# Initialization of the sparse matrix used as query point storage

Points = sparse.lil_matrix(csrCreate([range(imWidth) for k in range(imHeight)]), dtype=int)

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

#while True:
for k in range(1):
    
    mask0 = scipy.zeros((imWidth,imHeight),dtype=image.dtype)
    mask1 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
    mask2 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
    mask3 = scipy.zeros((imWidth,imHeight),dtype=image.dtype)
    mask = scipy.zeros((imWidth,imHeight),dtype=int)
    
    #p = int(imWidth*random()),int(imWidth*random())
    if Points.nnz != 0:
        p = pointSelection(Points)
    #p = 39,48
    else: break

    print "\np =",p
    if Points[p] != 0:
        # reverse p value, in order to show
        Visited[p] = image[p]
        pylab.imshow(Visited)
        pylab.show()

        # Computation of the sub-image of a point

        x0,y0 = p
        cross = x0,x0,x0,y0,y0,y0
        while cross != crossGrow(image,cross,):
            cross = crossGrow(image,cross,)
        xm,x0,xM,ym,y0,yM = cross
        xspan = xm,xM
        yspan = ym,yM
            
        # Draw the computed subimage of the point 
        
        for x in range(xm,xM+1):
            for y in range(ym,yM+1):
                mask0[x,y] = image[x,y]
        for x in range(xm,xM+1): mask0[x,y0] = image[x,y0]
        for y in range(ym,yM+1): mask0[x0,y] = image[x0,y]

        pylab.imshow(mask0)
        pylab.show()

