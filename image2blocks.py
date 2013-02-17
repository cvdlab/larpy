import numpy as np
import scipy
import pylab
import pymorph
import mahotas
from scipy import ndimage
from random import random
from copy import copy
from pyplasm import *

"""
   Decomposition of a 2D image into rectangular blocks.

"""
def crossGrow(cross):
    xm,x,xM,ym,y,yM = cross
    pointValue = image[x,y]
    if xm > 0 and image[xm-1,y] == pointValue: xm -= 1
    if xM < 49 and image[xM+1,y] == pointValue: xM += 1
    if ym > 0 and image[x,ym-1] == pointValue: ym -= 1
    if yM < 49 and image[x,yM+1] == pointValue: yM += 1
    return xm,x,xM,ym,y,yM

def xbar(span,xspan):
    xm,x0,xM,y0 = span
    xmin,xMAX = xspan
    pointValue = image[x0,y0]
    if xm > xmin and image[xm-1,y0] == pointValue: xm -= 1
    if xM < xMAX and image[xM+1,y0] == pointValue: xM += 1
    return xm,x0,xM,y0

def ybar(span,yspan):
    x0,ym,y0,yM = span
    ymin,yMAX = yspan
    pointValue = image[x0,y0]
    if ym > ymin and image[x0,ym-1] == pointValue: ym -= 1
    if yM < yMAX and image[x0,yM+1] == pointValue: yM += 1
    return x0,ym,y0,yM

def delta(p,q): # area of rectangle (q,p)
    vect = AA(abs)(VECTDIFF([p,q]))
    return vect[0] * vect[1]


# Input of image

image = mahotas.imread('test1.bmp')
imWidth,imHeight = image.shape
mask0 = copy(image)
mask1 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
mask2 = scipy.ones((imWidth,imHeight),dtype=image.dtype)
mask3 = scipy.zeros((imWidth,imHeight),dtype=image.dtype)
mask = scipy.zeros((imWidth,imHeight),dtype=int)

# Generation of a random point

p = int(imWidth*random()),int(imWidth*random())
# reverse p value, in order to show
image[p] = not(image[p])
pylab.imshow(image)
pylab.show()
# reverse back p value to its original value
image[p] = not(image[p])

# Computation of the sub-image of a point

x0,y0 = p
cross = x0,x0,x0,y0,y0,y0
while cross != crossGrow(cross):
    cross = crossGrow(cross)
xm,x0,xM,ym,y0,yM = cross
xspan = xm,xM
yspan = ym,yM

# Draw the computed subimage of the point (by value inversion)

for x in range(xm,xM+1):
    for y in range(ym,yM+1):
        mask0[x,y] = not(image[x,y])
for x in range(xm,xM+1): mask0[x,y0] = image[x,y0]
for y in range(ym,yM+1): mask0[x0,y] = image[x0,y]

pylab.imshow(mask0)
pylab.show()

# Computation of the 1D bar-charts of the point

pointValue = image[p]
x0,y0 = p
Y0,Y1 = yspan
for y in range(Y0,Y1+1):
    span = x0,x0,x0,y
    while span != xbar(span,xspan):
        span = xbar(span,xspan)
    xm,x,xM,y = span
    # drawing of a span
    for x in range(xm,xM+1): mask1[x,y] = not(mask1[x,y])
pylab.imshow(mask1)
pylab.show()

pointValue = image[p]
x0,y0 = p
X0,X1 = xspan
for x in range(X0,X1+1):
    span = x,y0,y0,y0
    while span != ybar(span,yspan):
        span = ybar(span,yspan)
    x,ym,y,yM = span
    # drawing of a span
    for y in range(ym,yM+1): mask2[x,y] = not(mask2[x,y])
pylab.imshow(mask2)
pylab.show()

# Visible set: intersection of 1D bar-charts of the point

for x in range(X0,X1+1):
    for y in range(Y0,Y1+1):
        mask3[x,y] = not(mask1[x,y]) and not(mask2[x,y])  # BUG (check ! )
pylab.imshow(mask3)
pylab.show()

# Computation of weights of visible set

x0,y0 = p
spans = []
for x in range(X0,X1+1):
    span = x,y0,y0,y0
    while span != ybar(span,yspan):
        span = ybar(span,yspan)
    spans += [span]

for span in spans:
    x,ym,y,yM = span
    for y in range(ym,yM+1):
        mask[x,y] = delta(p,(x,y))
print "\nmask =\n", mask[X0:(X1+1),Y0:(Y1+1)]


# Sorting of visible subregions

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

# block generation and possible join



