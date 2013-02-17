import numpy as np
import scipy
import pylab
import pymorph
import mahotas
from scipy import ndimage
from random import random
from copy import copy

"""
   Decomposition of a 2D image into rectangular blocks.

"""
def crossGrow(cross):
    xm,x,xM,ym,y,yM = cross
    pointValue = test[x,y]
    if xm > 0 and test[xm-1,y] == pointValue: xm -= 1
    if xM < 49 and test[xM+1,y] == pointValue: xM += 1
    if ym > 0 and test[x,ym-1] == pointValue: ym -= 1
    if yM < 49 and test[x,yM+1] == pointValue: yM += 1
    return xm,x,xM,ym,y,yM

def xbar(span,xspan):
    xm,x0,xM,y0 = span
    xmin,xMAX = xspan
    pointValue = test[x0,y0]
    if xm > xmin and test[xm-1,y0] == pointValue: xm -= 1
    if xM < xMAX and test[xM+1,y0] == pointValue: xM += 1
    return xm,x0,xM,y0

def ybar(span,yspan):
    x0,ym,y0,yM = span
    ymin,yMAX = yspan
    pointValue = test[x0,y0]
    if ym > ymin and test[x0,ym-1] == pointValue: ym -= 1
    if yM < yMAX and test[x0,yM+1] == pointValue: yM += 1
    return x0,ym,y0,yM


# Input of image

test = mahotas.imread('test1.bmp')
imWidth,imHeight = test.shape
mask0 = copy(test)
mask1 = scipy.ones((imWidth,imHeight),dtype=test.dtype)
mask2 = scipy.ones((imWidth,imHeight),dtype=test.dtype)

# Generation of a random point

p = int(imWidth*random()),int(imWidth*random())
# reverse p value, in order to show
test[p] = not(test[p])
pylab.imshow(test)
pylab.show()
# reverse back p value to its original value
test[p] = not(test[p])

# Computation of the sub-image of a point

x0,y0 = p
cross = x0,x0,x0,y0,y0,y0
while cross != crossGrow(cross):
    cross = crossGrow(cross)
xm,x0,xM,ym,y0,yM = cross
xspan = xm,xM
yspan = ym,yM
print "\n x0,y0, xspan, yspan =", (x0,y0, xspan, yspan)

# Draw the computed subimage of the point (by value inversion)

for x in range(xm,xM+1):
    for y in range(ym,yM+1):
        mask0[x,y] = not(test[x,y])
for x in range(xm,xM+1): mask0[x,y0] = test[x,y0]
for y in range(ym,yM+1): mask0[x0,y] = test[x0,y]

pylab.imshow(mask0)
pylab.show()

# Computation of the 1D bar-charts of the point

pointValue = test[p]
x0,y0 = p
Y0,Y1 = yspan
for y in range(Y0,Y1+1):
    span = x0,x0,x0,y
    while span != xbar(span,xspan):
        span = xbar(span,xspan)
    xm,x,xM,y = span
    print " x,y, xspan =", (x,y, (xm,xM))
    # drawing of a span
    for x in range(xm,xM+1): mask1[x,y] = not(mask1[x,y])
pylab.imshow(mask1)
pylab.show()

pointValue = test[p]
x0,y0 = p
X0,X1 = xspan
for x in range(X0,X1+1):
    span = x,y0,y0,y0
    while span != ybar(span,yspan):
        span = ybar(span,yspan)
    x,ym,y,yM = span
    print " x,y, yspan =", (x,y, (ym,yM))
    # drawing of a span
    for y in range(ym,yM+1): mask2[x,y] = not(mask2[x,y])

print "\n x0,y0, xspan, yspan =", (x0,y0, xspan, yspan)

pylab.imshow(mask2)
pylab.show()

