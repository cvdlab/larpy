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

insideBlocks = []
outsideBlocks = []

# Input of image

image = mahotas.imread('test1.bmp')
imWidth,imHeight = image.shape
Inside = copy(image)
Outside = copy(image)
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

# field value of p

outsidePoint = image[p]
print "\noutsidePoint =",outsidePoint

# Computation of the sub-image of a point

x0,y0 = p
cross = x0,x0,x0,y0,y0,y0
while cross != crossGrow(cross):
    cross = crossGrow(cross)
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

x,y = subregions[0][1:]
dx,dy = AA(abs)(VECTDIFF([(x,y),(x0,y0)]))
block00 = [x,y] + [dx,dy]
x,y = subregions[1][1:]
dx,dy = AA(abs)(VECTDIFF([(x,y),(x0,y0)]))
block10 = [x0,y] + [dx,dy]
x,y = subregions[2][1:]
dx,dy = AA(abs)(VECTDIFF([(x,y),(x0,y0)]))
block01 = [x,y0] + [dx,dy]
x,y = subregions[3][1:]
dx,dy = AA(abs)(VECTDIFF([(x,y),(x0,y0)]))
block11 = [x0,y0] + [dx,dy]

print "\nblock00,block01,block10,block11 =", (block00,block10,block01,block11)

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
print "\nblocks =", blocks

# write blocks to a file in SVG format

if outsidePoint: outsideBlocks += blocks
else: insideBlocks += blocks
blockImage = outsideBlocks + insideBlocks

def block(data):
	x,y,dx,dy = data
	return T([1,2])([x,y])(CUBOID([dx,dy]))

solid = AA(COLOR(BLACK))(AA(SKELETON(1))(AA(block)(insideBlocks)))
empty = AA(COLOR(WHITE))(AA(SKELETON(1))(AA(block)(outsideBlocks)))

VIEW(STRUCT(solid + empty))