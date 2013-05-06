"""
To import a stack of PNG images into a 3D array, with denoising and color quatization.

Return a scipy ndarray.
"""
import sys
import numpy as np
from numpy import reshape,array
from scipy.cluster.vq import kmeans,vq
import png
import matplotlib.pyplot as plt
from scipy import ndimage

	
def pngstack2array3d(path, MIN_SLICE, MAX_SLICE, colors):
	"""
	To import a stack of PNG images into a 3D array, with denoising and color quatization.
	
	path = the path to the directory of images
	MIN_SLICE, MAX_SLICE = the numbers of first and last slice.
	Return a scipy ndarray.
	"""

	# -----------------------------------------------------------------------------
	# import images in a 3D array -------------------------------------------------
	# -----------------------------------------------------------------------------
	image2d = []
	for slice in range(MIN_SLICE, MAX_SLICE):
		filename = path+str(slice)+'.png'
		print "filename =",filename
		r=png.Reader(filename)
		content = r.read()
		page = [list(row) for k,row in enumerate(content[2])]
		image2d.append(page)
	
	image3d = array(image2d,dtype='uint8')

	if __name__ == "__main__":
		image3d.shape, image3d.dtype
		image3d
		plt.imshow(image3d[0])
		plt.show()

	# -----------------------------------------------------------------------------
	# selecting colors for quantization, via clustering on first image ------------
	# -----------------------------------------------------------------------------
	# reshaping the pixels matrix
	pixel = reshape(image3d[0],(image3d[0].shape[0]*image3d[0].shape[1],1))
	# performing the clustering
	centroids,_ = kmeans(pixel,colors) # "colors" colors will be found

	# -----------------------------------------------------------------------------
	# -----------------------------------------------------------------------------
	# -----------------------------------------------------------------------------

	for page in range(image3d.shape[0]):

		# image denoising 
		# -------------------------------------------------------------------------
		image3d[page] = ndimage.median_filter(image3d[page], 10)

		# field quantization 
		# -------------------------------------------------------------------------
		# reshaping the pixels matrix
		pixel = reshape(image3d[page],(image3d[page].shape[0]*image3d[page].shape[1],1))
		# quantization
		qnt,_ = vq(pixel,centroids)
		# reshaping the result of the quantization
		centers_idx = np.reshape(qnt,image3d[page].shape)
		image3d[page] = centroids[centers_idx].reshape(image3d[page].shape)

		if __name__ == "__main__":
			# page show
			plt.imshow(image3d[page])
			plt.show()
			print image3d[page][1]

	# return a scipy ndarray 
	# -------------------------------------------------------------------------
	return image3d,colors
	
	

if __name__ == "__main__":
	
	# -----------------------------------------------------------------------------
	# import arguments from command line ------------------------------------------
	# -----------------------------------------------------------------------------
	if sys.argv[1] != 0: 
		IMAGE_PATH = sys.argv[1]
	else: IMAGE_PATH = './'
	if sys.argv[2] and sys.argv[3] != 0:
		MIN_SLICE, MAX_SLICE = int(sys.argv[2]),int(sys.argv[3])
	else: MIN_SLICE, MAX_SLICE = '00','01'
	if sys.argv[4] != 0: 
		colors = int(sys.argv[4])
	else: colors = 2
	
	pngstack2array3d(IMAGE_PATH, MIN_SLICE, MAX_SLICE,colors)
	
"""
== Conversione di un archivio jpeg in stack di png

Import all images to photoshop.
File > Scripts > Load Files into stack
Go to Layers panel and select all layers
Got to animation panel and Select "Create frame animation"
Go to Animation panel and select the small menu symbol on the top right of panel
Select "Make frames from layers"
File > Export > Render video
Select folder
Select Photoshop Image Sequence
Select PNG
Render
"""
