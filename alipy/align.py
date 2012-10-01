import star
import pysex
import quad
import os
import numpy as np
import scipy.ndimage
import pyfits



def affineremap(filepath, transform, shape, alifilepath=None, outdir = "alipy_out", makepng=False, verbose=True):
	"""
	Apply the simple affine transform to the image and saves the result as FITS, without using pyraf.
	
	:param filepath: FITS file to align
	:type filepath: string

	:param transform: as returned e.g. by alipy.ident()
	:type transform: SimpleTransform object
	
	:param shape: Output shape (width, height) 
	:type shape: tuple

	:param alifilepath: where to save the aligned image. If None, I put it in the default directory.
	:type alifilepath: string
	
	:param makepng: If True I make a png of the aligned image as well.
	:type makepng: boolean

	"""
	inv = transform.inverse()
	(matrix, offset) = inv.matrixform()
	#print matrix, offset
	
	data, hdr = fromfits(filepath, hdu = 0, verbose = verbose)
	data = scipy.ndimage.interpolation.affine_transform(data, matrix, offset=offset, output_shape = shape)
	
	basename = os.path.splitext(os.path.basename(filepath))[0]
	if not alifilepath:
		#(imgdir, filename) = os.path.split(self.file)
		#(common, ext) = os.path.splitext(filename)
		#filepath = os.path.join(imgdir, common + "_affineremap.fits")
		
		if not os.path.isdir(outdir):
			os.makedirs(outdir)
		alifilepath = os.path.join(outdir, basename + "_affineremap.fits")
	
	tofits(alifilepath, data, hdr = None, verbose = verbose)
	
	if makepng:
		try:
			import f2n
		except ImportError:
			print "Couldn't import f2n -- install it !"
			return
		myimage = f2n.f2nimage(numpyarray=data, verbose=False)
		myimage.setzscale("auto", "auto")
		myimage.makepilimage("log", negative = False)
		myimage.writetitle(basename + "_affineremap.fits")
		if not os.path.isdir(outdir):
				os.makedirs(outdir)
		myimage.tonet(os.path.join(outdir, basename + "_affineremap.png"))



def fromfits(infilename, hdu = 0, verbose = True):
	"""
	Reads a FITS file and returns a 2D numpy array of the data.
	Use hdu to specify which HDU you want (default = primary = 0)
	"""
	
	if verbose:
		print "Reading %s ..." % (os.path.basename(infilename))
	
	pixelarray, hdr = pyfits.getdata(infilename, hdu, header=True)
	pixelarray = np.asarray(pixelarray).transpose()
	
	pixelarrayshape = pixelarray.shape
	if verbose :
		print "FITS import (%i, %i) BITPIX %s / %s" % (pixelarrayshape[0], pixelarrayshape[1], hdr["BITPIX"], str(pixelarray.dtype.name))
		
	return pixelarray, hdr

def tofits(outfilename, pixelarray, hdr = None, verbose = True):
	"""
	Takes a 2D numpy array and write it into a FITS file.
	If you specify a header (pyfits format, as returned by fromfits()) it will be used for the image.
	You can give me boolean numpy arrays, I will convert them into 8 bit integers.
	"""
	pixelarrayshape = pixelarray.shape
	if verbose :
		print "FITS export (%i, %i) %s ..." % (pixelarrayshape[0], pixelarrayshape[1], str(pixelarray.dtype.name))

	if pixelarray.dtype.name == "bool":
		pixelarray = np.cast["uint8"](pixelarray)

	if os.path.isfile(outfilename):
		os.remove(outfilename)
	
	if hdr == None: # then a minimal header will be created 
		hdu = pyfits.PrimaryHDU(pixelarray.transpose())
	else: # this if else is probably not needed but anyway ...
		hdu = pyfits.PrimaryHDU(pixelarray.transpose(), hdr)

	hdu.writeto(outfilename)
	
	if verbose :
		print "Wrote %s" % outfilename

		
