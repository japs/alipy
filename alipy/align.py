import star
import os
import numpy as np
import math
import scipy.ndimage
import pyfits
import csv



def affineremap(filepath, transform, shape, alifilepath=None, outdir = "alipy_out", makepng=False, hdu=0, verbose=True):
	"""
	Apply the simple affine transform to the image and saves the result as FITS, without using pyraf.
	
	:param filepath: FITS file to align
	:type filepath: string

	:param transform: as returned e.g. by alipy.ident()
	:type transform: SimpleTransform object
	
	:param shape: Output shape (width, height) 
	:type shape: tuple

	:param alifilepath: where to save the aligned image. If None, I will put it in the outdir directory.
	:type alifilepath: string
	
	:param makepng: If True I make a png of the aligned image as well.
	:type makepng: boolean

	:param hdu: The hdu of the fits file that you want me to use. 0 is primary. If multihdu, 1 is usually science.


	"""
	inv = transform.inverse()
	(matrix, offset) = inv.matrixform()
	#print matrix, offset
	
	data, hdr = fromfits(filepath, hdu = hdu, verbose = verbose)
	data = scipy.ndimage.interpolation.affine_transform(data, matrix, offset=offset, output_shape = shape)
	
	basename = os.path.splitext(os.path.basename(filepath))[0]
	
	if alifilepath == None:
		alifilepath = os.path.join(outdir, basename + "_affineremap.fits")
	else:	
		outdir = os.path.split(alifilepath)[0]
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		
	
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
		myimage.writetitle(os.path.basename(alifilepath))
		if not os.path.isdir(outdir):
				os.makedirs(outdir)
		myimage.tonet(os.path.join(outdir, os.path.basename(alifilepath)+".png"))



def shape(filepath, hdu = 0, verbose=True):
	"""
	Returns the 2D shape (width, height) of a FITS image.
	
	:param hdu: The hdu of the fits file that you want me to use. 0 is primary. If multihdu, 1 is usually science.

	
	"""
	hdr = pyfits.getheader(filepath, hdu)
	if hdr["NAXIS"] != 2:
		raise RuntimeError("Hmm, this hdu is not a 2D image !")
	if verbose:
		print "Image shape of %s : (%i, %i)" % (os.path.basename(filepath), int(hdr["NAXIS1"]), int(hdr["NAXIS2"]))
	return (int(hdr["NAXIS1"]), int(hdr["NAXIS2"]))
	

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

		

def irafalign(filepath, uknstarlist, refstarlist, shape, alifilepath=None, outdir = "alipy_out", makepng=False, hdu=0, verbose=True):
	"""
	Uses iraf geomap and gregister to align the image. Three steps :
	 * Write the matched source lists into an input file for geomap
	 * Compute a geomap transform from these stars. 
	 * Run gregister
	
	:param filepath: FITS file to be aligned
	:type filepath: string
	
	:param uknstarlist: A list of stars from the "unknown" image to be aligned, that matches to ...
	:type uknstarlist: list of Star objects
	:param refstarlist: ... the list of corresponding stars in the reference image.
	:type refstarlist: list of Star objects
	
	:param shape: Output shape (width, height) 
	:type shape: tuple

	:param alifilepath: where to save the aligned image. If None, I put it in the default directory.
	:type alifilepath: string
	
	:param makepng: If True I make a png of the aligned image as well.
	:type makepng: boolean

	:param hdu: The hdu of the fits file that you want me to use. 0 is primary. If multihdu, 1 is usually science.


	
	"""

	try:
		from pyraf import iraf
	except ImportError:
		print "Couldn't import pyraf !"
		return

	assert len(uknstarlist) == len(refstarlist)
	if len(uknstarlist) < 2:
		if verbose:
			print "Not enough stars for using geomap !"
		return
		
	basename = os.path.splitext(os.path.basename(filepath))[0]
	geomapinpath = basename + ".geomapin"
	geodatabasepath = basename + ".geodatabase"
	if os.path.isfile(geomapinpath):
		os.remove(geomapinpath)
	if os.path.isfile(geodatabasepath):
		os.remove(geodatabasepath)

	
	# Step 1, we write the geomap input.
	table = []	
	for (uknstar, refstar) in zip(uknstarlist, refstarlist):
		table.append([refstar.x, refstar.y, uknstar.x, uknstar.y])
	geomap = open(geomapinpath, "wb") # b needed for csv
	writer = csv.writer(geomap, delimiter="\t")
	writer.writerows(table)
	geomap.close()
	
	
	# Step 2, geomap
		
	iraf.unlearn(iraf.geomap)	
	iraf.geomap.fitgeom = "rscale"		# You can change this to : shift, xyscale, rotate, rscale
	iraf.geomap.function = "polynomial"	# Surface type
	iraf.geomap.maxiter = 3			# Maximum number of rejection iterations
	iraf.geomap.reject = 3.0		# Rejection limit in sigma units
	
	# other options you could specify :
	#(xxorder=                    2) Order of x fit in x
	#(xyorder=                    2) Order of x fit in y
	#(xxterms=                 half) X fit cross terms type
	#(yxorder=                    2) Order of y fit in x
	#(yyorder=                    2) Order of y fit in y
	#(yxterms=                 half) Y fit cross terms type
	#(calctyp=                 real) Computation type

	iraf.geomap.transfo = "broccoli"	# keep it
	iraf.geomap.interac = "no"		# keep it
	iraf.geomap.verbose = "yes"		# keep it
	#iraf.geomap.results = "bla.summary" # The optional results summary files
	
	geomapblabla = iraf.geomap(input=geomapinpath, database=geodatabasepath, xmin = 1, xmax = shape[0], ymin = 1, ymax = shape[1], Stdout=1)
	
	# We read this output ...
	for line in geomapblabla:
		if "X and Y scale:" in line:
			mapscale = line.split()[4:6]
		if "Xin and Yin fit rms:" in line:
			maprmss = line.split()[-2:]
		if "X and Y axis rotation:" in line:
			mapangles = line.split()[-4:-2]
		if "X and Y shift:" in line:
			mapshifts = line.split()[-4:-2]
	
	geomaprms = math.sqrt(float(maprmss[0])*float(maprmss[0]) + float(maprmss[1])*float(maprmss[1]))
	geomapangle = float(mapangles[0])# % 360.0
	geomapscale = 1.0/float(mapscale[0])
	
	if mapscale[0] != mapscale[1]:
		raise RuntimeError("Error reading geomap scale")
	if verbose:
		print "IRAF geomap : Rotation %+11.6f [deg], scale %8.6f, RMS %.3f [pixel]" % (geomapangle, geomapscale, geomaprms)
	
	# Step 3
	
	if alifilepath == None:
		alifilepath = os.path.join(outdir, basename + "_gregister.fits")
	else:	
		outdir = os.path.split(alifilepath)[0]
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	if os.path.isfile(alifilepath):
		os.remove(alifilepath)
	

	iraf.unlearn(iraf.gregister)
	iraf.gregister.geometry = "geometric"	# linear, distortion, geometric
	iraf.gregister.interpo = "spline3"	# linear, spline3
	iraf.gregister.boundary = "constant"	# padding with zero
	iraf.gregister.constant = 0.0
	iraf.gregister.fluxconserve = "yes"
	
	if verbose:
		print "IRAF gregister ..."

	regblabla = iraf.gregister(input = '%s[%s]' % (filepath, hdu), output = alifilepath, database = geodatabasepath, transform = "broccoli", Stdout=1)

	if verbose:
		print "IRAF gregister done !"
	
	if os.path.isfile(geomapinpath):
		os.remove(geomapinpath)
	if os.path.isfile(geodatabasepath):
		os.remove(geodatabasepath)

	
	if makepng:
		try:
			import f2n
		except ImportError:
			print "Couldn't import f2n -- install it !"
			return
		myimage = f2n.fromfits(alifilepath, verbose=False)
		myimage.setzscale("auto", "auto")
		myimage.makepilimage("log", negative = False)
		myimage.writetitle(os.path.basename(alifilepath))
		if not os.path.isdir(outdir):
				os.makedirs(outdir)
		myimage.tonet(os.path.join(outdir, os.path.basename(alifilepath)+".png"))


