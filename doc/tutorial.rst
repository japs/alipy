Tutorial
========

Below is a highly commented demo script. Note that all the SExtractor and catalog identification stuff happens in **one** single line. The rest of this demo is a quick illustration of what can be done with the identifications in hand. See the API doc of these function and classes (search field on the left) for detailed information !

::
		
	import alipy
	import glob
	
	images_to_align = sorted(glob.glob("images/*.fits"))
	ref_image = "ref.fits"
	
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	# That's it !
	# Put visu=True to get visualizations in form of png files (nice but much slower)
	# On multi-extension data, you will want to specify the hdu (see API doc).
	
	# The output is a list of Identification objects, which contain the transforms :
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
			# id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
			# you can directly access its parameters :
			#print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			#print id.trans.matrixform()
			#print id.trans.inverse() # this returns a new SimpleTransform object
		
		else:
			print "%20s : no transformation found !" % (id.ukn.name) 
		
	# Minimal example of how to align images :
	
	outputshape = alipy.align.shape(ref_image)
	# This is simply a tuple (width, height)... you could specify any other shape.
	
	for id in identifications:
		if id.ok == True:
		
			# Variant 1, using only scipy and the simple affine transorm :
			alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)
			
			# Variant 2, using geomap/gregister, correcting also for distortions :
			alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)
			# id.uknmatchstars and id.refmatchstars are simply lists of corresponding Star objects.
			
			# By default, the aligned images are written into a directory "alipy_out".
	
	# To be continued ...

			
The important functions and classes (links take you to the API documentation) :
 * :py:func:`alipy.ident.run` : the function that returns the :py:class:`~alipy.ident.Identification` objects.
 * :py:class:`alipy.ident.Identification` : the objects returned by the above :py:func:`~alipy.ident.run`. Note that these objects also contain lists of the matched stars.
 * :py:class:`alipy.star.Star`
 * :py:class:`alipy.star.SimpleTransform`
 * :py:func:`alipy.align.affineremap`
 * :py:func:`alipy.align.irafalign`

