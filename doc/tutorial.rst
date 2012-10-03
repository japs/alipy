Tutorial
========

Minimalistic for now; below is a simple commented demo script. See the API doc of all these functions (search field on the left) for more information !

::
		
	import alipy
	import glob
	
	images_to_align = sorted(glob.glob("images/*.fits"))
	ref_image = "ref.fits"
	
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	# That's it -- all in one line.
	# Put visu=True to get visualizations in form of png files (nice but much slower)
	
	# The output contains the transforms :
	for id in identifications: # identifications is a list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			
			print "%20s %20s" % (id.ukn.name, id.trans) # This is a alipy.star.SimpleTransform object
			#print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			#print id.trans.matrixform()
			#print id.trans.inverse() # this returns a new SimpleTransform object
	
	
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

			
The most important functions (just testing sphinx links ...) :
 * :py:func:`alipy.ident.run`
 * :py:func:`alipy.align.affineremap`
 * :py:func:`alipy.align.irafalign`
 * :py:class:`alipy.star.Star`
 * :py:class:`alipy.star.SimpleTransform`


