Tutorial
========

Minimalistic for now ...

::
		
	import alipy
	import glob
	
	images_to_align = sorted(glob.glob("images/*.fits"))
	ref_image = "ref.fits"
	
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	# That's it -- all in one line.
	# Put visu=True to get visualizations (nice but much slower)
	
	# The output contains the transforms :
	for id in identifications: # identifications is a lists of objects as long as to_align_images.
		if id.ok == True: # i.e., if it worked
			
			print "%20s %20s" % (id.ukn.name, id.trans) # This is a alipy.star.SimpleTransform object
			#print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			#print id.trans.matrixform()
			#print id.trans.inverse() # this returns a new SimpleTransform object
	
	
	# Minimal example of how to align images :
	for id in identifications:
		if id.ok == True:
			alipy.align.affineremap(id.ukn.filepath, id.trans, shape=(1500,1500), makepng=True)
			

			
The most important functions (just testing sphinx links ...) :
 * :py:func:`alipy.ident.run`
 * :py:class:`alipy.star.SimpleTransform`


