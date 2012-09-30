Tutorial
========

Minimalistic for now ...

::
	
	import glob
	import alipy

	to_align_images = sorted(glob.glob("images/*.fits"))
	ref_image = "myref.fits"
	
	transforms = alipy.ident.run(ref_image, to_align_images, visu=False)
	# Put visu=True to get visualizations (nice but much slower)

	for transform in transforms: # transforms is a lists of SimpleTransform objects
		if transform != None: # i.e., if it could be identified :
			print transform
			print transform.v # [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			print transform.matrixform()
			print transform.inverse() # this is a new SimpleTransform object
	
