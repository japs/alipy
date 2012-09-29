import imgcat
import quad
	

class Identification:
	"""
	Represents the identification of a transform between two ImgCat objects
	"""

	def __init__(self, ref, ukn):
		"""
		
		:param ref: The reference image
		:type ref: ImgCat object
		:param ukn: The unknown image, whose transform will be adjusted to match the ref
		:type ukn: ImgCat object
		
			
		"""
		self.ref = ref
		self.ukn = ukn
		
		self.ok = False
		
	
	def __str__(self):
		return "Hello"
		
	
	def run(self):
		"""
		I try to find a match, increasing the number of quads incrementally if I don't succeed.
		"""
		
		self.ref.makequadlist()
		self.ukn.makequadlist()
	

		print "Finding cands ..."
		cands = quad.proposecands(self.ukn.quadlist, self.ref.quadlist)
		
		print "Refining ..."
		for cand in cands:
			trans = quad.quadtrans(*cand)
			nident = trans.teststars(self.ukn.starlist, self.ref.starlist)
			if nident > 10:
				#print "Match"
				trans.refinestars(self.ukn.starlist, self.ref.starlist)
				self.ukn.transform = trans
				self.ok = True
				break


def run(refpath, uknpathlist, visu=True):
	"""
	refpath is a reference FITS file
	uknpathlist is a list of FITS files to align
	"""
	
	print "Preparing ref ..."
	ref = imgcat.ImgCat(refpath)
	ref.makecat(rerun=False)
	ref.makestarlist()
	if visu:
		ref.showstars()
	
	
	for uknpath in uknpathlist:
		
		ukn = imgcat.ImgCat(uknpath)
		ukn.makecat(rerun=False)
		ukn.makestarlist()
		if visu:
			ukn.showstars()

		idn = Identification(ref, ukn)
		idn.run()
		if idn.ok:
			print "Done with %s :-P" % (uknpath)
			#print ukn.transform
			
			ukn.affineremap(shape=(2000, 2000))
	



