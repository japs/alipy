import imgcat
import quad
import star
import sys
import os
	

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
		
		
		self.uknmatchstars = []
		self.refmatchstars = []
		self.cand = None
		

		
	def findtrans(self, r = 5.0, verbose=True):
		"""
		Find the best trans given the quads, and tests if the match is sufficient	
		"""
		
		# First question : what is a good match ?
		if len(self.ukn.starlist) < 5: # Then we should simply try to get the smallest distance...
			minnident = 4
		else:
			minnident = max(4, min(8, len(self.ukn.starlist)/5.0)) # Perfectly arbitrary, let's see how it works
			
		
		# Let's start :
		if self.ref.quadlevel == 0:
			self.ref.makemorequads(verbose=verbose)
		if self.ukn.quadlevel == 0:
			self.ukn.makemorequads(verbose=verbose)
			
		while self.ok == False:
			
			# Find the best candidates
			cands = quad.proposecands(self.ukn.quadlist, self.ref.quadlist, n=4, verbose=verbose)
			for cand in cands:
				# Check how many stars are identified...
				
				nident = star.identify(self.ukn.starlist, self.ref.starlist, trans=cand["trans"], r=r, verbose=verbose, getstars=False)
				if nident >= minnident:
					self.ukn.transform = cand["trans"]
					self.cand = cand
					self.ok = True
					break # get out of the for
					
			if self.ok == False:
				# We add more quads...
				addedmorerefquads = self.ref.makemorequads(verbose=verbose)
				addedmoreuknquads = self.ukn.makemorequads(verbose=verbose)
				
				if addedmorerefquads == False and addedmoreuknquads == False:
					break # get out of the while, we failed.
	

		if self.ok: # we refine the transform
			# get matching stars :
			(self.uknmatchstars, self.refmatchstars) = star.identify(self.ukn.starlist, self.ref.starlist, trans=self.ukn.transform, r=r, verbose=False, getstars=True)
			# refit the transform on them :
			if verbose:
				print "Refitting transform (before/after) :"
				print self.ukn.transform
			self.ukn.transform = star.fitstars(self.uknmatchstars, self.refmatchstars)
			if verbose:
				print self.ukn.transform
			# Generating final matched star lists :
			(self.uknmatchstars, self.refmatchstars) = star.identify(self.ukn.starlist, self.ref.starlist, trans=self.ukn.transform, r=r, verbose=verbose, getstars=True)

			if verbose:
				print "I'm done !"
	
		else:
			if verbose:
				print "Failed to find transform !"
			
			
	
	def showmatch(self, show=False, verbose=True):
		"""
		A plot of the transformed stars and the candidate quad
		"""
		if verbose:
			print "Plotting match ..."
		import matplotlib.pyplot as plt
		#import matplotlib.patches
		#import matplotlib.collections
		
		plt.figure(figsize=(10, 10))
		
		# The ref in black
		a = star.listtoarray(self.ref.starlist, full=True)
		plt.scatter(a[:,0], a[:,1], s=2.0, color="black")
		a = star.listtoarray(self.refmatchstars, full=True)
		plt.scatter(a[:,0], a[:,1], s=10.0, color="black")
		
		# The ukn in red
		a = star.listtoarray(self.ukn.transform.applystarlist(self.ukn.starlist), full=True)
		plt.scatter(a[:,0], a[:,1], s=2.0, color="red")
		a = star.listtoarray(self.ukn.transform.applystarlist(self.uknmatchstars), full=True)
		plt.scatter(a[:,0], a[:,1], s=6.0, color="red")
		
		# The quad
		
		polycorners = star.listtoarray(self.cand["refquad"].stars)
		polycorners = imgcat.ccworder(polycorners)
		plt.fill(polycorners[:,0], polycorners[:,1], alpha=0.1, ec="none", color="red")

		plt.xlim(self.ref.xlim)
		plt.ylim(self.ref.ylim)
		plt.title("Match of %s" % (str(self.ukn.common)))
		plt.xlabel("ref x")
		plt.ylabel("ref y")
		ax = plt.gca()
		ax.set_aspect('equal', 'datalim')
	
		if show:
			plt.show()
		else:
			if not os.path.isdir("alipy_visu"):
				os.makedirs("alipy_visu")
			plt.savefig(os.path.join("alipy_visu", self.ukn.common + "_match.png"))


def run(refpath, uknpathlist, visu=True, skipsaturated=False, r = 5.0, verbose=True):
	"""
	refpath is a reference FITS file / asciidata catalog TODO
	uknpathlist is a list of FITS files to align
	
	:param skipsaturated: Should I skip saturated stars ?
	:type skipsaturated: boolean
	
	:param r: Identification radius in pixels of the reference image (default 5.0 should be fine).
	:type r: float
	
	"""
	
	if verbose:
		print 10*"#", " Preparing reference ..."
	ref = imgcat.ImgCat(refpath)
	ref.makecat(rerun=False, verbose=verbose)
	ref.makestarlist(skipsaturated=skipsaturated, verbose=verbose)
	if visu:
		ref.showstars(verbose=verbose)
	ref.makemorequads(verbose=verbose)
	
	
	transforms = []
	
	for uknpath in uknpathlist:
		
		if verbose:
			print 10*"#", "Processing %s" % (uknpath)
		
		ukn = imgcat.ImgCat(uknpath)
		ukn.makecat(rerun=False, verbose=verbose)
		ukn.makestarlist(skipsaturated=skipsaturated, verbose=verbose)
		if visu:
			ukn.showstars(verbose=verbose)

		idn = Identification(ref, ukn)
		idn.findtrans(verbose=verbose, r=r)
		
		if visu:
			ukn.showquads(verbose=verbose)
			idn.showmatch(verbose=verbose)
		
		if idn.ok:
			transforms.append(ukn.transform)
		else:
			transforms.append(None)
				
		#ukn.affineremap(shape=(2000, 2000), makepng=True)
	if visu:
		ref.showquads(verbose=verbose)
		
	
	return transforms



