import imgcat
import quad
import star
import sys
import os
	

class Identification:
	"""
	Represents the identification of a transform between two ImgCat objects.
	Regroups all the star catalogs, the transform, the quads, the candidate, etc.
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
		
		self.trans = None
		self.uknmatchstars = []
		self.refmatchstars = []
		self.cand = None
		

		
	def findtrans(self, r = 5.0, verbose=True):
		"""
		Find the best trans given the quads, and tests if the match is sufficient	
		"""
		
		# First question : how many stars should match ?
		if len(self.ukn.starlist) < 5: # Then we should simply try to get the smallest distance...
			minnident = 4
		else:
			minnident = max(4, min(8, len(self.ukn.starlist)/5.0)) # Perfectly arbitrary, let's see how it works
		
		# Hmm, arbitrary for now :
		minquaddist = 0.001
		
		# Let's start :
		if self.ref.quadlevel == 0:
			self.ref.makemorequads(verbose=verbose)
		if self.ukn.quadlevel == 0:
			self.ukn.makemorequads(verbose=verbose)
			
		while self.ok == False:
			
			# Find the best candidates
			cands = quad.proposecands(self.ukn.quadlist, self.ref.quadlist, n=4, verbose=verbose)
			if cands[0]["dist"] < minquaddist:
				for cand in cands:
					# Check how many stars are identified...					
					nident = star.identify(self.ukn.starlist, self.ref.starlist, trans=cand["trans"], r=r, verbose=verbose, getstars=False)
					if nident >= minnident:
						self.trans = cand["trans"]
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
			(self.uknmatchstars, self.refmatchstars) = star.identify(self.ukn.starlist, self.ref.starlist, trans=self.trans, r=r, verbose=False, getstars=True)
			# refit the transform on them :
			if verbose:
				print "Refitting transform (before/after) :"
				print self.trans
			self.trans = star.fitstars(self.uknmatchstars, self.refmatchstars)
			if verbose:
				print self.trans
			# Generating final matched star lists :
			(self.uknmatchstars, self.refmatchstars) = star.identify(self.ukn.starlist, self.ref.starlist, trans=self.trans, r=r, verbose=verbose, getstars=True)

			if verbose:
				print "I'm done !"
	
		else:
			if verbose:
				print "Failed to find transform !"
			
			
	
	def showmatch(self, show=False, verbose=True):
		"""
		A plot of the transformed stars and the candidate quad
		"""
		if self.ok == False:
			return
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
		a = star.listtoarray(self.trans.applystarlist(self.ukn.starlist), full=True)
		plt.scatter(a[:,0], a[:,1], s=2.0, color="red")
		a = star.listtoarray(self.trans.applystarlist(self.uknmatchstars), full=True)
		plt.scatter(a[:,0], a[:,1], s=6.0, color="red")
		
		# The quad
		
		polycorners = star.listtoarray(self.cand["refquad"].stars)
		polycorners = imgcat.ccworder(polycorners)
		plt.fill(polycorners[:,0], polycorners[:,1], alpha=0.1, ec="none", color="red")

		plt.xlim(self.ref.xlim)
		plt.ylim(self.ref.ylim)
		plt.title("Match of %s" % (str(self.ukn.name)))
		plt.xlabel("ref x")
		plt.ylabel("ref y")
		ax = plt.gca()
		ax.set_aspect('equal', 'datalim')
	
		if show:
			plt.show()
		else:
			if not os.path.isdir("alipy_visu"):
				os.makedirs("alipy_visu")
			plt.savefig(os.path.join("alipy_visu", self.ukn.name + "_match.png"))




def run(ref, ukns, visu=True, skipsaturated=False, r = 5.0, verbose=True):
	"""
	Top-level function to identify transorms between images.
	Returns a list of alipy.Identification objects that contain all the info to go further.
	TODO : MAKE THIS GUY ACCEPT EXISTING ASCIIDATA CATALOGS

	:param ref: path to FITS file that acts as the "reference".
	:type ref: string
	
	:param ukns: list of paths to FITS files to be "aligned" on the reference. **ukn** stands for unknown ...
	:type ref: list of strings
	
	:param visu: If yes, I'll draw some visualizations of the process (good to understand problems ...)
	:type visu: boolean
	
	:param skipsaturated: Should I skip saturated stars ?
	:type skipsaturated: boolean
	
	:param r: Identification radius in pixels of the reference image (default 5.0 should be fine).
	:type r: float
	
	"""
	
	if verbose:
		print 10*"#", " Preparing reference ..."
	ref = imgcat.ImgCat(ref)
	ref.makecat(rerun=False, verbose=verbose)
	ref.makestarlist(skipsaturated=skipsaturated, verbose=verbose)
	if visu:
		ref.showstars(verbose=verbose)
	ref.makemorequads(verbose=verbose)
	
	identifications = []
	
	for ukn in ukns:
		
		if verbose:
			print 10*"#", "Processing %s" % (ukn)
		
		ukn = imgcat.ImgCat(ukn)
		ukn.makecat(rerun=False, verbose=verbose)
		ukn.makestarlist(skipsaturated=skipsaturated, verbose=verbose)
		if visu:
			ukn.showstars(verbose=verbose)

		idn = Identification(ref, ukn)
		idn.findtrans(verbose=verbose, r=r)
		identifications.append(idn)
		
		if visu:
			ukn.showquads(verbose=verbose)
			idn.showmatch(verbose=verbose)
		
	if visu:
		ref.showquads(verbose=verbose)
		
	return identifications



