"""
Quads are asterisms of 4 stars, used to match the catalogs.
"""

import sys, os
import math
import numpy as np
import operator # For sorting
import copy
import itertools
import scipy.spatial



class Quad:
	"""
	A geometric "hash", or asterism, as used in Astrometry.net :
	http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0910.2233
	It is made out of 4 stars, and it is shift / scale / rotation invariant
	"""
	
	def __init__(self, fourstars):
		"""
		fourstars is a list of four stars
		
		We make the following attributes :
		self.hash
		self.stars (in the order A, B, C, D)
		
		"""
		assert len(fourstars) == 4
		
		tests = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
		other = [(2,3), (1,3), (1,2), (0,3), (0,2), (0,1)]
		dists = np.array([fourstars[i].distance(fourstars[j]) for (i,j) in tests])
		assert np.min(dists) > 1.0

		maxindex = np.argmax(dists)
		(Ai, Bi) = tests[maxindex] # Indexes of stars A and B
		(Ci, Di) = other[maxindex] # Indexes of stars C and D
		A = fourstars[Ai]
		B = fourstars[Bi]
		C = fourstars[Ci]
		D = fourstars[Di]
		
		# We look for matrix transform [[a -b], [b a]] + [c d] that brings A and B to 00 11 :
		x = B.x - A.x
		y = B.y - A.y
		b = (x-y)/(x*x + y*y)
		a = (1.0/x) * (1.0 + b*y)
		c = b*A.y - a*A.x 
		d = - (b*A.x + a*A.y)
		
		t = SimpleTransform((a, b, c, d))
		
		# Test
		#print t.apply((A.x, A.y))
		#print t.apply((B.x, B.y))
		
		(xC, yC) = t.apply((C.x, C.y))
		(xD, yD) = t.apply((D.x, D.y))
		
		# Normal case
		self.hash = (xC, yC, xD, yD)
		
		# Break symmetries :
		testa = xC > xD
		testb = xC + xD > 1
		
		if testa and not testb: # we switch C and D
			#print "a"
			self.hash = (xD, yD, xC, yC)
			(C, D) = (D, C)
		
		if testb and not testa: # We switch A and B
			#print "b"
			self.hash = (1.0-xD, 1.0-yD, 1.0-xC, 1.0-yC)
			(A, B) = (B, A)
			(C, D) = (D, C)
			
		if testa and testb:
			#print "a + b"
			self.hash = (1.0-xC, 1.0-yC, 1.0-xD, 1.0-yD)
			(A, B) = (B, A)
	
		# Checks :
		assert self.hash[0] <= self.hash[2]
		assert self.hash[0] + self.hash[2] <= 1
		
		self.stars = [A, B, C, D] # Order might be different from the fourstars !
		
		
	def __str__(self):
		return "Hash : %6.3f %6.3f %6.3f %6.3f / IDs : (%s, %s, %s, %s)" % (
			self.hash[0], self.hash[1], self.hash[2], self.hash[3],
			self.stars[0].name, self.stars[1].name, self.stars[2].name, self.stars[3].name)


def mindist(fourstars):
	"""
	Function that tests if 4 stars are suitable to make a good quad...
	"""
	tests = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
	dists = np.array([fourstars[i].distance(fourstars[j]) for (i,j) in tests])
	return np.min(dists)


def ccworder(a):
	"""
	Sorting a coordinate array CCW to plot polygons ...
	"""
	ac = a - np.mean(a, 0)
	indices = np.argsort(np.arctan2(ac[:, 1], ac[:, 0]))
	return a[indices]

	

def makequads(starlist, n=10, plot=False, verbose=True):
	"""
	Give me a list, I return a list of some quads. This is the magic kitchen recipe...
	"""

	quadlist = []
	sortedstars = sortstarlistbyflux(starlist)
	
	if verbose:
		print "Building quads for %i stars ..." % len(starlist)
	
	# We start by combis of the brightest ones :
	
	for fourstars in itertools.combinations(sortedstars[:8], 4):
		if mindist(fourstars) > 50.0:
				quadlist.append(Quad(fourstars))	
	
	
	# Bright quads in subareas :
	(xmin, xmax, ymin, ymax) = area(sortedstars)
	
	f = 3
	r = 1.5*max(xmax - xmin, ymax - ymin)/f
	for xc in np.linspace(xmin, xmax, f+2):
		for yc in np.linspace(ymin, ymax, f+2):
			cstar = Star(x=xc, y=yc)
			das = cstar.distanceandsort(sortedstars[:200])
			#closest = [s["star"] for s in das[0:4]]
			brightestwithinr = sortstarlistbyflux([s["star"] for s in das if s["dist"] <= r])[0:5]
			for fourstars in itertools.combinations(brightestwithinr, 4):
				if mindist(fourstars) > 10.0:
					quadlist.append(Quad(fourstars))
		
		
	f = 6
	r = 1.5*max(xmax - xmin, ymax - ymin)/f
	for xc in np.linspace(xmin, xmax, f+2):
		for yc in np.linspace(ymin, ymax, f+2):
			cstar = Star(x=xc, y=yc)
			das = cstar.distanceandsort(sortedstars[:200])
			#closest = [s["star"] for s in das[0:4]]
			brightestwithinr = sortstarlistbyflux([s["star"] for s in das if s["dist"] <= r])[0:5]
			for fourstars in itertools.combinations(brightestwithinr, 4):
				if mindist(fourstars) > 10.0:
					quadlist.append(Quad(fourstars))

	"""
	f = 12
	r = 2.0*max(xmax - xmin, ymax - ymin)/f
	for xc in np.linspace(xmin, xmax, f+2):
		for yc in np.linspace(ymin, ymax, f+2):
			cstar = Star(x=xc, y=yc)
			das = cstar.distanceandsort(sortedstars[:200])
			#closest = [s["star"] for s in das[0:4]]
			brightestwithinr = sortstarlistbyflux([s["star"] for s in das if s["dist"] <= r])[0:4]
			for fourstars in itertools.combinations(brightestwithinr, 4):
				if mindist(fourstars) > 10.0:
					quadlist.append(Quad(fourstars))
	"""
	
	
	if verbose:
		print "Done, %i quads" % (len(quadlist))

	if plot:
		import matplotlib.pyplot as plt
		import matplotlib.patches as patches
		from matplotlib.collections import PatchCollection
	
		a = listtoarray(starlist)
		plt.plot(a[:,0], a[:,1], marker=",", ls="none", color="black")
		ax = plt.gca()

		for quad in quadlist:
			polycorners = listtoarray(quad.stars)
			polycorners = ccworder(polycorners)
			plt.fill(polycorners[:,0], polycorners[:,1], alpha=0.1, ec="none")
	
		(xmin, xmax, ymin, ymax) = area(starlist)
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		ax.set_aspect('equal', 'datalim')
		plt.show()
	
	return quadlist
	

def proposecands(uknquadlist, refquadlist, n=5):
	"""
	Function that identifies similar quads between the unknown image and a reference.
	"""

	uknhashs = np.array([quad.hash for quad in uknquadlist])	
	refhashs = np.array([quad.hash for quad in refquadlist])
	#print "Unknown quads   : ", uknhashs.shape[0]
	#print "Reference quads : ", refhashs.shape[0]
	
	# Brute force...
	dists = scipy.spatial.distance.cdist(refhashs, uknhashs)
	uknmindistindexes = np.argmin(dists, axis=0) # For each ukn, the index of the closest ref
	uknmindist = np.min(dists, axis=0) # The corresponding distances
	uknbestindexes = np.argsort(uknmindist)
	
	candlist = []
	for i in range(n):
		cand = [uknquadlist[uknbestindexes[i]], refquadlist[uknmindistindexes[uknbestindexes[i]]]]
		candlist.append(cand)
		#print uknmindist[uknbestindexes[i]]
	
	return candlist
	

def quadtrans(uknquad, refquad):
	"""
	Quickly return a transform estimated from the stars A and B of two quads.
	"""
	t = SimpleTransform()
	t.fitstars(uknquad.stars[:2], refquad.stars[:2])
	return t





