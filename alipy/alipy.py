import imgcat

	


def run(refcat, cat, refn=200, catn=100):
	"""
	ref is a reference catalog
	cats is a list of catalogs to match to the ref
	-> For each cat in cats, finds a trans to make it match to ref
	"""
	print "Making ref quads"
	refstars = star.sortstarlistbyflux(star.readsexcat(refcat, verbose=False))[:refn]
	refquads = star.makequads(refstars, plot=False)

	print "Making cat quads"
	uknstars = star.sortstarlistbyflux(star.readsexcat(cat, verbose=False))[:catn]
	uknquads = star.makequads(uknstars, plot=False)

	print "Finding cands"
	cands = star.proposecands(uknquads, refquads)
	
	print "Refining"
	for cand in cands:
		trans = star.quadtrans(*cand)
		nident = trans.teststars(uknstars, refstars)
		if nident > 10:
			#print "Match"
			trans.refinestars(uknstars, refstars)
			print trans, "->", nident, "stars"
			break
		

