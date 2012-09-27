import star


def findtrans(refstars, refquads, cat):
	

def run(refcat, cats, refn=200, catn=100):
	"""
	ref is a reference catalog
	cats is a list of catalogs to match to the ref
	-> For each cat in cats, finds a trans to make it match to ref
	"""
	refstars = star.sortstarlistbyflux(star.readsexcat(refcat, verbose=False))[:refn]
	refquads = star.makequads(refstars, n=15)

	for cat in cats:
		
		uknstars = star.sortstarlistbyflux(star.readsexcat(cat, verbose=False))[:catn]
		uknquads = star.makequads(uknstars, n=15)

		cands = star.proposecands(uknquads, refquads, n=10)

		for cand in cands:
			trans = star.quadtrans(*cand)
			nident = trans.teststars(uknstars, refstars) 
			if nident > catn/3.0:
				#print "Match"
				trans.refinestars(uknstars, refstars)
				print trans
				break
			

