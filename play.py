
import pysex
import glob
#import numpy as np
#import matplotlib.pyplot as plt
#import star

import alipy


imagepaths = sorted(glob.glob("images/*.fits"))
print "\n".join(imagepaths)

cats = [pysex.run(ip, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'FWHM_IMAGE', 'FLAGS', 'ELONGATION', 'NUMBER'], keepcat=True, rerun=False, catdir="cats") for ip in imagepaths]

#refcat = cats[1]

#alipy.run(refcat, cats[0])
#alipy.run(refcat, cats[2])
#alipy.run(refcat, cats[3])

refcat = cats[0]
alipy.run(refcat, cats[1])
alipy.run(refcat, cats[2])
alipy.run(refcat, cats[3])


exit()

#refstars = star.readsexcat(refcat, verbose=False)
#uknstars = star.readsexcat(ukncat, verbose=False)

"""
fourstars = star.sortstarlistbyflux(stars)[0:4]
print star.Quad(fourstars)
print star.Quad(fourstars[::-1])
print star.Quad([fourstars[0], fourstars[1], fourstars[3], fourstars[2]])
print star.Quad([fourstars[1], fourstars[0], fourstars[2], fourstars[3]])
print star.Quad([fourstars[0], fourstars[2], fourstars[1], fourstars[3]])
print star.Quad([fourstars[2], fourstars[1], fourstars[3], fourstars[0]])
"""


#stars = star.sortstarlistbyflux(stars)[0:100]

uknquads = star.makequads(uknstars, n=10, plot=False)
refquads = star.makequads(refstars, n=15, plot=False)

cands = star.proposecands(uknquads, refquads)

for cand in cands:
	trans = star.quadtrans(*cand)
	#print trans
	#print trans.teststars(uknstars, refstars)
	print trans
	print trans.teststars(uknstars, refstars)
	trans.refinestars(uknstars, refstars)
	print trans
	print trans.teststars(uknstars, refstars)

#star.checktrans(uknstars, refstars, trans)
#star.checktrans(uknquad.stars, refquad.stars, trans)



#import matplotlib.pyplot as plt





#for s in fourstars:
#	#s.x *= 2.0
#	s.y *= -1.0

#print star.Quad(fourstars)

"""
plt.scatter(a[:,0], a[:,1])
plt.show()

plt.figure(figsize=(7, 7))
plt.scatter([q.hash[0], q.hash[2]], [q.hash[1], q.hash[3]])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.show()
"""


#print len(a)
#plt.scatter(a[:,0], a[:,1])
#plt.show()



def sourcearray(cat):
	"""
	Puts the source catalog into a 2D numpy array for fast manipulations.
	x y id
	"""

	#print cat.info()
	
	xs = cat["X_IMAGE"].tonumpy()
	ys = cat["Y_IMAGE"].tonumpy()
	nums = cat["NUMBER"].tonumpy()
	
	a = np.vstack((xs, ys, nums))#.transpose()
	return a
	
#a = sourcearray(cat)
