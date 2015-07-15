"""
Overhaul of cosmouline's star module, for alipy2.
This module contains stuff for geometric matching algorithms.
"""

import sys
import os
import math
import numpy as np
import operator  # For sorting
import copy
import itertools
import scipy.linalg
import scipy.spatial


class Star:
    """
    Simple class to represent a single source
    (usually stars, but not necessarily).
    In this module we often manipulate lists of such Star objects.
    """

    def __init__(self, x=0.0, y=0.0, name="untitled",
                 flux=-1.0, props={}, fwhm=-1.0, elon=-1.0):
        """
        flux : Some "default" or "automatic" flux, might be a just good guess.
               Used for sorting etc.
               If you have several fluxes, colours, store them in the
               props dict.
        props : A placeholder dict to contain other properties of your choice
                (not required nor used by the methods).
        """
        self.x = float(x)
        self.y = float(y)
        self.name = str(name)
        self.flux = float(flux)
        self.props = props
        self.fwhm = float(fwhm)
        self.elon = float(elon)

    def copy(self):
        return copy.deepcopy(self)

    def __getitem__(self, key):
        """
        Used for sorting list of stars.
        """
        if key == 'flux':
            return self.flux
        if key == 'fwhm':
            return self.fwhm
        if key == 'elon':
            return self.elon

    def __str__(self):
        """
        A string representation of a source.
        """
        return "%10s : (%8.2f,%8.2f) | %12.2f | %5.2f %5.2f" % (self.name,
                                                                self.x,
                                                                self.y,
                                                                self.flux,
                                                                self.fwhm,
                                                                self.elon)

    def coords(self, full=False):
        """
        Returns the coords in form of an array.

        :param full: If True, I also include flux, fwhm, elon
        :type full: boolean

        """
        if full:
            return np.array([self.x, self.y, self.flux, self.fwhm, self.elon])
        else:
            return np.array([self.x, self.y])

    def distance(self, otherstar):
        """
        Returns the distance between the two stars.
        """
        return math.sqrt(np.sum((self.coords() - otherstar.coords()) ** 2))

    def trigangle(self, otherstar):
        """
        returns the "trigonometric" angle of the vector that goes from
        self to the otherstar, in degrees
        """
        return math.atan2(otherstar.y - self.y, otherstar.x - self.x) * \
               (180.0 / math.pi) % 360.0

    def distanceandsort(self, otherstarlist):
        """
        Returns a list of dicts(star, dist, origpos),
        sorted by distance to self. The 0th star is the closest.

        otherstarlist is not modified.
        """
        import operator  # for the sorting

        returnlist = []
        for i, star in enumerate(otherstarlist):
            dist = self.distance(star)
            returnlist.append({'star': star, 'dist': dist, 'origpos': i})
        returnlist = sorted(returnlist, key=operator.itemgetter('dist'))
                            # sort stars according to dist

        return returnlist

# And now some functions to manipulate list of such stars ###


def printlist(starlist):
    """
    Prints the stars ...
    """
    for source in starlist:
        print(source)


def listtoarray(starlist, full=False):
    """
    Transforms the starlist into a 2D numpy array for fast manipulations.
    First index is star, second index is x or y

    :param full: If True, I also include flux, fwhm, elon
    :type full: boolean

    """
    return np.array([star.coords(full=full) for star in starlist])


def area(starlist, border=0.01):
    """
    Returns the area covered by the stars.
    Border is relative to max-min
    """
    if len(starlist) == 0:
        return np.array([0, 1, 0, 1])

    if len(starlist) == 1:
        star = starlist[0]
        return np.array([star.x - 0.5, star.x + 0.5,
                         star.y - 0.5, star.y + 0.5])

    a = listtoarray(starlist)
    (xmin, xmax) = (np.min(a[:, 0]), np.max(a[:, 0]))
    (ymin, ymax) = (np.min(a[:, 1]), np.max(a[:, 1]))
    xw = xmax - xmin
    yw = ymax - ymin
    xmin = xmin - border * xw
    xmax = xmax + border * xw
    ymin = ymin - border * yw
    ymax = ymax + border * yw
    return np.array([xmin, xmax, ymin, ymax])


def readmancat(mancatfilepath, verbose="True"):
    """
    Reads a "manual" star catalog --
    by manual, I mean "not written by sextractor", so this is typically a 
    *short* file.

    Comment lines start with #, blank lines are ignored.
    The format of a data line is

    starname xpos ypos [flux]

    The data is returned as a list of star objects.
    """

    if not os.path.isfile(mancatfilepath):
        print("File does not exist :")
        print(mancatfilepath)
        print("Line format to write : starname xpos ypos [flux]")
        sys.exit(1)

    myfile = open(mancatfilepath, "r")
    lines = myfile.readlines()
    myfile.close

    table = []
    knownnames = []  # We check for uniqueness of the names

    for i, line in enumerate(lines):
        if line[0] == '#' or len(line) < 4:
            continue
        elements = line.split()
        nbelements = len(elements)

        if nbelements != 3 and nbelements != 4:
            print("Format error on line", i + 1, "of :")
            print(mancatfilepath)
            print("The line looks like this :")
            print(line)
            print("... but we want : starname xpos ypos [flux]")
            sys.exit(1)

        name = elements[0]
        x = float(elements[1])
        y = float(elements[2])
        if nbelements == 4:
            flux = float(elements[3])
        else:
            flux = -1.0

        if name in knownnames:
            print("Error in %s" % (mancatfilepath))
            print("The name '%s' (line %i) is already taken." % (name, i + 1))
            print("This is insane, bye !")
            sys.exit(1)
        knownnames.append(name)

        # table.append({"name":name, "x":x, "y":y, "flux":flux})
        table.append(Star(x=x, y=y, name=name, flux=flux))

    if verbose:
        print("I've read", len(table), "sources from",
              os.path.split(mancatfilepath)[1])
    return table


def readsexcat(sexcat, hdu=0, verbose=True,
               maxflag=3, posflux=True, minfwhm=2.0, propfields=[]):
    """
    sexcat is either a string (path to a file),
    or directly an asciidata catalog object as returned by pysex

    :param hdu: The hdu containing the science data from which I should build
                the catalog. 0 will select the only available extension.
                If multihdu, 1 is usually science.

    We read a sextractor catalog with astroasciidata and return a list
    of stars.
    Minimal fields that must be present in the catalog :

        * NUMBER
        * EXT_NUMBER
        * X_IMAGE
        * Y_IMAGE
        * FWHM_IMAGE
        * ELONGATION
        * FLUX_AUTO
        * FLAGS

    maxflag : maximum value of the FLAGS that you still want to keep.
              Sources with higher values will be skipped.
        * FLAGS == 0 : all is fine
        * FLAGS == 2 : the flux is blended with another one;
                       further info in the sextractor manual.
        * FLAGS == 4    At least one pixel of the object is saturated
                        (or very close to)
        * FLAGS == 8    The object is truncated
                        (too close to an image boundary)
        * FLAGS is the sum of these ...

    posflux : if True, only stars with positive FLUX_AUTO are included.

    propfields : list of FIELD NAMES to be added to the props of the stars.

    I will always add FLAGS as a propfield by default.
    """
    returnlist = []

    if isinstance(sexcat, str):

        import asciidata
        if not os.path.isfile(sexcat):
            print("Sextractor catalog does not exist :")
            print(sexcat)
            sys.exit(1)

        if verbose:
            print("Reading %s " % (os.path.split(sexcat)[1]))
        mycat = asciidata.open(sexcat)

    else:  # then it's already a asciidata object
        mycat = sexcat

    # We check for the presence of required fields :
    minimalfields = ["NUMBER", "X_IMAGE", "Y_IMAGE", "FWHM_IMAGE",
                     "ELONGATION", "FLUX_AUTO", "FLAGS", "EXT_NUMBER"]
    minimalfields.extend(propfields)
    availablefields = [col.colname for col in mycat]
    for field in minimalfields:
        if field not in availablefields:
            print("Field %s not available in your catalog file !" % (field))
            sys.exit(1)

    if verbose:
        print("Number of sources in catalog : %i" % (mycat.nrows))

    extnumbers = np.unique(mycat['EXT_NUMBER'].tonumpy())
    if verbose:
        print("EXT_NUMBER values found in catalog : %s" %
              (", ".join(["%i" % val for val in extnumbers])))

    if len(extnumbers) > 1 and hdu == 0:
        print(("Looks like we have several FITS extensions. "
               "You have to specify which hdu to use !"))
        sys.exit(1)

    propfields.append("FLAGS")
    propfields = list(set(propfields))

    if mycat.nrows == 0:
        if verbose:
            print("No stars in the catalog :-(")
    else:
        for i, num in enumerate(mycat['NUMBER']):
            if mycat['FLAGS'][i] > maxflag:
                continue
            if hdu != 0 and mycat['EXT_NUMBER'][i] != hdu:
                continue
            flux = mycat['FLUX_AUTO'][i]
            if posflux and (flux < 0.0):
                continue
            fwhm = mycat['FWHM_IMAGE'][i]
            if float(fwhm) <= minfwhm:
                continue

            props = dict([[propfield, mycat[propfield][i]]
                         for propfield in propfields])

            newstar = Star(
                x=mycat['X_IMAGE'][i],
                y=mycat['Y_IMAGE'][i],
                name=str(num),
                flux=flux,
                props=props,
                fwhm=mycat['FWHM_IMAGE'][i],
                elon=mycat['ELONGATION'][i])
            returnlist.append(newstar)

    if verbose:
        print("I've selected %i sources" % (len(returnlist)))

    return returnlist


def findstar(starlist, nametofind):
    """
    Returns a list of stars for which name == nametofind
    """
    foundstars = []
    for source in starlist:
        if source.name == nametofind:
            foundstars.append(source)
    return foundstars


def sortstarlistbyflux(starlist):
    """
    We sort starlist according to flux : highest flux first !
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter('flux'))
    sortedstarlist.reverse()
    return sortedstarlist


def sortstarlistby(starlist, measure):
    """
    We sort starlist according to measure : lowest first !
    Where measure is one of flux, fwhm, elon
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter(measure))
    return sortedstarlist


class SimpleTransform:

    """
    Represents an affine transformation consisting of rotation,
    isotropic scaling, and shift.
    [x', y'] = [[a -b], [b a]] * [x, y] + [c d]
    """

    def __init__(self, v=(1, 0, 0, 0)):
        """
        v = (a, b, c, d)
        """
        self.v = np.asarray(v)

    def getscaling(self):
        return math.sqrt(self.v[0] * self.v[0] + self.v[1] * self.v[1])

    def getrotation(self):
        """
        The CCW rotation angle, in degrees
        """
        return math.atan2(self.v[1], self.v[0]) * (180.0 / math.pi)  # % 360.0

    def __str__(self):
        return "Rotation %+11.6f [deg], scale %8.6f" % (self.getrotation(),
                                                        self.getscaling())

    def inverse(self):
        """
        Returns the inverse transform !
        """

        # To represent affine transformations with matrices, we can use
        # homogeneous coordinates.
        homo = np.array([
                        [self.v[0], -self.v[1], self.v[2]],
                        [self.v[1],  self.v[0], self.v[3]],
                        [0.0, 0.0, 1.0]
                        ])

        inv = scipy.linalg.inv(homo)
        # print inv

        return SimpleTransform((inv[0, 0], inv[1, 0], inv[0, 2], inv[1, 2]))

    def matrixform(self):
        """
        Special output for scipy.ndimage.interpolation.affine_transform
        Returns (matrix, offset)
        """

        return (np.array([[self.v[0], -self.v[1]],
                          [self.v[1], self.v[0]]]),
                self.v[2:4])

    def apply(self, xy):
        """
        Applies the transform to a point (x, y)
        """
        x, y = xy  # backward compatibility hack, JN
        xn = self.v[0] * x - self.v[1] * y + self.v[2]
        yn = self.v[1] * x + self.v[0] * y + self.v[3]
        return (xn, yn)

    def applystar(self, star):
        transstar = star.copy()
        (transstar.x, transstar.y) = self.apply((transstar.x, transstar.y))
        return transstar

    def applystarlist(self, starlist):
        return [self.applystar(star) for star in starlist]


def fitstars(uknstars, refstars, verbose=True):
    """
    I return the transform that puts the unknown stars (uknstars)
    onto the refstars.
    If you supply only two stars, this is using linalg.solve() --
    perfect solution.
    If you supply more stars, we use linear least squares,
    i.e. minimize the 2D error.

    Formalism inspired by:
    http://math.stackexchange.com/questions/77462/
    """

    assert len(uknstars) == len(refstars)
    if len(uknstars) < 2:
        if verbose:
            print("Sorry I cannot fit a transform on less than 2 stars.")
        return None

    # ukn * x = ref
    # x is the transform (a, b, c, d)

    ref = np.hstack(listtoarray(refstars))  # a 1D vector of lenth 2n

    uknlist = []
    for star in uknstars:
        uknlist.append([star.x, -star.y, 1, 0])
        uknlist.append([star.y, star.x, 0, 1])
    ukn = np.vstack(np.array(uknlist))  # a matrix

    if len(uknstars) == 2:
        trans = scipy.linalg.solve(ukn, ref)
    else:
        trans = scipy.linalg.lstsq(ukn, ref)[0]

    return SimpleTransform(np.asarray(trans))


#     def teststars(self, uknstars, refstars, r=5.0, verbose=True):
#         """
#         We apply the trans to the uknstarlist, and check for correspondance with the refstarlist.
#         Returns the number of uknstars that could be matched to refstars within r [unit : reference image pixels !].
#         """
#
#         transuknstars = self.applystarlist(uknstars)
#
#         transukn = listtoarray(transuknstars)
#         ref = listtoarray(refstars)
# print "Unknown stars   : ", transukn.shape[0]
# print "Reference stars : ", ref.shape[0]
#
#         mindists = np.min(scipy.spatial.distance.cdist(ref, transukn), axis=0)
#         print " must become smarter !!!!"
#         nbmatch = np.sum(mindists < r)
#         if verbose:
#             print "Tested match on %4i references : OK for %4i/%4i unknown stars (r = %.1f)." % (len(refstars), nbmatch, len(uknstars), r)
#
#         return nbmatch
#
#
#     def refinestars(self, uknstars, refstars, r=5.0, verbose=True):
#         """
#         I refit myself to all matching stars.
#         """
#
#         transuknstars = self.applystarlist(uknstars)
#         transukn = listtoarray(transuknstars)
#         ref = listtoarray(refstars)
#
# Brute force...
#         dists = scipy.spatial.distance.cdist(ref, transukn)
# uknmindistindexes = np.argmin(dists, axis=0) # For each ukn, the index of the closest ref
# uknmindist = np.min(dists, axis=0) # The corresponding distances
#         uknkeepers = uknmindist < r
#
#         matchuknstars = []
#         matchrefstars = []
#         for i in range(len(uknkeepers)):
#             if uknkeepers[i] == True:
#                 matchuknstars.append(uknstars[i])
#                 matchrefstars.append(refstars[uknmindistindexes[i]])
#         if verbose:
#             print "Refining (before / after) :"
#             print self
#         self.fitstars(matchuknstars, matchrefstars)
#         if verbose:
#             print self
#


def identify(uknstars, refstars,
             trans=None, r=5.0, verbose=True, getstars=False):
    """
    Allows to:
     * get the number or matches, i.e. evaluate the quality of the trans
     * get corresponding stars from both lists (without the transform applied)

    :param getstars: If True, I return two lists of corresponding stars,
                     instead of just the number of matching stars
    :type getstars: boolean

    Inspired by the "formpairs" of alipy 1.0 ...
    """

    if trans != None:
        ukn = listtoarray(trans.applystarlist(uknstars))
    else:
        ukn = listtoarray(uknstars)
    ref = listtoarray(refstars)

    dists = scipy.spatial.distance.cdist(
        ukn, ref)  # Big table of distances between ukn and ref
    mindists = np.min(dists, axis=1)  # For each ukn, the minimal distance
    minok = mindists <= r  # booleans for each ukn
    minokindexes = np.argwhere(
        minok).flatten()  # indexes of uknstars with matches

    if verbose:
        print(("%i/%i stars with distance < r "
               "= %.1f (mean %.1f, "
               "median %.1f, std %.1f)") % (np.sum(minok), len(uknstars), r,
                                            np.mean(mindists[minok]),
                                            np.median(mindists[minok]),
                                            np.std(mindists[minok])))
    matchuknstars = []
    matchrefstars = []

    for i in minokindexes:  # we look for the second nearest ...
        sortedrefs = np.argsort(dists[i, :])
        firstdist = dists[i, sortedrefs[0]]
        seconddist = dists[i, sortedrefs[1]]
        if seconddist > 2.0 * firstdist:  # Then the situation is clear,
                                          # we keep it.
            matchuknstars.append(uknstars[i])
            matchrefstars.append(refstars[sortedrefs[0]])
        else:
            pass  # Then there is a companion, we skip it.

    if verbose:
        print("Filtered for companions, keeping %i/%i matches" %
              (len(matchuknstars), np.sum(minok)))

    if getstars == True:
        return (matchuknstars, matchrefstars)
    else:
        return len(matchuknstars)


#

#
#
#
#
# def formpairs(starlist1, starlist2, tolerance = 2.0, onlysingle = False, transform = False, scalingratio = 1.0, angle = 0.0, shift = (0.0, 0.0), verbose = True):
#     """
#     starlist1 and starlist2 are two lists of stars.
#     For each star in starlist1, we find the closest star of starlist2, if found within a given tolerance.
#     starlist1 = hand picked stars
#     starlist2 = large catalog of
#     We return a list of pairs of the corresponding stars (in form of a dict). See first lines to get that.
#
#     transform == True :
#         starlist2 is tranformed, using scalingration, angle, and shift, prior to the pair creation.
#         Nevertheless, the idlist 2 will be filled with the raw untransformed stars from starlist2 !!!
#
#
#     tolerance = maximum distance between identified stars. Set it high -> simply select the closest one.
#     onlysingle == False : closest star within tolerance
#     onlysingle == True : same, but only if no other star is within tolerance
#
#     """
# idlist1 = [] # Stars of starlist1 identified in starlist2
# idlist2 = [] # The corresponding starlist2 stars (same number, same order)
# iddists = [] # The list of distances between the stars of idlist1 and idlist2 (same number and order)
# nomatch = [] # Stars of starlist1 that could not be identified in starlist2
# notsure = [] # Stars of starlist1 that could not doubtlessly be identified in starlist2
#
#
# If required, we transform the starlist 2 :
#     if transform :
#         transtarlist2 = copy.deepcopy(starlist2)
#         zoomstarlist(transtarlist2, scalingratio)
#         rotatestarlist(transtarlist2, angle, (0, 0))
#         shiftstarlist(transtarlist2, shift)
# Remember : we will pick the stars to fill idlist2 from the raw starlist2 !
#
#     else:
#         transtarlist2 = starlist2
#
#     returndict = {"idlist1":idlist1, "idlist2":idlist2, "iddists":iddists, "nomatch":nomatch, "notsure":notsure}
#
#     if len(starlist1) == 0:
#         if verbose :
#             print "Your starlist1 is empty, nothing to do."
#         return returndict
#
#     if len(transtarlist2) == 0:
#         if verbose :
#             print "Your starlist2 is empty, no stars to identify."
#         nomatch.extend(starlist1)
#         return returndict
#
# Special treatment in the case there is only one star in starlist2
#     if len(transtarlist2) == 1:
#         if verbose :
#             print "Your starlist2 is quite small..."
#         for handstar in starlist1:
#             closest = handstar.distanceandsort(transtarlist2)
#             if closest[0]['dist'] > tolerance:
#                 if verbose :
#                     print "No match for star %s" % handstar.name
#                 nomatch.append(handstar)
#                 continue
#             else:
#                 idlist1.append(handstar)
#                 idlist2.append(starlist2[closest[0]['origpos']])
#                 iddists.append(closest[0]['dist'])
#
#         return returndict
#
# The usual case :
#     else:
#         for handstar in starlist1:
#             closest = handstar.distanceandsort(transtarlist2)
#             if closest[0]['dist'] > tolerance:
#                 if verbose :
#                     print "No match for star %s" % handstar.name
#                 nomatch.append(handstar)
#                 continue
#
# Ok, then it must be closer then tolerance. We check for other stars whose distance is less then tolerance different from the first ones distance :
#             elif onlysingle and (closest[1]['dist'] - closest[0]['dist'] < tolerance):
#                 if verbose :
#                     print "Multiple candidates for star %s, skipping" % handstar.name
#                 notsure.append(handstar)
#                 continue
#
# Finally, this means we have found our star
#             else:
#                 idlist1.append(handstar)
#                 idlist2.append(starlist2[closest[0]['origpos']])
#                 iddists.append(closest[0]['dist'])
#
#         return returndict
#
#
#
# def listidentify(starlist1, starlist2, tolerance = 2.0, onlysingle = False, transform = False, scalingratio = 1.0, angle = 0.0, shift = (0.0, 0.0), verbose = True):
#     """
#     Same as formpairs (we call it), but we return only the idlist2 (not transformed, even if you give a transform), but with names taken from idlist1.
#     Typical : starlist2 is a sextractor catalog with random names, starlist 1 is a handpicked catalog with special names,
#     and you want to get stars with sextractor properties but your own names.
#     """
#
#     formpairsdict = formpairs(starlist1, starlist2, tolerance = tolerance, onlysingle = onlysingle, transform = transform, scalingratio = scalingratio, angle = angle, shift = shift, verbose = verbose)
#
#     match = []
#
#     for (s1, s2, d) in zip(formpairsdict["idlist1"], formpairsdict["idlist2"], formpairsdict["iddists"]):
#         s2.name = s1.name
#         s2.props["iddist"] = d
#         match.append(s2)
#
#     nomatchnames = [s.name for s in formpairsdict["nomatch"]]
#     notsurenames = [s.name for s in formpairsdict["notsure"]]
#
#     return {"match":match, "nomatchnames":nomatchnames, "notsurenames":notsurenames}
#
#
#
#
# def findtrans(preciserefmanstars, autostars, scalingratio = 1.0, tolerance = 2.0, minnbrstars = 5, mindist = 100.0, nref = 10, nauto = 30, verbose=True):
#
#     """
#     Finds a rotation and shift between two catalogs (a big dirty one and a small handpicked one).
#     Both catalogs should be SORTED IN FLUX, and the second one should be smaller for max performance.
#
#     Only the first nref stars of preciserefmanstars are considered for searching the possible matches, and furthermore only
#     pairs farther then mindist are considered.
#
#     tolerance is used when looking if a match was found.
#
#     minnbrstars = as soon as this number of stars are identified, the algo stops, we look no further.
#
#     The scalingratio parameter is a float to multiply with a distance of the autostars to match the same distance between the preciserefmanstars.
#
#     We return a dict of 3 things :
#     - nbr of identified stars (-1 if failed)
#     - rotation angle (center = 0,0)
#     - shift
#
#     This should then be used to transform your autostars, and then run listidentify between the catalogs if you want ...
#     This is done with the function formpairs
#
#     """
#
# Think of a Hubble expansion with "origin" (0,0)
# We apply this to the image to align, so that it matches the distances in the reference image.
#     autostarscopy = copy.deepcopy(autostars)
#     zoomstarlist(autostarscopy, scalingratio)
#
# n = 0 # a counter for the number of tries
# indentlist = [] # only used in case of failure
#
#     for b, brightstar in enumerate(preciserefmanstars[:nref]):
#         for f, faintstar in enumerate(preciserefmanstars[:nref]):
#             if f == b: continue
#             stardistance = brightstar.distance(faintstar)
#             if stardistance < mindist : continue
#
# We have a pair of stars from the preciserefmancat.
# Let's see if we find to stars in the autocat with a similar distance.
#
#             for bc, brightcandidate in enumerate(autostarscopy[:nauto]):
#                 for fc, faintcandidate in enumerate(autostarscopy[:nauto]):
#                     if fc == bc: continue
#                     candidatedistance =  brightcandidate.distance(faintcandidate)
#                     if math.fabs(candidatedistance - stardistance)/stardistance > 0.05 :
# So if there is a disagreement larger then 5 percent...
#                         continue
#
# We now have a promising pair of pairs, let's check them out.
#
#                     n = n+1
#
#                     starangle = brightstar.trigangle(faintstar)
#                     candidateangle = brightcandidate.trigangle(faintcandidate)
# rotcandangle = (starangle - candidateangle) % 360.0    # value to "add" to cand to match the star
#
# We apply this rotation to the bright candidate, to determine the shift :
#                     testcand = copy.deepcopy(brightcandidate)
#                     testcand.rotate(rotcandangle, (0, 0))
#                     candshift = testcand.findshift(brightstar)
#
# We apply the rotation and this shift to the full zoomed autostarlist :
#                     testcandlist = copy.deepcopy(autostarscopy)
#                     rotatestarlist(testcandlist, rotcandangle, (0, 0))
#                     shiftstarlist(testcandlist, candshift)
#
# We evaluate the match between the transformed autostars and the ref stars :
#
#                     pairsdict = formpairs(preciserefmanstars, testcandlist, tolerance = tolerance, onlysingle = True, verbose = False)
#                     nbrids = len(pairsdict["idlist1"])
#                     indentlist.append(nbrids)
#
#                     if nbrids >= minnbrstars:
# We got it !
#
#                         if verbose :
#                             print "Number of tries : %i" % n
#                             print "Distance difference : %.2f pixels" % math.fabs(candidatedistance - stardistance)
#                             print "Candidate rotation angle : %.2f degrees" % rotcandangle
#
#                             print "Star pairs used :"
#                             print brightstar
#                             print faintstar
#                             print brightcandidate
#                             print faintcandidate
#
#                             print "Identified stars : %i / %i" % (nbrids, len(preciserefmanstars) )
#
#                         return {"nbrids":nbrids, "angle":rotcandangle, "shift":candshift}
#
#
#     if verbose :
#         print "I'm a superhero, but I failed"
#     if len(indentlist) > 0:
#         if verbose :
#             print "Maximum identified stars : %i" % max(indentlist)
#
#     return {"nbrids":-1, "angle":0.0, "shift":(0.0, 0.0)}
#
#
#
# def writeforgeomap(filename, pairs):
#     """
#     Writes an input catalog of corresponding star pairs, for geomap
#     Pair is a list of couples like (refstar, startoalign)
#     """
#
#
#     import csv
#
#     table = []
#     for pair in pairs:
#         table.append([pair[0].x, pair[0].y, pair[1].x, pair[1].y])
#
# geomap = open(filename, "wb") # b needed for csv
#     writer = csv.writer(geomap, delimiter="\t")
#     writer.writerows(table)
#     geomap.close()
#
