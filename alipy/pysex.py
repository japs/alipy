"""
Author: Nicolas Cantale - n.cantale@gmail.com

Small module wrapping sextractor. The idea is to have a single function taking an image and returning a sextractor catalog.

Dependencies:
 - sextractor (mandatory)
 - astroasciidata (mandatory)
 - numpy (optional, needed for the array support)
 - pyfits (optional, needed for the array support)


Usage:

	import pysex
	cat = pysex.run(myimage, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER'], conf_args={'PHOT_APERTURES':5})
	print cat['FLUX_APER']

"""

import os
import shutil
import asciidata


def _check_files(conf_file, conf_args, verbose=True):
    if conf_file is None:
        os.system("sex -d > .pysex.sex")
        conf_file = '.pysex.sex'
    if not conf_args.has_key('FILTER_NAME') or not os.path.isfile(conf_args['FILTER_NAME']):
        if verbose:
            print 'No filter file found, using default filter'
        f = open('.pysex.conv', 'w')
        print>>f, """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1"""
        f.close()
        conf_args['FILTER_NAME'] = '.pysex.conv'
    if not conf_args.has_key('STARNNW_NAME') or not os.path.isfile(conf_args['STARNNW_NAME']):
        if verbose:
            print 'No NNW file found, using default NNW config'
        f = open('.pysex.nnw', 'w')
        print>>f, """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	 9 for profile parameters + 1 for seeing.
# outputs:	  ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02	8.31957e-01	 2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00	3.34497e+00	 3.51885e+00  3.65570e+00  3.74856e+00	3.84541e+00	 4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01	 5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01	 1.29492e-01  1.42290e+00  2.90741e+00	2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01	2.18845e+00	 2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02	 9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00	5.59549e-01	 8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00	9.87912e-01	 6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00	1.45967e+00	 1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01	 1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00	 8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01	 1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01	1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01	1.27910e-01	 1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00	2.42038e+00	 1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00
 1.00000e+00"""
        f.close()
        conf_args['STARNNW_NAME'] = '.pysex.nnw'

    return conf_file, conf_args


def _setup(conf_file, params):
    try:
        shutil.copy(conf_file, '.pysex.sex')
    except:
        pass  # already created in _check_files

    f = open('.pysex.param', 'w')
    print>>f, '\n'.join(params)
    f.close()


def _setup_img(image, name):
    if not type(image) == type(''):
        import pyfits
        pyfits.writeto(name, image)


def _get_cmd(img, img_ref, conf_args):
    ref = img_ref if img_ref is not None else ''
    cmd = ' '.join(['sex', ref, img, '-c .pysex.sex '])
    args = [''.join(['-', key, ' ', str(conf_args[key])]) for key in conf_args]
    cmd += ' '.join(args)
    return cmd


def _read_cat(path='.pysex.cat'):
    cat = asciidata.open(path)
    return cat


def _cleanup():
    files = [f for f in os.listdir('.') if '.pysex.' in f]
    for f in files:
        os.remove(f)

# def run(image='', imageref='', params=[], conf_file=DEFAULT_CONF,
# conf_args={}):


def run(image='', imageref='', params=[], conf_file=None, conf_args={}, keepcat=True, rerun=False, catdir=None):
    """
    Run sextractor on the given image with the given parameters.

    image: filename or numpy array
    imageref: optional, filename or numpy array of the the reference image
    params: list of catalog's parameters to be returned
    conf_file: optional, filename of the sextractor catalog to be used
    conf_args: optional, list of arguments to be passed to sextractor (overrides the parameters in the conf file)


    keepcat : should I keep the sex cats ?
    rerun : should I rerun sex even when a cat is already present ?
    catdir : where to put the cats (default : next to the images)


    Returns an asciidata catalog containing the sextractor output

    Usage exemple:
            import pysex
            cat = pysex.run(myimage, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER'], conf_args={'PHOT_APERTURES':5})
            print cat['FLUX_APER']
    """

    # Preparing permanent catalog filepath :
    (imgdir, filename) = os.path.split(image)
    (common, ext) = os.path.splitext(filename)
    catfilename = common + \
        ".pysexcat"  # Does not get deleted by _cleanup(), even if in working dir !
    if keepcat:
        if catdir:
            if not os.path.isdir(catdir):
                os.makedirs(catdir)
                # raise RuntimeError("Directory \"%s\" for pysex cats does not
                # exist. Make it !" % (catdir))

    if catdir:
        catpath = os.path.join(catdir, catfilename)
    else:
        catpath = os.path.join(imgdir, catfilename)

    # Checking if permanent catalog already exists :
    if rerun == False and type(image) == type(''):
        if os.path.exists(catpath):
            cat = _read_cat(catpath)
            return cat

    # Otherwise we run sex :
    conf_args['CATALOG_NAME'] = '.pysex.cat'
    conf_args['PARAMETERS_NAME'] = '.pysex.param'
    if 'VERBOSE_TYPE' in conf_args and conf_args['VERBOSE_TYPE'] == 'QUIET':
        verbose = False
    else:
        verbose = True
    _cleanup()
    if not type(image) == type(''):
        import pyfits
        im_name = '.pysex.fits'
        pyfits.writeto(im_name, image.transpose())
    else:
        im_name = image
    if not type(imageref) == type(''):
        import pyfits
        imref_name = '.pysex.ref.fits'
        pyfits.writeto(imref_name, imageref.transpose())
    else:
        imref_name = imageref
    conf_file, conf_args = _check_files(conf_file, conf_args, verbose)
    _setup(conf_file, params)
    cmd = _get_cmd(im_name, imref_name, conf_args)
    res = os.system(cmd)
    if res:
        print "Error during sextractor execution!"
        _cleanup()
        return

    # Keeping the cat at a permanent location :
    if keepcat and type(image) == type(''):
        shutil.copy('.pysex.cat', catpath)

    # Returning the cat :
    cat = _read_cat()
    _cleanup()
    return cat
