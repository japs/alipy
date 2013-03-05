Installation
============


Requirements
------------

We need python, numpy, scipy, matplotlib, and :
 * `Source Extractor <http://www.astromatic.net/software/sextractor>`_ by Bertin & Arnouts. The executable must be available as ``sex``, not as ``sextractor``. Make an alias if required.
 * `astroasciidata <http://www.stecf.org/software/PYTHONtools/astroasciidata/>`_ : a package to read SExtractor catalogs
 * `f2n <http://obswww.unige.ch/~tewes/f2n_dot_py/>`_ : **optional**, but very useful to make control visualizations
 * `PyRAF/IRAF <http://www.stsci.edu/institute/software_hardware/pyraf>`_ : **optional**, only needed if you want to use IRAF's geomap/gregister image alignment

We use `pysex <http://pypi.python.org/pypi/pysex/>`_ by Nicolas Cantale to interact with SExtractor, but pysex comes bundled with alipy -- no need to install it.


Download
--------

Get alipy directly from its public repository, by typing::

	svn checkout https://svn.epfl.ch/svn/mtewes-public/trunk/alipy2 ./alipy

And then,
::

	cd alipy
	python setup.py install

or maybe
::

	python setup.py install --user

... if you don't have write access to the global site-packages directory of your machine.
More info about the installation : `<http://docs.python.org/install/>`_.


To generate this documentation
------------------------------

You'll need `sphinx <http://sphinx.pocoo.org/>`_.

::
	
	cd doc
	make apidoc
	make html

And then open _build/html/index.html