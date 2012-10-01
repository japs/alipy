Installation
============


Requirements
------------

We need python, numpy, scipy, matplotlib, and :
 * `Source Extractor <http://www.astromatic.net/software/sextractor>`_ by Bertin & Arnouts.
 * `astroasciidata <http://www.stecf.org/software/PYTHONtools/astroasciidata/>`_
 * `f2n <http://obswww.unige.ch/~tewes/f2n_dot_py/>`_ (optional, but very useful to make control visualizations)

We use `pysex <http://pypi.python.org/pypi/pysex/>`_ to interact with SExtractor, but pysex comes bundled with alipy -- no need to install.


Download
--------

Get alipy directly from its public repository, by typing::

	svn checkout https://svn.epfl.ch/svn/mtewes-public/trunk/alipy2 ./alipy

And then,
::

	cd alipy
	python setup.py install

	
More info about the installation : `<http://docs.python.org/distutils/>`_.


To generate this documentation
------------------------------

You'll need `sphinx <http://sphinx.pocoo.org/>`_.

::
	
	cd doc
	make apidoc
	make html

And then open _build/html/index.html