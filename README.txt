About
=====

This is a python package to quickly and automatically identify geometrical transforms between conventional optical astronomical images of a given patch of sky, using only field stars. The images can have different pixel sizes, orientations, pointings and filters.

http://obswww.unige.ch/~tewes/


Installation
============

Quick :
python setup.py install

To create a source distribution :
python setup.py sdist

More info :
http://docs.python.org/distutils/


To generate the documentation
=============================

cd doc
make apidoc
make html

--> _build/index.html