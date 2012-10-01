About
=====

This is a python package to quickly, automatically, and robustly identify geometrical transforms between optical astronomical images, using only field stars. The images can have different pixel sizes, orientations, pointings and filters.

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