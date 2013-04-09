Notes
=====


Multi-hdu support
-----------------

When SExtractor runs on a multi-extension ("multi-hdu") image, it puts the sources from all extension into the same catalog.
The param EXT_NUMBER in the catalog allows to separate these sources afterwards.
Doing this when using alipy is crucial, otherwise very wrong sources (bad pixels !) pollute the catalogs and completely prevent alipy from finding correct identifications.

Therefore, when using images with several extensions, you will have to specify, e.g., ``hdu = 1``.
Alipy will then only consider sources extracted with this particular EXT_NUMBER.
The default setting is ``hdu = 0``, and in this case alipy will use the `only` EXT_NUMBER which it finds in the catalog (not necessarily 0 !), and exit if it finds more then one EXT_NUMBER.

