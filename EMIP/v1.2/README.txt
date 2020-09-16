EMIP - EMEP Input Processor
===========================

Created for conversion of satellite data for use in EMEP-DAS.


Quick start
-----------

Start the main script with a settings file:

  ./bin/emip  rc/emip.rc


Documenation
------------

Generated from source files, run:

  make docu

and browse to:

  doc/build/html/index.html
  

Changes
-------

- Download using 'wget' first to temporary file.
- Support conversion of TROPOMI NO2 files that are still in tar file.
- Extract TROPOMI NO2 kernel data.
Frozen into v1.2.





