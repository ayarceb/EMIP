
.. Documentation description.

.. Label between '.. _' and ':' ; use :ref:`text <label>` for reference
.. _tropomi-processing:

***********************
TROPOMI data processing
***********************

This chapter describes the tasks performed for processing OMI data.


Configuration
=============

An example configuration could be found in:
  
* `rc/emip.rc <../../../rc/emip_omi.rc>`_

The examples bellow assume that a task '``emip.omi``' has been
configured as a list of sub tasks::

  ! initialize:
  emip.tropomi.class        :  emip.EmipTaskList
  emip.tropomi.args         :  '%{rcfile}'

  ! task list:
  emip.tropomi.tasks        :  download convert catalogue



Downloading original data files
===============================

OMI data is available from the `TEMIS <http://temis.nl/>`_ website.
Files are collected per day in files named::

  http://www.temis.nl/airpollution/no2col/data/tropomi/2018/07/tropomi_no2_20180728.tar

The EMIP task '``emip.tropomi.download``' could be used to download and unpack
a series of these files.
The work is done by by the :py:class:`.UtopyaInstall` class from the :py:mod:`utopya` module;
the task initalized with::

  emip.tropomi.download.class        :  utopya.UtopyaInstall
  emip.tropomi.download.args         :  '%{rcfile}'

The class requires rcfile settings that specify the files to be downloaded and
where to store them.
The first setting is a list of data sets nicknames that should be download;
here just a single set of 'daily' files is requested::

  emip.tropomi.download.sets                :  daily

For each data the details should be specified.
One file per day should be installed; specify the time range and time step with::

  emip.tropomi.download.daily.timerange   :  2012-01-01 upto 2012-01-02 by 1 day

Specify the archive location as an url; this path should not include templates for time values::

  emip.tropomi.download.daily.arch        :  http://temis.nl/airpollution/no2col/data/tropomi

The filename(s) to be download could include templates for the time
which are expanded by the :py:meth:`datetime.datetime.strftime` method.
The extension '``.tar``' will enforce automatic unpacking::

  emip.tropomi.download.daily.files       :  %Y/%m/tropomi_%Y%m%d.tar

Specify the target location with::

  emip.tropomi.download.daily.dir         :  /scratch/TEMIS/airpollution/no2col/data/tropomi

Eventually enable a flag to show debug messages::

  emip.tropomi.download.daily.verbose     :  True


Convert and select
==================

The '``emip.tropomi.convert``' task creates netCDF files with selected pixels,
for example only those within some region or cloud free pixels.
The selection criteria are defined in the settings, and added
to the '``history``' attribute of the created files to not be forgotten.

The work is done by the :py:class:`.EmipTropomiConvert` class,
which is initialized using the settings file::

  ! task initialization:
  emip.tropomi.convert.class                  :  emip_omi.EmipTropomiConvert
  emip.tropomi.convert.args                   :  '%{rcfile}'
  
See the class documentation for the configuration.



Pixel selection
---------------
        
The :py:meth:`.TROPOMI_NC4_File.SelectPixels` method is called
to create a pixel selection mask.
The selection is done specification of one or more filters.
First provide a list of filter names::

  emip.tropomi.convert.filters   :  lons lats valid tcflag albedo

Then provide for each filter the details of the input variable,
the assumed units (safety check!), the type of filter,
and (depending on the type) other settings that specify the allowed values.
For example, the following ``lons`` filter selects a longitude range::

  emip.tropomi.convert.filter.lons.var                :  Geolocation Fields/Longitude
  emip.tropomi.convert.filter.lons.units              :  deg
  emip.tropomi.convert.filter.lons.type               :  minmax
  emip.tropomi.convert.filter.lons.minmax             :  -30.0 45.0
  
See the the description of the method for the supported filter types
and their required settings.


Variable specification
----------------------

The converted data is stored in an :py:class:`.OmiExtract` object
which will take care of writing to a netCDF file.
The :py:meth:`.OmiExtract.AddSelection` method is used to process
the selected pixels, e.g. selecting the variables, apply conversions, etc.

The first setting that is read is a list with variable names to be 
created in the target file::

  emip.tropomi.convert.output.vars    :  longitude corner_longitudes \
                                         latitude corner_latitudes \
                                         vcd_trop  ...

For each variable settings should be specified that describe 
how to obtain the values and the target units (for automatic conversion if possible).

For most variables it is sufficient to provide only the name of the original
variable from which the data should be read::

  emip.tropomi.convert.output.var.longitude.from    :   PRODUCT/longitude
  emip.tropomi.convert.output.var.longitude.units   :   degrees_east

For some variables some special processing needs to be done.
For these variables a key '``special``' is used which will enable the 
correct code for conversion. For example, the following setting will
ensure that the variable '``image_number``' will be filled with the scan
number within the track::

  emip.tropomi.convert.output.var.image_number.special            :   scan_number
  emip.tropomi.convert.output.var.image_number.units              :   1

The special conversions are implemented in the :py:meth:`.TropomiExtract.AddSelection` method.
If new variables require special processing, just insert a new '``special``' keyword
and wait for the method to complain about an unsupported value.


Output files
------------

The converted data is written to a file specified by directory
and filename templates.

The output directory could include templates for time values::

  emip.tropomi.convert.output.dir    :  /data/OMI-selection/%Y/%m

Also the filenames could include time values,
and in addition a template '``%{orbit}``' in which the orbit id
from the OMI file is inserted::

  emip.tropomi.convert.output.filename    :  OMI-Aura_NO2_%Y%m%d_%{orbit}.nc

Specify global attributes and their value with::

  emip.tropomi.convert.output.attrs  :  format Conventions author institution email

  emip.tropomi.convert.output.attr.format         :  1.0
  emip.tropomi.convert.output.attr.Conventions    :  CF-1.6
  emip.tropomi.convert.output.attr.author         :  Arjo Segers
  emip.tropomi.convert.output.attr.institution    :  MetNorway, Oslo, Norway
  emip.tropomi.convert.output.attr.email          :  Arjo.Segers@met.no




Regridding
==========

The '``emip.tropomi.regrid``' task could be used to resample the
pixels onto a regular grid.

The work is done by the :py:class:`.EmipOmiRegrid` class,
which is initialized using the settings file::

  emip.tropomi.regrid.class            :  emip_omi.EmipOmiRegrid
  emip.tropomi.regrid.args             :  '%{rcfile}'

Input files are expected to be produced by the
:py:class:`.EmipOmiConvert` class.
For each input file, an output file with similar format is created
which has 'pixels' with a footprint equal to a grid cell.
The regular grid (centers, corners) is saved as the 'track'.

The remapping is done by distributing the footprint polygon over the grid cells.
Each grid cell is therefore filled with a weighed sum of contributions from pixels
that (partly) overlap the cell:

.. math::
    y ~=~ \sum\limits_{i=1}^{npix} w_i x_i

The weights are relative to the area covered by a pixel, thus the more area of a cell
is covered by a pixel the more weight that pixel has.

Variables that start with '``sigma_``' are assumed to be error estimates.
The errors are assumed to be uncorrrelated between pixels, and therefore
the combined error could be computed as weighted sum over variances:

.. math::
    \sigma_y ~=~ \sqrt{ \sum\limits_{i=1}^{npix} w_i\ \sigma_{x,i}^2 }

See the class documentation for detailed configuration settings.



Catalogues
==========

The '``emip.tropomi.catalogue``' task could be used to create a catalogue
of figures extracted from the converted files.

.. figure:: figs/omi-catalogue.png
   :scale: 50 %
   :align: center
   :alt: OMI image catalogue
   
   *Example of image catalogue produced from converted OMI files.*

The work is done by the :py:class:`.EmipOmiCatalogue` class,
which is initialized using the settings file::

  emip.tropomi.catalogue.class            :  emip_omi.EmipOmiCatalogue
  emip.tropomi.catalogue.args             :  '%{rcfile}'

In the settings, specify a time range for which images should be created::

  emip.tropomi.catalogue.timerange.start  :  2012-01-01 00:00
  emip.tropomi.catalogue.timerange.end    :  2012-12-31 23:59

The images as well as an html index are written to a single directory
specified with::

  emip.tropomi.catalogue.dir    :  /data/OMI-selection/catalgoue

The location of the converted OMI files could be specified with
time templates and a filename filter::

  ! converted OMI files, absolute path or relative to catalogue:
  emip.tropomi.catalogue.input.filenames        :  ../%Y/%m/OMI-Aura_NO2_*.nc

Specify a list of variables to be plotted; usually this is the tropospheric
vertical column density that is the main product, but also variables might be
of interest::

  emip.tropomi.catalogue.vars                   :  vcd_trop

Per variable the maximum value for the color bar could be specified;
if not defined, the color bar is simply stretched to the maximum value
present in the data:

  emip.tropomi.catalogue.var.vcd_trop.vmax      :  20.0

Specify the domain of the map, projection is regular longitude/latitude:

  ! map domain (west east south north):
  emip.tropomi.catalogue.domain       :  -30 45 35 75

Enable the following flag to re-create existing files,
by default only non-existing files are created::

  emip.tropomi.catalogue.renew                  :  False

When finished, the :py:mod:`catalogue` module is used create index pages.
The url that should be loaded in a browser is shown::

  Point your browser to :
    file:///data/OMI-selection/catalgoue/index.html

