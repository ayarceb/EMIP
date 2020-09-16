
********
Tutorial
********

This chapter describes step by start how to run EMIP.


Run script
==========

Use the main run script to start the processing.
Withouot arguments, the default configuration is used::

  ./bin/emip

This is actually the same as passing the following arguments::

  ./bin/emip  rc/emip.rc  --rcbase='emip'

The first argument '`rc/emip.rc <../../../rc/emip.rc>`_' specifies a configuration file.
The content is formatted similar to an `X-resource` file,
and therefore has extension ``.rc``.
For a description of the format, see the section on :ref:`rcfile formatting <rc-formatting>`
in the :py:mod:`rc` module.

The keywords in the rcfile all start with '``emip``', 
as defined by the optional '``rcbase``' argument.
The user could pass a different rcfile or different base if necessary.


Task tree
=========

The configuration defines a series of tasks to be perfomed.
For example, the following tasks might be defined to process OMI satellite data::

  emip.omi.download
  emip.omi.convert
  emip.omi.regrid

This list is actually defined as tree, using lists in which the elements could be lists too::

  emip                 # list with elements "omi"
      .omi             # list with elements "download", "convert", and "regrid"
          .download    # first task
          .convert     # second task
          .regrid      # etc

For each element in the tree, the configuration file should specify a python class name
and the arguments that should be used to initalize it.
These classes could be available by default in the EMIP modules already,
but could also be user defined.
The first element in the task tree is the trunk '``emip``', 
which is configured to be a list of tasks using standard EMIP class::

  emip.class        :  emip.EmipTaskList
  emip.args         :  '%{rcfile}'

The :py:class:`.EmipTaskList` class that should be used is preceeded by the :py:mod:`emip` module
name in which it is implemented. Eventually a path to the module could be added in case that
it is not on the python search path::

  emip.class        :  /path/to/emip.EmipTaskList

The only required arguments to initialize the :py:class:`.EmipTaskList` is the name of a
configuration file; if this is the same file as the one that contains the class/args definition,
then the '``%{rcfile}``' template could be used.

The :py:class:`.EmipTaskList` class will try to read a list of task names from the configuration.
The key should start with the name of the list ('``emip``'), 
and defines in this example only a single element::

  emip.tasks        :  omi

The '``omi``' task is configured using settings with the full path in the task tree,
thus '``emip.omi``'. In this example, the settings should define a task list again::

  emip.omi.class        :  emip.EmipTaskList
  emip.omi.args         :  '%{rcfile}'
  emip.omi.tasks        :  download convert regrid

The 3 sub-tasks defined here are configured using classes that actually do some real work::

  emip.omi.download.class        :  emip_omi.Download
  emip.omi.download.args         :  '%{rcfile}'

  emip.omi.convert.class         :  emip_omi.Convert
  emip.omi.convert.args          :  '%{rcfile}'

  emip.omi.regrid.class          :  emip_omi.Regrid
  emip.omi.regrid.args           :  '%{rcfile}'

In this example, the classes that do the work are implemented in the :py:mod:`emip_omi` module.
The work to be done is defined by settings read from the configuration file.

See the chapter on :ref:`omi-processing` for details of the configuration.


Documentation
=============

A '``Makefile``' is present to (re)create the documentation::

  make docu
  
To remove the created documentation and other temporary files, use::

  make clean


