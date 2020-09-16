.. Documentation description.

.. Label between '.. _' and ':' ; use :ref:`text <label>` for reference
.. _documentation:


*************
Documentation
*************


Generation
==========

The documentation of EMIP that you are reading right now is generated
from EMIP source files.
To (re)generate it, run in the main directory::

  make docu

The documentation files are formatted for use with 'Sphinx'; see:

* `Sphinx homepage <http://sphinx-doc.org/index.html>`_

The 'Sphink' python package should be available to be able to (re)generate the documentation.
If not available on your system yet, see below for installation instructions.


Source files
============

The main files of the documentation are located in::

  ./doc/source/

Parts of the documentation are included in the various Python modules and scripts.
The documentation is configured in such a way that these are automatically incorporated.

The source tree of the documentation files is:

* `index.rst <../../source/index.rst>`_
  
  * `documentation.rst <../../source/documentation.rst>`_

Installation of Sphinx
======================

To check if Sphinx is already installed, try to run the quick-start script::

  sphinx-quickstart --version
  
  sphinx-quickstart 1.8.1

If not available yet, try to install using one of the following methods.

Install Sphinx
--------------

Try if it is possible to download and install Sphinx using 
one of the standard installation commands.

* For an `Anaconda <https://www.anaconda.com/>`_ distribution, try::

    conda install -c anaconda sphinx
    
* For other distributions, try::

    pip install sphinx

Build Spinx from source 
-----------------------

Download the latest version from the
`Sphinx homepage <http://sphinx-doc.org/index.html>`_,
for example::

  Sphinx-1.2.tar.gz

To let the package be installed at a local destination, 
define the following environment variable::

  export PYTHON_PREFIX="${HOME}/opt/Python-2.7.6"

Build and install in the user-defined location using::

  cd Sphinx-1.2
  python setup.py install --user

Extend the search path with::

  export PATH="${PYTHON_PREFIX}/bin:${PATH}"


Configuration
=============

The source of the documentation is located in (subdirectories of):

  ./doc/

In this directory, the Sphinx quick start command has been called::

  sphinx-quickstart
  
The following settings were entered 
(lines marked with ``<--`` are non-default)::

  Welcome to the Sphinx 1.2 quickstart utility.
  > Root path for the documentation [.]:
  > Separate source and build directories (y/N) [n]: y                          <--
  > Name prefix for templates and static dir [_]: 
  > Project name: EMIP                                                          <--
  > Author name(s): Arjo Segers                                                 <--
  > Project version: trunk                                                      <--
  > Project release [trunk]:     
  > Source file suffix [.rst]: 
  > Name of your master document (without suffix) [index]: 
  > Do you want to use the epub builder (y/N) [n]: 
  > autodoc: automatically insert docstrings from modules (y/N) [n]: y          <--
  > doctest: automatically test code snippets in doctest blocks (y/N) [n]: 
  > intersphinx: link between Sphinx documentation of different projects (y/N) [n]: 
  > todo: write "todo" entries that can be shown or hidden on build (y/N) [n]: 
  > coverage: checks for documentation coverage (y/N) [n]: 
  > pngmath: include math, rendered as PNG images (y/N) [n]: 
  > mathjax: include math, rendered in the browser by MathJax (y/N) [n]: y      <--
  > ifconfig: conditional inclusion of content based on config values (y/N) [n]: 
  > viewcode: include links to the source code of documented Python objects (y/N) [n]: 
  > Create Makefile? (Y/n) [y]: 
  > Create Windows command file? (Y/n) [y]: n                                   <--

The following entities have been created now:

* ``source`` directory to hold the documenation source files;
  initially the following files are created:
  
  * `conf.py <../../source/conf.py>`_ : configuration file;
  * `index.rst <../../source/index.rst>`_ : source of the main page of the documentation;

* ``build`` directory to hold the generated documenation;

* ``Makefile`` : make commands.

In the ``./doc/source/conf.py`` file created in this way,
the following changes were made manually:

* The location of the python modules was added to the search path::

    # If extensions (or modules to document with autodoc) are in another directory,
    # add these directories to sys.path here. If the directory is relative to the
    # documentation root, use os.path.abspath to make it absolute, like shown here.
    sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,os.pardir,'py')) )

* Added default options for `Autodoc <http://sphinx-doc.org/ext/autodoc.html>`_ entries, 
  these are used to create documentation out of Python modules::

    autodoc_default_flags = ['show-inheritance','members']
    autodoc_member_order  = 'bysource'
    intersphinx_mapping   = { 'python' : ('http://docs.python.org/2.7',None) }

* The html theme was changed for a less default layout;
  this requires that the '``html_sidebars``' section is commented::

    html_theme = 'nature'

The initial documentation could be created using:

    (cd doc; make html)


How to add a new module to the documentation?
=============================================

Example taken from introduction of the 'emip_tropomi' module.

* Create new module file based on template::

    cp  py/empi_omi.py  py/emip_tropomi.py
      
  Replace names from the template ("OMI") by the new names ("TROPOMI").
   
* Create a documentation file specific for the new module::

    doc/source/python-module-emip_tropomi.rst

  which only refers to the documentation included in the module::

    .. automodule:: emip_tropomi

* Add a reference to the new module in the :ref:`modules-and-classes` page
  in the file `python-modules.rst <../../source/python-modules.rst>`_::

    .. toctree::
       :maxdepth: 1

       python-module-emip_omi
       python-module-emip_tropomi              <---
       :


* Add a reference to the module documentation to the top-level 
  table-of-content in `index.rst <../../source/index.rst>`_::

    .. toctree::
       :maxdepth: 2

       tutorial
       emip_omi
       emip_tropomi                            <---
       python-modules
       documentation

   
