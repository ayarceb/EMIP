EMIP Documentation
==================

Using the Sphinx system for python documentation.

To build the documenation from here, use:

  make html


Installation of Sphinx
----------------------

The 'Sphink' python package should be installed.
Try if 'sphinkx-quickstart' is available.

For installation in user directory.
Environment variable defines destination:
  export PYTHON_PREFIX="${HOME}/opt/Python-2.7.4"
Download and unpack:
  Sphinx-1.2.tar.gz
Build:
  cd Sphinx-1.2
  python setup.py install --user
Extend search path:
  export PATH="${PYTHON_PREFIX}/bin:${PATH}"


Configuration of documentation
------------------------------

To create a new documentation, use the Sphinx Quick Start tool:

  spinkx-quickstart
  
The following settings were entered 
(lines marked with '<--' are non-default):

  Welcome to the Sphinx 1.2 quickstart utility.
  > Root path for the documentation [.]:
  > Separate source and build directories (y/N) [n]: y                          <--
  > Name prefix for templates and static dir [_]: 
  > Project name: EMIP                                                          <--
  > Author name(s): Arjo Segers                                                 <--
  > Project version: trunk                                                      <--
  > Project release [trunk]:
  > Project language [en]: 
  > Source file suffix [.rst]: 
  > Name of your master document (without suffix) [index]: 
  > Do you want to use the epub builder (y/N) [n]: 
  > autodoc: automatically insert docstrings from modules (y/N) [n]: y          <--
  > doctest: automatically test code snippets in doctest blocks (y/N) [n]: 
  > intersphinx: link between Sphinx documentation of different projects (y/N) [n]: y  <--
  > todo: write "todo" entries that can be shown or hidden on build (y/N) [n]: 
  > coverage: checks for documentation coverage (y/N) [n]: 
  > pngmath: include math, rendered as PNG images (y/N) [n]: 
  > mathjax: include math, rendered in the browser by MathJax (y/N) [n]: y      <--
  > ifconfig: conditional inclusion of content based on config values (y/N) [n]: 
  > viewcode: include links to the source code of documented Python objects (y/N) [n]: 
  > Create Makefile? (Y/n) [y]: 
  > Create Windows command file? (Y/n) [y]: n                                   <--

The following entities have been created now:

* 'source' directory to hold the documenation source files;
  initially the following files are created:
  
  * conf.py    : configuration file;
  * index.rst  : source of the main page of the documentation;

* 'build' directory to hold the generated documenation;

* 'Makefile' : make commands.

In the 'source/conf.py' file created in this way,
the following changes were made manually:

* The location of the python modules was added to the search path:

    # If extensions (or modules to document with autodoc) are in another directory,
    # add these directories to sys.path here. If the directory is relative to the
    # documentation root, use os.path.abspath to make it absolute, like shown here.
    sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,os.pardir,'py')) )

* Added default options for 'Autodoc' entries, 
  these are used to create documentation out of Python modules:

    autodoc_default_flags = ['show-inheritance','members']
    autodoc_member_order  = 'bysource'
    intersphinx_mapping   = { 'python' : ('http://docs.python.org/2.7',None) }

* The html theme was changed for a less default layout;
  this requires that the 'html_sidebars' section is commented:

    html_theme = 'bizstyle'

The initial documentation could be created using:

    make html

