
#-------------------------------------------------
# help
#-------------------------------------------------

"""

.. Label, use :ref:`text <label>` for reference
.. _utopya-module:

*****************
``utopya`` module
*****************

The ``utopya`` module is the top level access to the various 
:ref:`UTOPyA classes <utopya-classes>` and routines.
    
Actual implementations can be found in submodules:

.. The following are names of '.rst' files in the 'doc/sources' directory ;
   this ensures that each sub-module is given a seperate page in the documentation.
    
.. toctree::
   :maxdepth: 1

   python-module-utopya_base
   python-module-utopya_rc
   python-module-utopya_jobscript
   python-module-utopya_jobtree
   python-module-utopya_runscript
   python-module-utopya_install
   python-module-utopya_tools


Other modules
=============
    
Tools are available from other modules.

.. The following are names of '.rst' files in the 'doc/sources' directory ;
   this ensures that each sub-module is given a seperate page in the documentation.
    
.. toctree::
   :maxdepth: 1

   python-module-rc
   python-module-gss


.. Label, use :ref:`text <label>` for reference
.. _utopya-classes:

Utopya classes
==============

The UTOPyA python classes could be imported from the 
:ref:`utopya module <utopya-module>`,
and are defined according to the following hierchy:

* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

    * :py:class:`.UtopyaJobStep`

      * :py:class:`.UtopyaJobTree`

        * :py:class:`.UtopyaJobIteration`
      
    * :py:class:`.UtopyaInstall`

    
  * :py:class:`.UtopyaJobTask`

    * :py:class:`.UtopyaJobTaskSubmit`

  * :py:class:`.UtopyaRunScript`

  * :py:class:`.UtopyaJobScript`

    * :py:class:`.UtopyaJobScriptForeground`
    * :py:class:`.UtopyaJobScriptRedirect`
    * :py:class:`.UtopyaJobScriptBackground`
    * :py:class:`.UtopyaJobScriptBatch`

      * :py:class:`.UtopyaJobScriptBatchTest`
      * :py:class:`.UtopyaJobScriptBatchLoadLeveler`
      * :py:class:`.UtopyaJobScriptBatchLSF`

"""

#-------------------------------------------------
# imports from sub modules
#-------------------------------------------------

from utopya_base       import *
from utopya_rc         import *
from utopya_jobscript  import *
from utopya_jobtree    import *
from utopya_runscript  import *
from utopya_install    import *
from utopya_tools      import *


#-------------------------------------------------
# end
#-------------------------------------------------
