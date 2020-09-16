
#-------------------------------------------------
# help
#-------------------------------------------------

"""
***************
``emip`` module
***************

The :py:mod:`emip` module provides classes
for preprocessing EMEP input data.


Class hierchy
=============

The classes and are defined according to the following hierchy:

* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

    * :py:class:`.Emip`

    * :py:class:`.EmipTask`

      * :py:class:`.EmipTaskList`

      * :py:class:`.EmipMessage`


Classes
=======


"""


#-------------------------------------------------
# modules
#-------------------------------------------------

# tools:
import utopya_rc


#-------------------------------------------------
# main
#-------------------------------------------------

class Emip( utopya_rc.UtopyaRc ) :

    """
    Main class to process EMIP task list.
    
    Initialization by arguments:
     
    * ``rcfile``  : name of settings file
     
    Optional arguments:
    
    * ``rcbase``, ``env`` : arguments passed to rc base class
        
    """
    
    def __init__( self, rcfile, rcbase='', env={} ) :
    
        """
        Initialize EMIP object.
        """
                        
        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # info ...
        self.logger.info( 'EMIP main task' )
        
        # load class:
        cls = self.ImportClass( 'class' )
        # initialize:
        self.main = cls( rcfile, rcbase=rcbase, env=env )        
        
    #enddef __init__

#endclass Emip


#-------------------------------------------------
# main
#-------------------------------------------------

class EmipTask( utopya_rc.UtopyaRc ) :

    """
    Base class for EMIP tasks.
    The name of an rcfile with configuration settings should be passed as argument.
    """
    
    def __init__( self, rcfile, rcbase='', env={} ) :
    
        """
        Initialize EMIP task object.
        """
                        
        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
    #enddef __init__

#endclass EmipTask


#-------------------------------------------------
# main
#-------------------------------------------------

class EmipTaskList( EmipTask ) :

    """
    Class to define task list.
    
    The name of the sub-tasks to be performed is read from the setting::
    
      <rcbase>.tasks     :   download convert archive
      
    For each of these sub-tasks, the settings should define the class
    and initialization arguments, for example::
    
      <rcbase>.download.class    :  mymodule.Download
      <rcbase>.download.args     :  '%{rcfile}'
        
    The sub-task could be a list too, for example::

      <rcbase>.convert.class    :  emip.EmipTaskList
      <rcbase>.convert.args     :  '%{rcfile}'
      <rcbase>.convert.tasks    :  csv2netcdf regrid
        
    """
    
    def __init__( self, rcfile, rcbase='', env={}, indent='' ) :
    
        """
        Initialize EMIP task list object.
        """
                        
        # init base object:
        EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # info ...
        self.logger.info( indent+'run task list ...' )
        
        # read list:
        self.tasks = self.GetSetting( 'tasks' ).split()
        # loop:
        for task in self.tasks :
            # info ...
            self.logger.info( indent+'  task "%s" ...' % task )
        
            # info ...
            self.logger.info( indent+'    load class ...' )
            # load class:
            cls = self.ImportClass( task+'.class' )

            # info ...
            self.logger.info( indent+'    arguments ...' )
            # arguments as string:
            args = self.GetSetting( task+'.args' )
            # replace templates:
            args = args.replace( '%{rcfile}', rcfile )
            
            # add indent argument for EMIP classes:
            if 'Emip' in cls.__name__ :
                args = args+", indent=indent+'    '"
            #endif

            # info ...
            self.logger.info( indent+'    evaluate ...' )
            # initialize:
            obj = eval( "cls( %s, rcbase='%s.%s', env=env )" % (args,rcbase,task) ) 

        #endfor
        
    #enddef __init__

#endclass EmipTaskList


#-------------------------------------------------
# message only
#-------------------------------------------------

class EmipMessage( EmipTask ) :

    """
    EMIP task that only displays a message;
    useful for testing.
    
    The message to be displayed is read from the rcfile::
    
       [<rcbase>].msg       :  Sorry, task '%{name}' not implemented yet ...
       
    The following templates might be used:
    
    * '``%{name}``' inserts the name of the task (the rcfile base)

    """
    
    def __init__( self, rcfile, rcbase='', env={}, indent='' ) :
    
        """
        Display message.
        """
                        
        # init base object:
        EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # read message:
        msg = self.GetSetting( 'msg' )
        # replace:
        msg = msg.replace( '%{name}', self.rcbase )
        # info ...
        self.logger.info( indent+msg )
        
    #enddef __init__

#endclass EmipMessage


#-------------------------------------------------
# end
#-------------------------------------------------



