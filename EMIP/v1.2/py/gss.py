#! /usr/bin/env python

"""
.. Label, use :ref:`text <label>` for reference
.. _gss-module:

**************
``gss`` module
**************

The gss module provides routines to access various types of
file systems using a single interface.

Usage as Python module
======================

The module is imported using::

  import gss
  
Use the :py:meth:`List` routine to list the content of a directory::

  gss.List( '/data' )
  
In this example, the path specification describes a local file system.
To list the files on the ECMWF tape archive using the EcAccess tools, use::

  gss.List( 'ec:ec:/you/data' )

See the the section on :ref:`File specifications <filespecs>` 
for the supported file systems.

The following operations on files are supported:

* List content of path, return list with :py:class:GSS_Element' objects::

    gss.List( 'ec:scratch:/you/data/' )

* Copy file from the local file system to for example the ECMWF tape archive::

    gss.Put( 'test.txt', 'ec:ec:/you/data/test.txt' )

* Copy file from a remote location to the local file system::

    gss.Get( 'ec:scratch:/you/data/test.txt', 'test.txt' )

See the individual methods for further details.

.. _filespecs:

File specifications
===================

Examples of file specifications:

* Standard file system::

    [file[<options>]:][<path>/]<file>

  Supported options::
  
    %umask=022    # permission mask for created files
      
  Examples::
  
    data/file.txt
    file:data/file.txt
    file%umask=022:data/file.txt
    
* Remote file systems:

    <url>

  where the URL could be:
  
     https://some.archive.net/data/file.txt

* ECMWF file system accessed via EcAccess tools::

    ec[<options>]:[<domain>:][<path>/][<file>]

  Supported domains:
  
  * 'home'       : HOME partition on member state server (default)
  * 'scratch'    : SCRATCH partition on member state server
  * 'ec'         : tape archive
  * 'ectmp'      : temporary tape archive 

  Supported options::
  
    umask=022    # permission mask for created files and directories
      
  Examples::

    ec:data/file.txt                  # file on HOME partition
    ec:home:data/file.txt             # "    "  "    "
    ec:scratch:data/file.txt          # file on SCRATCH partition
    ec:ec:data/file.txt               # file on tape archive
    ec:ectmp:data/file.txt            # file on temporary tape archive
    ec:ec:/TMP/you/data/file.txt      # "    "  "         "    "
    ec%umas=022:ec:data/file.txt      # file on tape archive, read permissions


"""

# ======================================================================
# ===
# === GSS_Element
# ===
# ======================================================================

class Element( object ) :

    """
    Base class for GSS elements, which is a file or directory.
    
    Attributes:
    
    * ftype   :  character: f=file, d=directory
    * name
    
    """
    
    def __init__( self, name, ftype='f' ) :
    
        """
        Define element.
        """
        
        # store:
        self.name  = name
        self.ftype = ftype
        
    #enddef __init__
    
    # *
    
    def Print( self, long=False, indent='' ) :
   
        """
        Pretty print of content.
        """
        
        # format:
        if long :
            line = indent+'%s %s' % (self.ftype,self.name)
        else :
            line = indent+self.name
        #endif
        
        # show:
        print( line )
        
    #enddef Print

    
#endclass Element


# ======================================================================
# ===
# === GSS_Base
# ===
# ======================================================================

class GSS_Base( object ) :

    """
    Base class for GSS objects.
    """
    
    def __init__( self, umask=None, verbose=False, indent='' ) :
    
        """
        Initialize base object for GSS classes.
        
        Optional arguments:
        
        * verbose  : (bool) set to True for messages
        * indent   : (str of whitespace) initial indent in messages
        """
        
        # store:
        self.verbose = verbose
        self.indent  = indent
        self.umask   = umask
        
    #enddef __init__
    
    # *
    
    def info( self, msg ) :
    
        """
        Print message if verbose mode is on.
        Optional indent is preceeded to message.
        """

        # show messages?
        if self.verbose :
            # write message:
            print( self.indent+msg )
        #endif
        
    #enddef info
    
    # *
    
    def GetMode( self, directory=False ) :
    
        """
        Convert umask attribute to mode, or None if no umask was defined.
        The mode is formed by subtracting the umask from default creation modes:
        
        * For files the default creation mode is '666' (readible/writable for all).
          With umask '022', the actual creation mode is then '644' (readible for all,
          writable for user only).
        
        * For files the default creation mode is '777' (readible/writable/executable for all).
          With umask '022', the actual creation mode is then '755' (readible and
          executable for all, writable for user only).
          
        By default the mode for a file is returned, unless the optional
        argument 'directory' is True.
        """
        
        # not defined?
        if self.umask is None :

            # no mode:
            mode = None

        else :

            # default mode:
            if directory :
                default = 7
            else :
                default = 6
            #endif

            # init result:
            mode = ''
            # loop over mask elements:
            for c in self.umask :
                # convert:
                u = int(c)
                # add:
                mode = mode+str(default-u)
            #endfor
            
        #endif
        
        # ok
        return mode
        
    #enddef GetMode
    
#endclass GSS_Base


# ======================================================================
# ===
# === GSS_File
# ===
# ======================================================================

class GSS_File( GSS_Base ) :

    """
    GSS class to access standard file system.
    """
    
    def List( self, path ) :
    
        """
        List directory content or specific file.
        Returns a list with :py:class:GSS_Element' objects.
        """
        
        # modules:
        import os
        
        ## check ...
        #if not os.path.isdir( path ) :
        #    print( 'ERROR - path not found: %s' % path )
        #    raise Exception
        ##endif
        
        # init result:
        elements = []

        # info ...
        self.info( 'file list: %s' % path )
        # get files:
        if os.path.isdir( path ) :
            fnames = os.listdir( path )
        elif os.path.isfile( path ) :
            fnames = [path]
        else :
            fnames = []
        #endif
        # loop:
        for fname in fnames :
            # full path:
            filename = os.path.join( path, fname )
            # set type:
            if os.path.isdir( filename ) :
                ftype = 'd'
            else :
                ftype = 'f'
            #endif
            # store:
            elements.append( Element( filename, ftype=ftype ) )
        #endfor
        
        # ok
        return elements
        
    #enddef List
    
    # *
    
    def IsFile( self, path ) :
      
        """
        Return True if path describes a file.
        """
        
        # modules:
        import os
        
        # check:
        return os.path.isfile(path)
      
    #enddef IsFile
    
    # *
    
    def MakeDirs( self, path ) :
    
        """
        Make directory and subdirectories.
        """
        
        # modules:
        import os

        # create if necessary:
        if not os.path.isdir(path) :
            # info ..
            self.info( 'create directory: %s' % path )
            # get creation mode:
            mode = self.GetMode( directory=True )
            # with mode?
            if mode is not None :
                # create recursively, convert mode from octal to integer:
                os.makedirs( path, mode=eval('0o'+mode) )
            else :
                # create:
                os.makedirs( path )
            #endif
        #endif
        
    #enddef MakeDirs

    # *

    def Put( self, source, target ) :

        """
        Copy source file from local file system to target location.
        """

        # modules:
        import os
        import subprocess
        import shutil
        
        # directory part:
        dname = os.path.dirname( target )
        # any path in target?
        if len(dname) > 0 :
            # create if necessary:
            self.MakeDirs( dname )
        #endif
        
        # info ...
        self.info( 'copy %s to %s ...' % (source,target) )
        # copy:
        shutil.copy( source, target )

        # get creation mode, None if not defined:
        mode = self.GetMode()
        # change mode?
        if mode is not None :
            # set mode:
            subprocess.check_call(['chmod',mode,target])
        #endif
        
    #enddef Put

    # *

    def Get( self, source, target ) :

        """
        Copy source file from local file system to target location.
        """
        
        # copy:
        self.Put( source, target )
        
    #enddef Put

#endclass GSS_File


# ======================================================================
# ===
# === WGet
# ===
# ======================================================================

class GSS_WGet( GSS_Base ) :

    """
    GSS class to access remote file system through the :py:mod:`wget` module.
    """
    
    def __init__( self, protocol, **kwargs ) :
    
        """
        Initialize remote file, protocol ('http', 'https') needed to format file names.
        """
        
        # init base class:
        GSS_Base.__init__( self, **kwargs )
        
        # store protocl, add colon for quick combination with path into url:
        self.protocol = protocol+':'
    
    #enddef __init__
    
    # *
    
    def IsFile( self, path ) :
      
        """
        Return True if path describes a file.
        For this class the method simply checks if the path can be openend,
        if an HTTP error is raised this is interpreted as that the path is not present.
        The method can also not distinguish between a file and a directory yet.
        """
        
        # modules:
        import urllib.request
        import urllib.error
        
        # try to open file, fails if not present:
        try :
            # try to open:
            urllib.request.urlopen( self.protocol+path )
            # ok:
            return True
        except urllib.error.HTTPError :
            # something wrong, assume that means file not present ..
            return False
        except :
            # unknown error, raise it and break:
            raise
        #endtry
        
    #enddef IsFile
    
    # *
    
    def Get( self, source, target ) :

        """
        Get file from remote file system and store under target name.
        """

        # modules:
        import os
        import subprocess
        
        # no standard pacakge ..
        try :
            import wget
            external_wget = False
        except :
            #print( 'ERROR - ' )
            #print( 'ERROR - Could not import "wget" module from python distribution.' )
            #print( 'ERROR - When using an Anaconda distribution, try:' )
            #print( 'ERROR - ' )
            #print( 'ERROR -    conda install -c bioconda python-wget' )
            #print( 'ERROR - ' )
            #raise Exception
            external_wget = True
        #endtry
        
        # info ..
        self.info( 'download ..' )
        self.info( '  source : %s' % self.protocol+source )
        self.info( '  target : %s' % target )
        
        # directory part:
        dname = os.path.dirname( target )
        # any path in target?
        if len(dname) > 0 :
            # create if necessary:
            if not os.path.isdir(dname) :
                # info ..
                self.info( '    create directory: %s' % dname )
                # create:
                os.makedirs( dname )
            #endif
        #endif
        
        # external command or module?
        if external_wget :
            # info ...
            self.info( '  download using external "wget" command ...' )
            # subdir?
            if len(dname) > 0 :
                # current location:
                cwd = os.getcwd()
                # change:
                os.chdir( dname )
            #endif
            # long messages with error bar ...
            wout = 'wget.out'
            werr = 'wget.err'
            self.info( '    (output redirected to %s, error to %s)' % (wout,werr) )
            # use wget from os, trap errors:
            try :
                # redirect progres bar ...
                subprocess.check_call( 'wget %s > %s 2> %s' % (self.protocol+source,wout,werr), shell=True )
            except :
                # remove (partly) file if present:
                if os.path.isfile(target) : os.remove( target )
                # show output/error:
                for wfile in [wout,werr] :
                    if os.path.isfile(wfile) :
                        f = open( wfile, 'r' )
                        lines = f.readlines()
                        f.close()
                        for line in lines : self.logger.error( line )
                        os.remove( wfile )
                    #endif
                #endfor
                # raise error:
                raise
            #endtry
            # cleanup:
            if os.path.isfile(wout) : os.remove(wout)
            if os.path.isfile(werr) : os.remove(werr)
            # subdir?
            if len(dname) > 0 :
                # back:
                os.chdir( cwd )
            #endif
        else :
            # info ...
            self.info( '  download ...' )
            # use the 'wget' module to download,
            # do not show a progres bar:
            wget.download( self.protocol+source, out=target+'.tmp', bar=None )
            # rename:
            os.rename( target+'.tmp', target )
        #endif
        
    #enddef Get

#endclass GSS_WGet


# ======================================================================
# ===
# === EcAccess
# ===
# ======================================================================

class GSS_EcAccess( GSS_Base ) :

    """
    GSS class to access ECMWF file systems through the EcAccess tools.
    """
    
    def _Call( self, command, retry=5, wait=5, **kwargs ) :
    
        """
        Perform ecaccess command.
        Optionally retry until succeed, wait for a number of seconds between two attempts.
        
        Return values:
        * stdout line (incl newline characters)
        * stderr line (incl newline characters)
        * return code, 0 if no error, <0 indicates warning, >0 error
        
        Known error codes and messages;
        some are translated into warning code:
        * 255 "No such file or directory" ; warning code -1
        * 255 "Unknown host: xscratch"
        * 255 "ECFS not available (sleep)"
        
        """
        
        # modules:
        import os
        import subprocess
        import time
        
        # info ...
        self.info( 'call command: %s' % str(command) )
        
        # run:
        attempt = 1
        while True :
        
            # info ...
            self.info( '  attempt %i ...' % attempt )

            # run and obtain output and error:
            p = subprocess.Popen( command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs )
            stdout,stderr = p.communicate()
            # info ...
            self.info( '    std. out.   : %s' % str(stdout) )
            self.info( '    std. err.   : %s' % str(stderr) )
            self.info( '    return code : %i' % p.returncode )
            # check return code:
            if p.returncode == 0 :
                # info ...
                self.info( '    return code ok, leaving ...' )
                # ok, leave repeat loop:
                retcode = 0
                break
            #
            elif stderr.startswith('No such file or directory') :
                # info ...
                self.info( '    leave with warning code ...' )
                # leave with warning code:
                retcode = -1
                break
            #
            elif stderr.startswith('ECFS not available') :
                # info ...
                self.info( '    ECFS not available, try again ...' )
            #
            elif attempt >= retry :
                # info ...
                self.info( '    return code not ok, but number of attempts reached maximum ...' )
                # extend error message:
                stderr = stderr+('\nGSS ERROR - tried %i times to run command:' % attempt)
                stderr = stderr+('\nGSS ERROR -   %s' % str(command))
                # return with error code:
                retcode = p.returncode
                break
            #endif
            
            # info ..
            self.info( '    try gain, first wait for %i seconds ...' % wait )
            # wait for a number of seconds ..
            time.sleep(wait)
            # try again ...
            attempt = attempt + 1
            
        #endwhile # attempts
    
        # ok
        return stdout,stderr,retcode
        
    #enddef _Call
    
    # *
    
    def List( self, path ) :
    
        """
        List directory content.
        Returns a list with :py:class:`Element` objects,
        or None if the path does not exist.
        
        Examples of accepted input paths:
        
        * [home:]path/to/data
        * scratch:path/to/data
        * ec:path/to/data
        * ectmp:path/to/data
        
        """
        
        # modules:
        import os
        
        # split in domain and remaining path:
        if ':' in path :
            domain,dpath = path.split(':',1)
        else :
            domain,dpath = 'home',path
        #endif
        
        # info ...
        self.info( 'ecaccess list: %s:%s' % (domain,dpath) )
        # listing command:
        command = ['ecaccess-file-dir','-long',path]
        # run, return code -1 for 'file not found':
        stdout,stderr,retcode = self._Call( command )
        # check return code:
        if retcode == -1 :
            # file or directory not found:
            return None
        elif retcode != 0 :
            print( stdout )
            print( stderr )
            print( 'ERROR from command: %s' % command )
            raise Exception
        #endif
        # init result:
        elements = []
        # loop:
        for line in stdout.split('\n') :
            # cleanup:
            line = line.strip()
            # skip empty:
            if len(line) == 0 : continue
            # format of listing:
            #   -rw-r--r--   1 nlx      nl       44077       Feb 23 15:07 file.txt
            # split:
            fields = line.split()
            # extract type:
            ftype = fields[0][0]
            if ftype == '-' : ftype = 'f'
            # extract name:
            fname = fields[-1]
            # full path on domain:
            filename = os.path.join( dpath, fname )
            # store:
            elements.append( Element( filename, ftype=ftype ) )
        #endfor
        
        # ok
        return elements
        
    #enddef List
    
    # *
    
    def IsFile( self, path ) :
      
        """
        Return True if path describes a file.
        """
        
        # list file; return None in case path is not valid (or other error):
        elements = self.List( path )
        # return values:
        #   None        : some error (file not found?)
        #   []          : empty directory
        #   [elem]      : single file
        #   [elem,...]  : directory content
        # not found, or other error?
        if elements is None :
            # no file:
            isfile = False
        else :
            # check length:
            isfile = len(elements) == 1
        #endif
        
        # ok
        return isfile
      
    #enddef IsFile
    
    # *
    
    def MakeDirs( self, path ) :
    
        """
        Make directory and subdirectories.
        """
        
        # modules:
        import os
        
        # split in domain and remaining path:
        if ':' in path :
            domain,dpath = path.split(':',1)
        else :
            domain,dpath = 'home',path
        #endif
        # remove trailing '/' if present:
        dpath = dpath.rstrip('/')
        
        # first try to list full path, if this fails try parent, etc.
        # keep list of subdirs to be created:
        subdirs = []
        # init parents directory as full path:
        ppath = dpath
        # loop:
        while len(ppath) > 0 :
            # loop until directory can be listed or if it is not present;
            # list command:
            command = ['ecaccess-file-dir','%s:%s' % (domain,ppath)]
            # try to list, return code -1 for 'file not found':
            stdout,stderr,retcode = self._Call( command )
            # no failure if exists ...
            if retcode == 0 :
                # directory was found, no new mkdirs needed; leave:
                break
            elif retcode == -1 :
                # any subdirs left?
                if '/' in ppath :
                    # split into new parent directory
                    ppath,subdir = ppath.rsplit('/',1)
                    # prepend in list:
                    subdirs = [subdir]+subdirs
                else :
                    # no subdirs left, the current ppath is a subdir to be created too:
                    subdirs = [ppath]+subdirs
                    # empty:
                    ppath = ''
                #endif
            else :
                # some error ...
                print( stdout )
                print( stderr )
                print( 'ERROR from command: %s' % command )
                raise Exception
            #endif
        #endwhile
        # create subdirs?
        if len(subdirs) > 0 :
            # get creation mode, None if not defined:
            mode = self.GetMode( directory=True )
            # loop:
            for subdir in subdirs :
                # new path:
                if len(ppath) > 0 : ppath = ppath+'/'
                ppath = ppath+subdir
                # create:
                command = ['ecaccess-file-mkdir','%s:%s' % (domain,ppath)]
                stdout,stderr,retcode = self._Call( command )
                # error?
                if retcode != 0 :
                    print( stdout )
                    print( stderr )
                    print( 'ERROR from command: %s' % command )
                    raise Exception
                #endif
                # change mode?
                if mode is not None :
                    # set mode:
                    command = ['ecaccess-file-chmod',mode,'%s:%s' % (domain,ppath)]
                    stdout,stderr,retcode = self._Call( command )
                    # error?
                    if retcode != 0 :
                        print( stdout )
                        print( stderr )
                        print( 'ERROR from command: %s' % command )
                        raise Exception
                    #endif
                #endif
            #endfor # subdirs
        #endif # create subdirs
        
    #enddef MakeDirs

    # *

    def Put( self, source, target ) :

        """
        Copy source file to target location at ECMWF file system.
        """

        # modules:
        import os
        
        # directory part:
        dname = os.path.dirname( target )
        # create if necessary:
        self.MakeDirs( dname )
        
        # info ...
        self.info( 'copy %s to ecaccess location %s ...' % (source,target) )
        # put:
        command = [ 'ecaccess-file-put', source, target ]
        stdout,stderr,retcode = self._Call( command )
        # error?
        if retcode != 0 :
            print( stdout )
            print( stderr )
            print( 'ERROR from command: %s' % command )
            raise Exception
        #endif

        # get creation mode, None if not defined:
        mode = self.GetMode()
        # change mode?
        if mode is not None :
            # set mode:
            command = ['ecaccess-file-chmod',mode,target]
            stdout,stderr,retcode = self._Call( command )
            # error?
            if retcode != 0 :
                print( stdout )
                print( stderr )
                raise Exception
            #endif
        #endif
        
    #enddef Put

    # *

    def Get( self, source, target ) :

        """
        Get file from ECMWF file system and store under target name.
        """

        # modules:
        import os
        
        # directory part:
        dname = os.path.dirname( target )
        # any path in target?
        if len(dname) > 0 :
            # create if necessary:
            if not os.path.isdir(dname) :
                # info ..
                self.info( 'create directory: %s' % dname )
                # create:
                os.makedirs( dname )
            #endif
        #endif
        
        # info ...
        self.info( 'copy file from ecaccess location %s to %s ...' % (source,target) )
        # temporary target:
        tmpfile = target+'.tmp'
        # run:
        command = ['ecaccess-file-get',source,tmpfile]
        stdout,stderr,retcode = self._Call( command )
        # error?
        if retcode != 0 :
            print( stdout )
            print( stderr )
            print( 'ERROR from command: %s' % command )
            raise Exception
        #endif
        # not found? something wrong ...
        if not os.path.isfile( tmpfile ) :
            print( 'ERROR - temporary copy "%s" not found, but no error was catched ...' % tmpfile )
            print( 'ERROR - command used: %s' % command )
            raise Exception
        #endif
        # rename:
        os.rename( tmpfile, target )
        
    #enddef Get

#endclass GSS_EcAccess


# ======================================================================
# ===
# === wrapper
# ===
# ======================================================================


def GSS( path, **kwargs ) :

    """
    Return GSS object that corresponds with the file system specification in the path,
    as well as the remaining path.
    
    The object is an instance of one of the following classes:
    
    * :py:class:`GSS_File` for standard file systems;
    * :py:class:`GSS_EcAccess` for an ECMWF file system.

    See the the section on :ref:`File specifications <filespecs>` for the supported 
    formats of the path.
    
    This class is only used internally by module routines such as
    :py:meth:`List` and :py:meth:`Put`.

    """
    
    # any prefix?
    if ':' in path :

        # split:
        domain,rpath = path.split(':',1)
        
        # options? add them to the keyword arguments:
        if '%' in domain :
            # split:
            elements = domain.split('%')
            # first is domain, rest are options:
            domain = elements[0]
            opts = elements[1:]
            # convert options to dictionairy:
            for opt in opts :
                # check ...
                if '=' not in opt :
                    print( 'ERROR - options should be key=value, found "%s"' % opt )
                    raise Exception
                #endif
                # split:
                key,val = opt.split('=')
                # store:
                kwargs[key] = val
            #endfor
        #endif

        # switch:
        if domain == 'file' :
            # standard file system:
            obj = GSS_File( **kwargs )
        #
        elif domain in ['http','https'] :
            # remote file system, provide protocol as first argument:
            obj = GSS_WGet( domain, **kwargs )
        #
        elif domain == 'ec' :
            # standard file system:
            obj = GSS_EcAccess( **kwargs )
        #
        else :
            print( 'ERROR - unsupported domain "%s" in "%s"' % (domain,path) )
            raise Exception
        #endif
    #
    else :
        # standard file system:
        obj = GSS_File( **kwargs )
        # same path:
        rpath = path
    #endif
    
    # return object and remaining path:
    return obj,rpath

#enddef GSS


# ======================================================================
# ===
# === routines
# ===
# ======================================================================

def List( path, **kwargs ) :

    """
    List files on directory path.
    
    Arguments:
    
    * path    : source path, eventually including gss prefixes
    * kwargs  : keyword arguments passed to :py:class:`GSS` class.
    
    Return values:
    
    * list of :py:class:`Element` objects.
    """
    
    # obtain GSS ojbect and remaining path:
    gss,rpath = GSS( path, **kwargs )

    # list:
    return gss.List( rpath )
    
#enddef List

# *

def IsFile( path, **kwargs ) :
  
    """
    Return True if source describes a file.
    """
    
    # obtain GSS ojbect and remaining path:
    gss,rpath = GSS( path, **kwargs )

    # test:
    return gss.IsFile( rpath )
    
#enddef IsFile

# *

def MakeDirs( path, **kwargs ) :
  
    """
    Create directory including subdirectories.
    """
    
    # obtain GSS ojbect and remaining path:
    gss,rpath = GSS( path, **kwargs )

    # create:
    gss.MakeDirs( rpath )
    
#enddef MakeDirs

# *

def Put( source, target, **kwargs ) :

    """
    Copy source file to target location.
    """
    
    # modules:
    import os
    
    # check ..
    if not os.path.isfile(source) :
        print( 'ERROR - file not found: %s' % source )
        raise Exception
    #endif
    
    # obtain GSS object for target file system and return and remaining path:
    gss,rtarget = GSS( target, **kwargs )
    
    # copy to storage system:
    gss.Put( source, rtarget )
    
#enddef Put

# *

def Get( source, target, **kwargs ) :

    """
    Copy source file from source location to local target.
    """
    
    # modules:
    import os
    
    # obtain GSS object for source file system and return and remaining path:
    gss,rsource = GSS( source, **kwargs )
    
    # copy from storage system:
    gss.Get( rsource, target )
    
#enddef Get


# ======================================================================
# ===
# === end
# ===
# ======================================================================

