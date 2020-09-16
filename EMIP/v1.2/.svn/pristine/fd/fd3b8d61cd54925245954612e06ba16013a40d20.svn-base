
"""
.. Label, use :ref:`text <label>` for reference
.. _catalogue-module:

********************
``catalogue`` module
********************

Python module to create html catalogue to browse through a collection of images.

Tutorial
========

In this example, a catalogue is build for an archive
of meteorological images and statistical data::

   /data/fig/world-pressure-map.png
             world-pressure-time.png
             world-pressure-stats.txt
             :

To browse through the images, a catalogue in html format
should be created. Since the image names consist of 3 parts
(here seperated by dashes), the catalogue should have 3 levels too.
The ordering of the levels determines how a user could browse
through the archive.

In this example, the order of the levels is:

1. the variables (pressure, temperature, humidity)
   presented as a bullet list; each item is a link to a new page;
2. the region (world, europe)
   presented as a row in a table;
3. the image types (map, time series, ...)
   presented as a column in the table.

To create a catlogue, first define the main object.
The first argument describes the image filenames,
other information is optional::

  cata = catalogue.CatalogueCreator( '<region>-<variable>-<figure>.png',
                                     imgdir='/data/fig', 
                                     title='Example catalogue',
                                     longnames={ 'key' : 'value', ..},
                                     verbose=False )

The supplied filename is either a single template::

 '<region>-<variable>-<figure>.png'

or a list of templates::

 ['<region>-<variable>-<figure>.png','<region>-all.png',...]

Filetypes that are supported are currently:

* images: '``*.png``', '``*.gif``', ..

  These will be displayed by the html img tag.

* text files:  '``*.txt``'

  The content of the file is displayed as text.

* html files: '``*.html``'

  A link to the file is displayed.
  This is typically used to create a main page with links to various catalogues.

Then define the levels in the catalogue.
Each level is defined by 4 variables.

* A '``name``' ; 
  this name is used in the 'imgname' supplied above by enclosing it with '<..>'

* The '``values``' are a list with all possible values
  that a level could take.

* The '``form``' describes how the values in a level are listed
  in the catalogue, e.g. the html tags used to display them.
  Supported forms are :

  * '``ul``'        : unordered list
  * '``table-tr``'  : row within table, level value is in first column
  * '``table-td``'  : column within table with single row only
  * '``tr``'        : row within table, level value is not displayed
  * '``td``'        : column within table row; should follow '``table-tr``' or '``tr``'

* The '``newpage``' is a boolean that should be set to '``True``' to have each 
  item in this level point to a new html page.

Add the levels via::

    # define the levels in the right order:
    #              name        values                              form  newpage 
    cata.AddLevel( 'variable', ['pressure','temperature','wind'] , 'ul', True   )
    cata.AddLevel( 'region'  , ['world','europe']                , 'tr', False  )
    cata.AddLevel( 'figure'  , ['map.png','time.png','stats.txt'], 'td', False  )

    # optional arguments:
    #    longnames = { 'key' : 'value', ... }  # used for long item names

Optionally, define some html templates to fine-tune the formatting of
the images and/or other data.
Specify one html template for each value of the final level, 
(usually a 'figure' level as used here),
or use a filename pattern.
Include in the template the the special keys:

* '``%(filename)s``'   to insert the filename;
* '``%(filetext)s``'   to insert the content of a text file

Add templates by::

    cata.AddTemplate( 'map.png' , '<img src="%(filename)s" width="400" border="0">' )
    cata.AddTemplate( 'time.png', '<img src="%(filename)s" width="300" border="0">' )
    cata.AddTemplate( '*.txt'   , '<pre>%(filetext)s</pre>' )

Eventually add extra longnames that are used in listings and headers::

    cata.AddLongnames( {...} )

Finaly create the html pages using::

    cata.Create( basename='index', silent=True )

Multiple archives could be made by defining a different
level order and/or a different presentation formats
  

Command line tool
=================

Use the script::

    catalogue-index

to create a new catalogue using a configuration file.
Provide '``--help``' as argument to see the options and
how to generate example settings.


Classes
=======


"""

# *

class _Level( object ) :

    """
    Storage for calogue level definition, e.g. the level name,
    posible values, html tag form, etc.
    """

    def __init__( self, name, values, form, newpage, longnames={}, logger=None ) :

        """
        Level name name and possible values.
        """
        
        # external:
        import logging

        # logger:
        if logger is None :
            self.logger = logging.getLogger()
        else :
            self.logger = None
        #endif
    
        # store:
        self.name      = name
        self.values    = values
        self.form      = form
        self.newpage   = newpage
        self.longnames = longnames
        
        # info ...
        self.logger.debug( 'added catalogue level "%s" with values : %s' % (name,str(values)) )
    
    #enddef __init__

    # *
    
    def Get( self, ivalue ) :
    
        """
        Return value and longname (None if not defined).
        """
        
        # current:
        value = self.values[ivalue]
        # extra?
        if value in self.longnames.keys() :
            longname = self.longnames[value]
        else :
            longname = None
        #endif
        
        # ok
        return value,longname
    
    #enddef Get

#endclass _Level

# *

class CatalogueCreator( object ) :

    """
    Catalogue object.

    Arguments:

    * '``imgnames``'   : (list of) template image name(s)

    Optional arguments:

    * '``imgdir``'
    * '``title``'
    * '``longnames``'  : dictionairy with long names for catalogue values
    * '``logger``'      :  logger from calling method, otherwise new logger is configured
    * '``verbose``'    :  if new logger defined, used to set the logging level
    """

    def __init__( self, imgnames, imgdir='', title='catalogue', longnames={}, 
                          logger=None, verbose=False ) :

        # modules:
        import sys
        import os
        import logging
        
        # setup logger:
        if logger is None :
            # logger:
            logging.basicConfig( stream=sys.stdout, level=logging.INFO,
                                   format='[%(levelname)-8s] %(message)s' )
            self.logger = logging.getLogger()
            # logging level:
            if verbose : self.logger.setLevel( logging.DEBUG )
        else :
            self.logger = logger
        #endif
        
        # store title:
        self.title = title
        # store long names:
        self.longnames = longnames
        
        # location of images:
        self.imgdir = imgdir
        # create?
        if len(self.imgdir) > 0 :
            if not os.path.isdir(self.imgdir) : os.makedirs(self.imgdir)
        #endif
        
        # image names, always store as list:
        if type(imgnames) == list :
            self.imgnames = imgnames
        else :
            self.imgnames = [imgnames]
        #endif
        
        # init storage:
        self.levels = []
        self.templates = {}

    #enddef
    
    # *
    
    def AddLevel( self, name, *args, **kwargs ) :
        """
        Add a new level definition.
        """
        # add:
        self.levels.append( _Level( name, *args, **kwargs ) )
    #enddef
    
    # *
    
    def AddTemplate( self, name, template ) :
        """
        Add a new template for formatting images etc.
        """
        # add:
        self.templates[name] = template
    #enddef

    # *
    
    def AddLongNames( self, longnames ) :
    
        """
        Extend longnames dictionairy.
        """
        
        # loop:
        for key in longnames.keys() :
            self.longnames[key] = longnames[key]
        #endfor
        
    #enddef AddLongNames

    # *
    
    def WritePage( self, filename, lines ) :

        """
        Write a new html page.
        """

        # external:
        import os

        # full path:
        fname = os.path.join(self.imgdir,filename)

        # create directory if necessary:
        dname = os.path.dirname( fname )
        if len(dname) > 0 :
            if not os.path.isdir(dname) : os.makedirs(dname)
        #endif

        # write:
        f = open( fname, 'w' )
        for line in lines : f.writelines( line+'\n' )
        f.close()

    #enddef WritePage
    
    # *
    
    def Create( self, basename='index', silent=False, indent='', footer=None ) :

        """
        Create catalogue.
        
        Optional arguments:
        
        * '``basename``' : Base name of main index file, probably 'index'.
        * '``silent``' : If True, do not shout location of created index file.
        """

        # external:
        import sys
        import os

        # store:
        self.indent = indent

        # generate list for first keyword:
        lines = []
        lines.append('<html>')
        lines.append('<head><title>%s</title></head>' % self.title )
        lines.append('<body>')
        lines.append('<h2>%s</h2>' % self.title )
        lines = lines + self.CreateList( [], self.levels, basename, indent=self.indent )
        if footer is not None : lines.extend( footer )
        lines.append('</body>')
        lines.append('</html>')

        # write:
        indexfile = '%s.html' % basename
        self.WritePage( indexfile, lines )
        
        # absolute path to main file:
        filename = os.path.join( os.getcwd(), self.imgdir, indexfile )
        ## use $PWD to have logical path instead of physical:
        #filename = os.path.join( os.environ['PWD'], self.imgdir, indexfile )
        

        # info ...
        if not silent :
            self.logger.info( self.indent+'' )
            self.logger.info( self.indent+'Point your browser to :' )
            self.logger.info( self.indent+'    file://%s' % filename )
            self.logger.info( self.indent+'' )
        #endif

    #enddef Create
    
    # *
    
    def CreateList( self, pvals, levels, itembase_up, indent='' ) :

        """
        Return list for this keyword.
        """

        # external:
        import os
        import fnmatch

        # first:
        level = levels[0]

        # number of values:
        nvalue = len(level.values)

        # storage for items in this level:
        items = []

        # loop over possible keyword values:
        for ivalue in range(nvalue) :
            # current:
            value,longname = level.Get( ivalue )
            # basename, used for name of html file:
            itembase = ''
            for pval in pvals+[value] : itembase = itembase+'-'+str(pval)
            if itembase.startswith('-') : itembase = itembase[1:]
            # info ...
            self.logger.debug( indent+itembase )
            # resolve inner loops:
            if len(levels) == 1 :
                # at this level the image will be included,
                # no lower levels: 
                sublines = []
            else :
                # current html file:
                if level.newpage :
                    itembase_up_curr = itembase
                else :
                    itembase_up_curr = itembase_up
                #endif
                # recursive call:
                sublines = self.CreateList( pvals+[value], levels[1:], itembase_up_curr, indent=indent+'  ' )
                # empty ? not necessary to add to index, try next value:
                if len(sublines) == 0 : continue
            #endif
            # store:
            items.append( { 'value'    : value, 
                            'longname' : longname,
                            'sublines' : sublines, 
                            'itembase' : itembase } )
        #endfor

        # init output:
        lines = []
        endlines = []
        anyfound = False
        # loop over possible keyword values:
        for iitem in range(len(items)) :

            # current:
            value    = items[iitem]['value']
            longname = items[iitem]['longname']
            sublines = items[iitem]['sublines']
            itembase = items[iitem]['itembase']

            # previous:
            if iitem > 0 :
                value_prev = items[iitem-1]['value']
            else :
                value_prev = None
            #endif
            # next:
            if iitem < len(items)-1 :
                value_next = items[iitem+1]['value']
            else :
                value_next = None
            #endif

            # other base names:
            if value_prev != None :
                itembase_prev = itembase_up
                for pval in pvals+[value_prev] : itembase_prev = itembase_prev+'-'+str(pval)
                if itembase_prev.startswith('-') : itembase_prev = itembase_prev[1:]
            #endif
            if value_next != None :
                itembase_next = itembase_up
                for pval in pvals+[value_next] : itembase_next = itembase_next+'-'+str(pval)
                if itembase_next.startswith('-') : itembase_next = itembase_next[1:]
            #endif

            # text in index: either image or level values:
            if len(levels) == 1 :

                # current value list:
                vals = pvals+[value]

                # try all image names until first match:
                tried = []
                found = False
                for imgfile in self.imgnames :
                    # replace keywords:
                    for ilev in range(len(vals)) :
                        imgfile = imgfile.replace( '<%s>' % self.levels[ilev].name, str(vals[ilev]) )
                    #endfor
                    # store:
                    tried.append( imgfile )
                    # match?
                    found = os.path.exists(os.path.join(self.imgdir,imgfile))
                    # found ? then leave:
                    if found : break
                #endfor
                # not found ?
                if found :
                    # info ...
                    self.logger.debug( indent+'    found     : %s' % imgfile )
                    # reset flag:
                    anyfound = True
                else :
                    # show all:
                    for tfile in tried :
                        self.logger.debug( indent+'    not found : %s' % tfile )
                    #endfor
                    # reset:
                    imgfile = None
                #endif

                # check for filename match:
                xvalue = value
                if imgfile is not None :
                    for pat in self.templates.keys() :
                        if fnmatch.fnmatch(imgfile,pat) :
                            xvalue = pat
                            break
                        #endif
                    #endfor
                #endif   # with imgfile
                
                # html template present ?
                if xvalue in self.templates.keys() :
                    # extract:
                    template = self.templates[xvalue]
                else :
                    # empty?
                    if imgfile is None :
                        # space:
                        template = '&nbsp;'
                    else :
                        # set a default template:
                        if imgfile.endswith('.txt') :
                            # preformatted text, file content is inserted below:
                            template = '<pre>%(filetext)s</pre>'
                        elif imgfile.endswith('.html') :
                            # insert root of filename:
                            template = '<pre>%(rootname)s</pre>'
                        else :
                            # image:
                            template = '<img src="%(filename)s" border="0">'
                        #endif
                    #endif
                #endif  # template defined
                # insert filename if necessary:
                if ('%(filename)s' in template) and (imgfile is not None) :
                    # insert:
                    itemtext = template.replace('%(filename)s',imgfile)
                #endif
                # insert root of filename if necessary:
                if ('%(rootname)s' in template) and (imgfile is not None) :
                    # longname?
                    if value in self.longnames.keys() :
                        itemtext = template.replace('%(rootname)s',self.longnames[value])
                    else :
                        root,ext = os.path.splitext( os.path.basename(imgfile) )
                        itemtext = template.replace('%(rootname)s',root)
                    #endif
                #endif
                # insert file content if necessary:
                if '%(filetext)s' in template :
                    # read text file:
                    f = open(os.path.join(self.imgdir,imgfile),'r')
                    filetext = f.read()
                    f.close()
                    # insert:
                    itemtext = template.replace('%(filetext)s',filetext )
                #endif
                # convert into a link:
                if imgfile is None :
                    itemtext = None
                else :
                    # link tag:
                    itemtext = '<a href="%s">%s</a>' % (imgfile,itemtext)
                    # add longname?
                    if longname is not None : itemtext = itemtext+' - '+longname
                #endif
                # a new page for a single figure does not make sense ..
                newpage = False

            else :

                # form text :  <level1> / <level2> / ... / <level>
                itemtext_long = ''
                for pval in pvals+[value] :
                    if len(itemtext_long) > 0 : itemtext_long = itemtext_long+' / '
                    long_pval = pval
                    if pval in self.longnames.keys() : long_pval = self.longnames[pval]
                    itemtext_long = itemtext_long+str(long_pval)
                #endfor
                # ... or just the current level:
                itemtext = value
                # reset flag:
                anyfound = True
                # copy flag
                newpage = level.newpage

            #endif  # last level

            # write to subfile ?
            if newpage :
                # name of sub file:
                subfile = '%s-%s.html' % (itembase_up,itembase)
                # content:
                slines = []
                slines.append('<html>')
                slines.append('<body>')
                # navigation:
                slines.append( '<a href="%s.html">[up]</a>&nbsp;&nbsp;' % itembase_up )
                if value_prev == None :
                    slines.append('[prev]' )
                else :
                    slines.append('<a href="%s.html">[prev]</a>' % itembase_prev )
                #endif
                if value_next == None :
                    slines.append('[next]')
                else :
                    slines.append('<a href="%s.html">[next]</a>\n' % itembase_next )
                #endif
                slines.append( '<hr>' )
                # header:
                slines.append( '<h2>%s</h2>' % itemtext_long )
                # content:
                slines = slines + sublines
                # navigation:
                slines.append( '<hr>' )
                slines.append( '<a name="footer">' )
                slines.append( '<a href="%s.html#footer">[up]</a>&nbsp;&nbsp;' % itembase_up )
                if value_prev == None :
                    slines.append('[prev]' )
                else :
                    slines.append('<a href="%s.html#footer">[prev]</a>' % itembase_prev )
                #endif
                if value_next == None :
                    slines.append('[next]')
                else :
                    slines.append('<a href="%s.html#footer">[next]</a>\n' % itembase_next )
                #endif
                # footer:
                slines.append('</body>')
                slines.append('</html>')
                # write:
                self.WritePage( subfile, slines )
                # make content empty:
                sublines = []
                # reset item text to link:
                if itemtext is not None :
                    # link tag:
                    itemtext = '<a href="%s">%s</a>' % (subfile,itemtext)
                    # add longname?
                    if longname is not None : itemtext = itemtext+' - '+longname
                #endif
            #endif

            # how to list ?
            if level.form == 'ul' :
                if iitem == 0 : lines.append('<ul>')
                if itemtext is not None :
                    lines.append( '<li> %s' % itemtext )
                #endif
                lines = lines + sublines
                endlines = ['</ul>']
            elif level.form == 'table-tr' :
                if iitem == 0 : lines.append('<table border="0" cellpadding="5">')
                lines.append( '<tr>' )
                if itemtext is None :
                    lines.append( '<td> &nbsp; </td>' )
                else :
                    lines.append( '<td> %s </td>' % itemtext )
                #endif
                lines = lines + sublines
                lines.append( '</tr>' )
                endlines = ['</table>']
            elif level.form == 'table-td' :
                if iitem == 0 : lines.append('<table border="0" cellpadding="5"><tr>')
                if itemtext is None :
                    lines.append( '<td> &nbsp; </td>' )
                else :
                    lines.append( '<td> %s </td>' % itemtext )
                #endif
                lines = lines + sublines
                endlines = ['</tr></table>']
            elif level.form == 'tr' :
                if iitem == 0 : lines.append('<table border="0" cellpadding="5">')
                lines.append( '<tr>' )
                lines = lines + sublines
                lines.append( '</tr>' )
                endlines = ['</table>']
            elif level.form == 'td' :
                if itemtext is None :
                    lines.append( '<td> &nbsp;' )
                else :
                    lines.append( '<td> %s' % itemtext )
                #endif
                lines = lines + sublines
                lines.append( '</td>' )
            else :
                self.logger.error( 'unsupported list form : %s' % level.form )
                raise ValueError
            #endif

        #endfor

        # add end lines:
        lines = lines+endlines
        
        # reset to empty if no content found:
        if not anyfound : lines = []

        # ok
        return lines

    #enddef
    
#endclass CatalogueCreator


