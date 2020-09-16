

"""
Plotting tools.

Class hierchy:

  Figure
    MapFigure
    ColorbarFigure
      ColorbarMapFigure

"""



########################################################################
###
### figure
###
########################################################################


class Figure( object ) :

    """
    Base class for figures.
    """

    def __init__( self, figsize=None ) :

        """
        Setup figure.
        
        Optional arguments:
          figsize=(8,6)   : size of the figure in inches

        """
        
        # modules:
        import matplotlib
        import matplotlib.pyplot

        # setup figure:
        self.fig = matplotlib.pyplot.figure( figsize=figsize )
        
        # no axes yet:
        self.ax  = None
        self.bmp = None
        
    #enddef __init__

    # *

    def Close( self ) :

        """
        Destroy figure.
        """
        
        # modules:
        import matplotlib
        import matplotlib.pyplot

        # destroy:
        matplotlib.pyplot.close( self.fig )

    #enddef  # Close

    # *

    def Export( self, filename ) :

        """
        Export figure.
        """

        # modules
        import os

        # create directory if necessary:
        dname = os.path.dirname(filename)
        if (len(dname) > 0) and (not os.path.isdir(dname)) : os.makedirs(dname)

        # export:
        self.fig.savefig( filename )

    #enddef  # Export
    
    # *
    
    def AddAxes( self, position=(0.02,0.02,0.96,0.96) ) :
    
        """
        Add ax to figure.
        
        Optional arguments:
            position=(left,bottom,width,height)
                Relative position of lower left corner and size;
                defaults values ensure maximum axis.
            
        """
        
        # add ax to figure:
        self.ax = self.fig.add_axes(position)

    #enddef AddAxes
    
    # *
    
    def Set_Domain( self, domain=None ) :
    
        """
        Set axis domain.
        Optional domain is (xmin,xmax,ymin,ymax) ;
        if undefined, no limits are set.
        """

        # axis defined ?
        if self.ax is not None :
            # domain defined ?
            if domain is not None :
                xmin,xmax,ymin,ymax = domain
                self.ax.set_xlim((xmin,xmax))
                self.ax.set_ylim((ymin,ymax))
            #endif
        #endif
        
    #enddef Set_Domain
    
    
    # *
    
    
    def Rectangle( self, domain, **kwargs ) :
    
        """
        Draw rectangle with corners specified by:
           domain = west,east,south,north
        """
        
        # axis defined ?
        if self.ax is not None :
            # extract:
            west,east,south,north = domain
            # add line:
            self.ax.plot( [west ,east ,east ,west ,west ],
                          [south,south,north,north,south], **kwargs )
        #endif
                      
    #enddef

#endclass Figure



########################################################################
###
### figure with geographical map
###
########################################################################


def Setup_BaseMap( ax, title=None, domain=None, raster=None,
                     coastlines=True, countries=False, rivers=False, 
                     **kwargs ) :

    """
    Setup geographical basemap in the 'map' field of the class.
    If no domain was specified, a global map is setup.

    Arguments:
      ax           :  figure axes to hold the map
      
    Optional arguments:
      title        :  title above map
      domain       :  (west,east,south,north), default globe
      raster       :  (meridian_distance,parallel_distance)
      coastlines   :  default True
      countries    :  True or False
      rivers       :  True or False

    Other kewyord arguments are passed directly to Basemap, e.g.:
      projection   :  'cyl'=cylindrical, 
      resolution   :  'c' crude, 'l' low, 'i' intermediate, 'h' high

    If no 'projection' is defined, a default cylindrical map is setup.

    Return basemap object.

    """

    # modules:
    try :
        import mpl_toolkits.basemap as bmp
    except :
        # From:
        #   https://github.com/conda-forge/basemap-feedstock/issues/30
        # Hack to fix missing PROJ4 env var
        import os
        # machine specific ...
        if os.path.isdir('/home/dirac/miniconda2/share/proj') :
            proj_lib = '/home/dirac/miniconda2/share/proj'
        else :
            import conda
            conda_file_dir = conda.__file__
            conda_dir = conda_file_dir.split('lib')[0]
            proj_lib = os.path.join( os.path.join(conda_dir,'share'), 'proj' )
        #endif
        # store in environment:
        os.environ["PROJ_LIB"] = proj_lib
        # try again ...
        import mpl_toolkits.basemap as bmp
    #endtry
    import numpy
    
    # set default domain ?
    if domain is None : domain=[-180.0,180.0,-90.0,90.0]

    # extract domain:
    west,east,south,north = domain

    # basemap settings:
    bkwargs = dict( projection='cyl' )
    bkwargs.update( kwargs )

    # setup map:
    projection = bkwargs['projection']
    if projection in ['cyl'] :
        # standard cylindrical projection:
        m = bmp.Basemap( ax=ax,
                            llcrnrlon=west, llcrnrlat=south,
                            urcrnrlon=east, urcrnrlat=north,
                            **bkwargs )
    elif projection in ['merc'] :
        # more equidistant
        m = bmp.Basemap( ax=ax,
                            llcrnrlon=west, llcrnrlat=south,
                            urcrnrlon=east, urcrnrlat=north,
                            lat_ts=0.5*(north+south),
                            **bkwargs )
    elif projection in ['tmerc'] :
        # more equidistant
        m = bmp.Basemap( ax=ax,
                            llcrnrlon=west, llcrnrlat=south,
                            urcrnrlon=east, urcrnrlat=north,
                            lon_0=0.5*(east+west),
                            lat_0=0.5*(north+south),
                            **bkwargs )
    else :
        # assume everything is defined in the kwargs:
        m = bmp.Basemap( ax=ax, **kwargs )
    #endif

    # add title:
    if title is not None : ax.set_title( title )

    # add coast lines:
    if coastlines : m.drawcoastlines()
    # add rivers ?
    if rivers : m.drawrivers()
    # add country borders ?
    if countries : m.drawcountries()

    # add raster ?
    if raster is not None :
        dmeri,dpara = raster
        meri0 = int( west/dmeri)*dmeri
        para0 = int(south/dpara)*dpara
        m.drawmeridians( meri0 + (numpy.arange((east -meri0)/dmeri+1))*dmeri, labels=[0,0,0,1] )
        m.drawparallels( para0 + (numpy.arange((north-para0)/dpara+1))*dpara, labels=[1,0,0,0] )
    #endif

    # ok
    return m

#enddef  # Setup_BaseMap

# *

class MapFigure( Figure ) :

    """
    Base class for figure with map.
    """

    def __init__( self, figsize=None,
                    position=None, bmp={} ) :

        """
        Setup figure.
        
        Optional arguments passed to 'Figure' initialization:
            figsize
          
        Optional arguments to setup data axes:
            position=(left,bottom,width,height)
                position of axes; defaults are guessed
            bmp
                define geograpical map, passed to 'Setup_BaseMap'

        """
        
        # init figure:
        Figure.__init__( self, figsize=figsize )
        
        # data axis:
        if position is None :
            # start with maximum ax position:
            left,bottom,width,height = 0.02,0.02,0.96,0.96
            # setup ax for basemap; raster labels needed ?
            if 'raster' in bmp.keys() :
                # remove part from left and bottom:
                dx,dy = 0.10,0.10
                left,width = left+dx,width-dx
                bottom,height = bottom+dy,height-dy
            #endif
            # some extra space for title:
            if 'title' in bmp.keys() : height = height - 0.06
            # combine:
            position = (left,bottom,width,height)
        #endif
        # add ax to figure:
        self.AddAxes( position )
        
        # replace ax by map:
        self.ax0 = self.ax
        self.ax = Setup_BaseMap( self.ax0, **bmp )
        
    #enddef __init__

#endclass MapFigure



########################################################################
###
### figure with colorbar
###
########################################################################


class ColorbarFigure( Figure ) :

    """
    Base class for figure with colorbar.
    """

    def __init__( self, figsize=None, 
                    position=None, domain=None, title=None,
                    cmap={}, cbar={} ) :

        """
        Setup figure with colorbar.
        
        Optional arguments passed to 'Figure' initialization:
            figsize
          
        Optional arguments to setup data axes:
            position=(left,bottom,width,height)
                position of axes; defaults is guesed
            domain=(xmin,xmax,ymin,ymax)
                axes limits passed to 'Figure.Set_Domain'
            title
                axes title
        
        Optional arguments that define colorbar:
            cmap   : defines colormap, passed to 'Set_ColorMap'
            cbar   : defines colorbar, passed to 'Set_ColorBar' (False to skip)
            
        """
        
        # modules:
        import matplotlib
        import matplotlib.pyplot
        
        # init figure:
        Figure.__init__( self, figsize=figsize )

        # with colorbar?
        self.with_colorbar = True
        if cbar == False : self.with_colorbar = False
        
        # with colorbar?
        if self.with_colorbar :
        
            # default direction:
            self.cbar_orientation = 'horizontal'
            # replace ?
            if 'orientation' in cbar.keys() : self.cbar_orientation = cbar['orientation']
        
            # color bar axes:
            if self.cbar_orientation == 'horizontal' :
                # colorbar position:
                cbar_rect = [0.10, 0.10, 0.80, 0.07]
            else :
                # colorbar position:
                cbar_rect = [0.80, 0.10, 0.07, 0.80]
            #endif
            # add axes for  colorbar:
            self.cbar_ax = self.fig.add_axes( cbar_rect )

        #endif # with colorbar

        # data axes:
        if position is None :
            # start with maximum ax positions, no colorbar:
            left,bottom,width,height = 0.02,0.02,0.96,0.94
            # which direction ?
            if self.with_colorbar :
                if self.cbar_orientation == 'horizontal' :
                    # remove part from bottom:
                    dy = 0.15
                    bottom,height = bottom+dy,height-dy
                else :
                    # remove part from left side:
                    dx = 0.20
                    width = width-dx
                #endif
            #endif
            # data ax, need space for labels:
            dx,dy = 0.10,0.10
            left,width = left+dx,width-dx
            bottom,height = bottom+dy,height-dy
            # some extra space for title:
            if title is not None : height = height - 0.08
            # also some whitespace at right hand side:
            width = width - 0.05
            # combine:
            position = (left,bottom,width,height)
        #endif
        # add ax to figure:
        self.AddAxes( position )
        # add title:
        if title is not None : self.ax.set_title( title )
        # set data limits if necessary:
        self.Set_Domain( domain )
        
        # setup colormap:
        self.Set_ColorMap( **cmap )
        # no color normalization defined yet,
        # will be done by the first routine that shows values by colors:
        self.cnorm = None
        
        # enabled?
        if self.with_colorbar :
            # init colorbar; will be drawn at the first call to PColor or
            # other function that shows values by colors:
            self.Set_ColorBar( **cbar )
        #endif
        
    #enddef  # __init__
    
    # *

    def Set_ColorMap( self, colors=None,
                      color_under=None, color_over=None, color_bad=None,
                      ncolor=256 ) :

        """
        Define color map with N values.
        The rgb values in the colormap are interpolated between the specified colors.
        The 'under' and 'over' colors for lower and higher values,
        and 'bad' for masked data.
        The colormap is stored as the 'cmap' field.

        Arguments:
          colors        : either:
                            - list with color description: ['white','yellow','red']
                            - own definitions:
                              - 'wbrb' (white-blue-red-brown)
                              - 'brb'  (blue-red-brown)
                              - 'bwr'  (blue-white-red)  
                            - predefined colormaps : 'jet', ...
          color_under   : color descriptions for underflow
          color_over    : color descriptions for overflow
          color_bad     : color for masked values in numpy.ma (masked array)
          ncolor        : (int) number of slots in color table
        
        Attributes set by this routine:
          cmap          : colormap instance

        See also:
          matplotlib.colorbar  
        """

        # external:
        import numpy
        import matplotlib
        
        # default:
        if colors    is None : colors = 'wbrb'
        if color_bad is None : color_bad = '0.80'  # light gray
        
        # name or values?
        if type(colors) == str :
            # our favorite colors:
            if colors == 'wbrb' :
                colors = ['white','cyan','blue','green','yellow','orange','red','magenta','brown']
            elif colors == 'brbw' :
                colors = ['brown','magenta','red','orange','yellow','green','blue','cyan','white']
            elif colors == 'brb' :
                colors = ['cyan','blue','green','yellow','orange','red','magenta','brown']
            elif colors == 'bwr' :
                colors = ['blue','cyan','white','yellow','red']
            #endif
        #endif
        
        # list defined?
        if type(colors) == list :

            # defaults:
            if color_under is None : color_under = colors[0]
            if color_over  is None : color_over  = colors[-1]

            # number of colors:
            nslot = len(colors)

            # check ...
            if nslot < 2 :
                print( 'ERROR - at least 2 colors should be provided' )
                raise
            #endif

            # get red/green/blue arrays for extensions:
            red_under,green_under,blue_under = matplotlib.colors.colorConverter.to_rgb(color_under)
            red_over ,green_over ,blue_over  = matplotlib.colors.colorConverter.to_rgb(color_over)

            # initialise color dictionary:
            cdict = { 'red' : [], 'green' : [], 'blue' : [] }

            # intersection values in [0,1] :
            slots = numpy.arange(nslot)/float(nslot-1)

            # loop over colors:
            for islot in range(nslot) :

                # convert to rgb value:
                red,green,blue = matplotlib.colors.colorConverter.to_rgb(colors[islot])

                # current intersection:
                x = slots[islot]

                # fill colors:
                if islot == 0 :
                    red_y0,red_y1     = red_under,red
                    green_y0,green_y1 = green_under,green
                    blue_y0,blue_y1   = blue_under,blue
                elif islot == nslot-1 :
                    red_y0,red_y1     = red,red_over
                    green_y0,green_y1 = green,green_over
                    blue_y0,blue_y1   = blue,blue_over
                else :
                    red_y0,red_y1     = red,red
                    green_y0,green_y1 = green,green
                    blue_y0,blue_y1   = blue,blue
                #endif

                # add tupple:
                cdict['red'  ].append( (x,red_y0  ,red_y1  ) )
                cdict['green'].append( (x,green_y0,green_y1) )
                cdict['blue' ].append( (x,blue_y0 ,blue_y1 ) )

            #endfor

            #for key in cdict.keys() :
            #    print key
            #    for i in range(len(cdict[key])) : print cdict[key][i]
            ##endfor

            # create colormap:
            self.cmap = matplotlib.colors.LinearSegmentedColormap( 'mymap', cdict, N=ncolor )
            
        else :

            # get predefined colormap:
            self.cmap = matplotlib.pyplot.get_cmap( colors )
            
        #endif
            
        # boundaries:
        self.cmap.set_under(color_under)
        self.cmap.set_over(color_over)
        self.cmap.set_bad(color_bad)
        

    #enddef  # Set_ColorMap

    # *

    def Set_ColorNorm( self, values, 
                         vmin=None, vmax=None,
                         vbounds=None ) :
    
        """
        Setup a color normalization to map data values to an index in [0,1] .
        If an array 'vbounds' is specified, it should contain an increasing range
        of values that define the boundaries of interfals which have the same
        number of slots in the color map.
        Otherwise, the values are mapped linear from values in the interval 
        '[vmin,vmax]' to the slots in the color map.linear to [0,1] ;
        if these are not provided, the range will be the min/max of the
        suplied values.
    
        This routine sets the fields: cnorm vmin vmax
    
        Arguments:
          values          :  data to be displayed; use for min/max
          
        Optional keywords:
          vmin, vmax      :  data range
    
        See also:
          matplotlib.colors.normalize  
        """
    
        # external:
        import matplotlib
        
        # store value range:
        self.vmin = vmin
        self.vmax = vmax

        # bounds or not ?
        if vbounds is None :
    
            # set value range if not explicitly done:
            if self.vmin is None : self.vmin = values.min()
            if self.vmax is None : self.vmax = values.max()

            # trap constant value ...
            if self.vmin == self.vmax :
                # zero or values ?
                if self.vmin == 0.0 :
                    # dummy ...
                    self.vmin = 0.0
                    self.vmax = 1.0
                else :
                    # 10% around:
                    self.vmin = 0.9 * self.vmin
                    self.vmax = 1.1 * self.vmax
                #endif
            #endif
            
            # linear scaling:
            self.cnorm = matplotlib.colors.Normalize( vmin=self.vmin, vmax=self.vmax, clip=False )

        else :
    
            # set value range if not explicitly done:
            if self.vmin is None : self.vmin = min(vbounds)
            if self.vmax is None : self.vmax = max(vbounds)
            
            # number of colors:
            ncolor = self.cmap.N
        
            # number of intervals:
            ncell = len(vbounds) - 1
            # number of colors per cell:
            ncol = int(float(ncolor)/float(ncell))
            # init list of values in bounds:
            boundaries = []
            # add contributions for each interval:
            for icell in range(ncell) :
                # distribute the colors within the cell:
                for icol in range(ncol) :
                    # value within bounds:
                    boundaries.append( vbounds[icell] + icol*(vbounds[icell+1]-vbounds[icell])/float(ncol) )
                #endfor  # colors within cell:
            #endfor  # cells
            # add final boundary:
            boundaries.append( vbounds[-1] )

            # define normalization function to map from intervals defined
            # by high-res boundary list to indices in color map :
            self.cnorm = matplotlib.colors.BoundaryNorm( boundaries, ncolor, clip=False )
            ## use the fixed class defined in the top of this file:
            #self.cnorm = BoundaryNormX( boundaries, ncolor, clip=False )
            
            # testing ...
            #self.anorm = matplotlib.colors.Normalize( vmin=self.vmin, vmax=self.vmax, clip=False )

            
        #endif
    
    #enddef  # Set_ColorNorm
    
    # *
    
    def Set_ColorBar( self, ticks=None, ntick=None, label=None, tickprops=None, **kwargs ) :
    
        """
        Set colorbar properties.
        
        Keyword arguments:
          ticks     : tick locations (explicitly set)
          ntick     : number of tickmarks (distributed within vrange)
          label     : text below colorbar
          tickprops : properties for ticklabels, e.g. fontsize=20 etc ;
                      see "matplotlib.text.Text" for all supported properties
          
        Keyword arguments passed directly to 'matplotlib.pyplot.colorbar', for example:
          orientation : default 'horizontal', also 'vertical' possible
          extend      : add left or right arrows: 'neither' | 'both' | 'min' | 'max'
          drawedges   : lines between colors in colorbar ? True or False
          format      : used for labels, e.g. '%3.1f'
          
        """
        
        # store for use by Add_ColorBar :
        self.cbar_ticks     = ticks
        self.cbar_ntick     = ntick
        self.cbar_label     = label
        self.cbar_tickprops = tickprops
        self.cbar_kwargs    = kwargs
        # no colorbar added yet:
        self.cbar = None

    #enddef
    
    # *
    
    def Add_ColorBar( self, p ) :
    
        """
        Fill the colorbar for the mappable p.
        Uses settings stored by Set_ColorBar.
        
        Arguments:
          p         : mappable object, output of 'pcolor', 'scatter', ...
          
        """
        
        # not defined yet ?
        if self.with_colorbar and (self.cbar is None) :
        
            # modules:
            import matplotlib
            import matplotlib.pyplot
            import numpy

            # 'ticks' contains the location where colors should be labeled;
            # not explicitly defined ?
            if self.cbar_ticks is None :
                # number is defined ?
                if self.cbar_ntick is not None :
                    # get value range:
                    vmin,vmax = p.get_clim()
                    # linear from min to max:
                    self.cbar_ticks = vmin + numpy.arange(self.cbar_ntick+1)/float(self.cbar_ntick)*(vmax-vmin)
                #endif
            #endif

            # arguments:
            colorbar_kwargs = {}
            colorbar_kwargs.update( { 'mappable'    : p } )
            colorbar_kwargs.update( { 'orientation' : self.cbar_orientation } )
            colorbar_kwargs.update( { 'ticks'       : self.cbar_ticks } )
            colorbar_kwargs.update( self.cbar_kwargs )

            # add colorbar:
            if self.cbar_ax is None :
                # steele space from the data ax:
                self.cbar = matplotlib.pyplot.colorbar( ax=self.ax, **colorbar_kwargs )
            else :
                # plot it in the the speial colorbar ax:
                self.cbar = matplotlib.pyplot.colorbar( cax=self.cbar_ax, **colorbar_kwargs )
            #endif

            # add label ?
            if self.cbar_label is not None : self.cbar.set_label( self.cbar_label )

            # extra properites for ticklabels ?
            if self.cbar_tickprops is not None :
                # init list to collect text values from labels:
                labels = []
                # need to extract from x or y axis depending on orientation:
                if self.cbar_orientation == 'horizontal' :
                    # extract:
                    for ticklabel in self.cbar.ax.get_xticklabels() : labels.append( ticklabel.get_text() )
                    # reset with new properties:
                    self.cbar.ax.set_xticklabels( labels, **self.cbar_tickprops )
                else :
                    # extract:
                    for ticklabel in self.cbar.ax.get_yticklabels() : labels.append( ticklabel.get_text() )
                    # reset with new properties:
                    self.cbar.ax.set_yticklabels( labels, **self.cbar_tickprops )
                #endif
            #endif
            
        #endif  # no cbar yet

    #enddef  Add_ColorBar
    
    # *
    
    def PColor( self, xx, yy, cc, 
                  vmin=None, vmax=None, vbounds=None,
                  **kwargs ) :
    
        """
        Add pseudo color plot.
        If no grid keywords are provided, the domain specified in the
        initilization might be used; if no domain was provided, the
        axes are simply the dimension indices.
        
        Arguments:
          xx, yy     :  2D corner fields, shapes (ny+1,nx+1)
          cc         :  data values, shape (ny,nx)
        Optional grid definitions:
        Optional color mapping:
          vmin,vmax  :  data range for coloring, min and max values by default,
                        or the settings for Scatter if that was called before
          vbounds    :  instead of linear mapping of values to color map,
                        specify with this variable a list of boundary values,
                        each interval will take the same number of slots in
                        the color map; typically use this with explicit definition
                        of the ticks in the colorbar:
                            cbar=dict(ticks=vbounds)
        Other keywords are passed directly to 'matplotlib.pyplot.pcolormesh'.

        Return value: pcolormesh object
          
        """
        
        # modules:
        import numpy

#>>> not needed anymore ?        
#        # explicitly set value range, otherwise problems with no-data colors:
#        if vmin is None :
#            try :
#                vmin = numpy.nanmin(cc)
#            except :
#                vmin = numpy.nanmin(cc.data)
#            #endtry
#        #endif
#        if vmax is None :
#            try :
#                vmax = numpy.nanmax(cc)
#            except :
#                vmax = numpy.nanmax(cc.data)
#            #endtry
#        #endif
#<<<

        # Not needed anymore to work with color_bad for masked values
        ## extract 2d field ; first try something that is needed for masked arrays:
        #try :
        #    # replace the masked values values by NaN,
        #    # otherwise it is not colored correctly:
        #    cc = field.filled(fill_value=numpy.NaN)
        #except :
        #    # just copy:
        #    cc = field.copy()
        ##endtry
        ## copy:
        #cc = field.copy()
        
        # set mapping from data to color if not done yet:
        if self.cnorm is None :
            # set norm function:
            self.Set_ColorNorm( cc, vmin=vmin, vmax=vmax, vbounds=vbounds )
        #endif
        
        # pcolor keywords:
        pcolor_kwargs = {}
        pcolor_kwargs.update( { 'edgecolors' : 'none' } )
        pcolor_kwargs.update( { 'cmap' : self.cmap } )
        pcolor_kwargs.update( { 'norm' : self.cnorm } )
        pcolor_kwargs.update( kwargs )
        
        # add to map:
        p = self.ax.pcolormesh( xx, yy, cc, **pcolor_kwargs )
    
        ## freeze axes:
        #self.ax.set_autoscale_on(False)
        
        # add colorbar if not done yet:
        self.Add_ColorBar( p )

        # ok
        return p
        
    #enddef PColor
    
    # *
    
    def PColorRGG( self, rgg, values, vmin=None, vmax=None, vbounds=None, **kwargs ) :
    
        """
        Add pseudocolor fieds for values defined on reduced gaussian grid.
        
        Arguments:
            rgg     :  grid_rgg.ReducedGrid object
            values  :  1D list with values, length is rgg.npoint
            
        Optional color mapping (see PColor):
          vmin,vmax
          vbounds
                        
        Other keywords are passed directly to 'matplotlib.pyplot.pcolormesh'.
        """
        
        # modules:
        import numpy

        # set mapping from data to color if not done yet:
        if self.cnorm is None :
            # set norm function:
            self.Set_ColorNorm( values, vmin=vmin, vmax=vmax, vbounds=vbounds )
        #endif
        
        # pcolor keywords:
        pcolor_kwargs = {}
        pcolor_kwargs.update( { 'edgecolors' : 'none' } )
        pcolor_kwargs.update( { 'cmap' : self.cmap } )
        pcolor_kwargs.update( { 'norm' : self.cnorm } )
        pcolor_kwargs.update( kwargs )
        
        # loop over latitude bands:
        for j in range(rgg.nlat) :
            # number of longitude points in this band:
            nlon = rgg.band_nlon[j]
            # range of indices in 1D array with values for this band:
            i0 = rgg.band_i0[j]
            i1 = i0 + nlon
            # longitude corners:
            xx = numpy.hstack((rgg.band_lons_bnds[j][0,0],rgg.band_lons_bnds[j][:,1]))
            yy = rgg.band_lats_bnds[j,:]
            # values, convert to 2D array:
            cc = values[i0:i1].reshape((1,nlon))
            # add:
            p = self.ax.pcolormesh( xx, yy, cc, **pcolor_kwargs )
        #endfor

        # add colorbar if not done yet:
        self.Add_ColorBar( p )
        
    #enddef PColorRGG
    
    # *
    
    def Scatter( self, xx, yy, vmin=None, vmax=None, vbounds=None, **kwargs ) :
    
        """
        Plot cloud of scattered symbols.
        
        Arguments:
          xx, yy     :  location lists
          vmin,vmax  :  data range for coloring, min and max values by default,
                        or the settings for PColor if that was called before
          vbounds    :  value boundaries (alternative for vmin,vmax)

        Optional keywords passed to 'pyplot.scatter', for example:
          c          :  (list of) data value(s), mapped to colors
          s          :  (list of) symbol size(s)
          marker     :  marker symbol
          linewidths :  edge widths

        """
        
        # set mapping from data to color if not done yet:
        if (self.cnorm is None) and ('c' in kwargs.keys()) :
            # set norm function:
            self.Set_ColorNorm( kwargs['c'], vmin=vmin, vmax=vmax, vbounds=vbounds )
        #endif
        
        # scatter keywords:
        scatter_kwargs = {}
        scatter_kwargs.update( { 'cmap' : self.cmap } )
        scatter_kwargs.update( { 'norm' : self.cnorm } )
        scatter_kwargs.update( kwargs )
        
        # normal ax or map ?
        if self.bmp is None :
            # add to ax:
            p = self.ax.scatter( xx, yy, **scatter_kwargs )
        else :
            # add to map:
            p = self.bmp.scatter( xx, yy, **scatter_kwargs )
        #endif
    
        # add colorbar if not done yet:
        if 'c' in kwargs.keys() :
            self.Add_ColorBar( p )
        #endif
        
    #enddef Scatter
    
    # *
        
    def Density( self, x, y, bins=10, normed=False, **kwargs ) :
    
        """
        Add density plot (2D histogram).
        
        Arguments are passed to 'plt.hist2d' :
          x, y   : values
          bins   : number, array, or list of 2 arrays
        """
        
        # modules:
        import numpy

#        # set mapping from data to color if not done yet:
#        if self.cnorm is None :
#            # set norm function:
#            self.Set_ColorNorm( cc, vmin=vmin, vmax=vmax, vbounds=vbounds )
#        #endif
        
        # hist2d keywords:
        hist2d_kwargs = {}
        hist2d_kwargs.update( { 'bins'   : bins } )
        hist2d_kwargs.update( { 'cmap'   : self.cmap } )
        hist2d_kwargs.update( { 'normed' : normed } )
        hist2d_kwargs.update( kwargs )
        
        # add to map:
        counts,xedges,yedges,p = self.ax.hist2d( x, y, **hist2d_kwargs )
        
        # 1-1 line:
        #self.ax.plot( self.ax.get_xlim(), self.ax.get_ylim(), 'k--' )
        
        # add colorbar if not done yet:
        self.Add_ColorBar( p )
        
    #enddef Density

#endclass   # ColorbarFigure



########################################################################
###
### map figure with colorbar
###
########################################################################


class ColorbarMapFigure( ColorbarFigure ) :

    """
    Base class for figure with colorbar.
    """

    def __init__( self, figsize=None, 
                    position=None, bmp={}, 
                    cmap={}, cbar={} ) :

        """
        Setup figure with colorbar.
        
        Optional arguments passed to 'Figure' initialization:
            figsize
          
        Optional arguments to setup data axes:
            position=(left,bottom,width,height)
                position of axes; defaults is guesed
            bmp
                define geograpical map, passed to 'Setup_BaseMap'
        
        Optional arguments that define colorbar:
            cmap   : defines colormap, passed to 'Set_ColorMap'
            cbar   : defines colorbar, passed to 'Set_ColorBar' (False to skip)
            
        """

        # with colorbar?
        self.with_colorbar = True
        if cbar == False : self.with_colorbar = False
        
        # data axis:
        if position is None :
            # start with maximum ax position:
            left,bottom,width,height = 0.02,0.02,0.96,0.96
            # which direction ?
            if self.with_colorbar :
                if 'orientation' in cbar.keys() :
                    cbar_orientation = cbar['orientation']
                else :
                    cbar_orientation = 'horizontal'
                #endif
                # switch:
                if cbar_orientation == 'horizontal' :
                    # remove part from bottom:
                    dy = 0.20
                    bottom,height = bottom+dy,height-dy
                else :
                    # remove part from right side:
                    dx = 0.20
                    width = width-dx
                #endif
            #endif
            # setup ax for basemap; raster labels needed ?
            if 'raster' in bmp.keys() :
                # remove part from left and bottom:
                dx,dy = 0.10,0.10
                left,width = left+dx,width-dx
                bottom,height = bottom+dy,height-dy
            #endif
            # some extra space for title:
            if 'title' in bmp.keys() : height = height - 0.06
            # combine:
            position = (left,bottom,width,height)
        #endif
        
        # init figure:
        ColorbarFigure.__init__( self, figsize=figsize, 
                                   position=position, 
                                   cmap=cmap, cbar=cbar )
        
        # replace ax by map:
        self.ax0 = self.ax
        self.ax = Setup_BaseMap( self.ax0, **bmp )
        
    #enddef  # __init__
    
#endclass ColorbarMapFigure



########################################################################
###
### tools
###
########################################################################


def GetColorMap( colors=None,
                  color_under=None, color_over=None, color_bad=None,
                  ncolor=256 ) :

    """
    Define color map with N values.
    The rgb values in the colormap are interpolated between the specified colors.
    The 'under' and 'over' colors for lower and higher values,
    and 'bad' for masked data.
    The colormap is stored as the 'cmap' field.

    Arguments:
      colors        : either:
                        - list with color description: ['white','yellow','red']
                        - 'wbrb' (white-blue-red-brown)
                        - 'brb'  (blue-red-brown)
                        - 'bwr'  (blue-white-red)  
      color_under   : color descriptions for underflow
      color_over    : color descriptions for overflow
      color_bad     : color for masked values in numpy.ma (masked array)
      ncolor        : (int) number of slots in color table

    Return values:
      cmap          : colormap instance

    See also:
      matplotlib.colorbar  
    """

    # external:
    import numpy
    import matplotlib

    # default colors:
    if (colors is None) or (colors == 'wbrb') :
        colors = ['white','cyan','blue','green','yellow','orange','red','magenta','brown']
    elif colors == 'brb' :
        colors = ['cyan','blue','green','yellow','orange','red','magenta','brown']
    elif colors == 'bwr' :
        colors = ['blue','cyan','white','yellow','red']
    #endif

    # defaults:
    if color_under is None : color_under = colors[0]
    if color_over  is None : color_over  = colors[-1]
    if color_bad   is None : color_bad = '0.80'  # light gray

    # number of colors:
    nslot = len(colors)

    # check ...
    if nslot < 2 :
        print( 'ERROR - at least 2 colors should be provided' )
        raise
    #endif

    # get red/green/blue arrays for extensions:
    red_under,green_under,blue_under = matplotlib.colors.colorConverter.to_rgb(color_under)
    red_over ,green_over ,blue_over  = matplotlib.colors.colorConverter.to_rgb(color_over)

    # initialise color dictionary:
    cdict = { 'red' : [], 'green' : [], 'blue' : [] }

    # intersection values in [0,1] :
    slots = numpy.arange(nslot)/float(nslot-1)

    # loop over colors:
    for islot in range(nslot) :

        # convert to rgb value:
        red,green,blue = matplotlib.colors.colorConverter.to_rgb(colors[islot])

        # current intersection:
        x = slots[islot]

        # fill colors:
        if islot == 0 :
            red_y0,red_y1     = red_under,red
            green_y0,green_y1 = green_under,green
            blue_y0,blue_y1   = blue_under,blue
        elif islot == nslot-1 :
            red_y0,red_y1     = red,red_over
            green_y0,green_y1 = green,green_over
            blue_y0,blue_y1   = blue,blue_over
        else :
            red_y0,red_y1     = red,red
            green_y0,green_y1 = green,green
            blue_y0,blue_y1   = blue,blue
        #endif

        # add tupple:
        cdict['red'  ].append( (x,red_y0  ,red_y1  ) )
        cdict['green'].append( (x,green_y0,green_y1) )
        cdict['blue' ].append( (x,blue_y0 ,blue_y1 ) )

    #endfor

    #for key in cdict.keys() :
    #    print key
    #    for i in range(len(cdict[key])) : print cdict[key][i]
    ##endfor

    # create colormap:
    cmap = matplotlib.colors.LinearSegmentedColormap( 'mymap', cdict, N=ncolor )
    cmap.set_under(color_under)
    cmap.set_over(color_over)
    cmap.set_bad(color_bad)
    
    # ok
    return cmap

#enddef  # GetColorMap

# *

def GetGrid( shp, xx=None, yy=None, 
               x=None, y=None, xm=None, ym=None, xxm=None, yym=None,
               domain=None ) :

    """
    Return 2D grid arrays with corner points.
    
    Arguments:
      shp        :  shape (ny,nx) of values

    Optional grid definitions:
      xx, yy     :  2D corner fields, shapes (ny+1,nx+1) ;
                    if defined, these are the results
      xxm, yym   :  2D mid of cell fields, shapes (nx,ny) ;
      x, y       :  1D corner points, shapes (nx+1) and (ny+1); used to define 2D corner fields
      xm, ym     :  1D mid of cell, shapes (nx) and (ny); used to define 2D corners
   
    If none of these pairs is defined, a regular grid within the domain is assumed;
    if the domain is not specified, the default in (0,nx,0,ny).
    """    
    
    # modules:
    import numpy
    
    # extract:
    ny,nx = shp

    # grid not defined ?
    if (xx is None) or (yy is None) :
        # check:
        if (xx is not None) or (yy is not None) :
            print( 'ERROR - specify either both or none of xx,yy' )
            raise Exception
        #endif
        # 1D axes provided ?
        if (x is not None) or (y is not None) :
            # check ...
            if (x is None) or (y is None) :
                print( 'ERROR - specify both x and y' )
                raise Exception
            #endif
            # 2D corner fields:
            xx,yy = numpy.meshgrid(x,y)
        elif (xm is not None) or (ym is not None) :
            # check ...
            if (xm is None) or (ym is None) :
                print( 'ERROR - specify both xm and ym' )
                raise Exception
            #endif
            # create corner axes:
            x = mid2bounds( xm )
            y = mid2bounds( ym )
            # 2D corner fields:
            xx,yy = numpy.meshgrid(x,y)
        elif (xxm is not None) or (yym is not None) :
            # interpolate/extrapolate to corners:
            xx = mid2corners( xxm )
            yy = mid2corners( yym )
        else :
            # domain definition:
            if domain is None :
                left,right,bottom,top = 0,nx,0,ny
            else :
                left,right,bottom,top = domain
            #endif
            # regular axes:
            x = left   + numpy.arange(nx+1)*(right-left  )/float(nx)
            y = bottom + numpy.arange(ny+1)*(top  -bottom)/float(ny)
            # 2D corner fields:
            xx,yy = numpy.meshgrid(x,y)
        #endif
    #endif
    
    # ok
    return xx,yy
            
#enddef GetGrid


# *


def QuickPat( cc, xx=None, yy=None, x=None, y=None, xm=None, ym=None, xxm=None, yym=None, domain=None,
               vmin=None, vmax=None, vbounds=None, **kwargs ) :

    """
    
    Plot 2D data as colored pattern.
    
    Arguments passed to 'GetGrid' :
       xm, ym     :  1D mid of cell, shapes (nx) and (ny); used to define 2D corners
       x, y       :  1D corner points, shapes (nx+1) and (ny+1); used to define 2D corner fields
       xxm, yym   :  2D mid of cell, shapes (ny,nx)
       xx, yy     :  2D corner fields, shapes (ny+1,nx+1)
    
    Arguments passed to 'ColorbarFigure.PColor':
       cc         :  data values, shape (ny,nx)
       vmin,vmax  :  data range for coloring

    Other arguments are passed to 'ColorbarFigure' :
       figsize = (width,height)
       domain  = (left,right,bottom,top)
       cmap    = dict( 'colors'=['blue','white','red'], ncolor=9 )
       cbar    = dict( 'ntick'=9, extend='both', label='entity [units]' )
                 False   # skip colorbar
       
    Return value:    
      fig    :  Figure instance
    
    """
    
    # corner points:
    xx,yy = GetGrid( cc.shape, xx=xx, yy=yy, x=x, y=y, xm=xm, ym=ym, xxm=xxm, yym=yym, domain=domain )
    
    # set domain if not present yet:
    if domain is None : domain = [xx.min(),xx.max(),yy.min(),yy.max()]
    
    # setup figure:
    fig = ColorbarFigure( domain=domain, **kwargs )
    
    # add image:
    fig.PColor( xx, yy, cc, vmin=vmin, vmax=vmax, vbounds=vbounds )

    # ok
    return fig
    
#enddef # QuickPat


# *


def QuickScatter( xx, yy, domain=None, 
                   vmin=None, vmax=None, vbounds=None, 
                   scatter=dict(), **kwargs ) :

    """
    
    Plot cloud of scattered symbols.
        
    Arguments:
      xx, yy     :  location lists
      domain     :  [xmin,xmax,ymin,ymax]
      vmin,vmax  :  data range for coloring, min and max values by default,
                    or the settings for PColor if that was called before
      vbounds    :  value boundaries (alternative for vmin,vmax)
      
      scatter    :  keyword arguments passed to 'pyplot.scatter', for example:
        marker     :  marker symbol
        s          :  (list of) symbol size(s)
        c          :  (list of) data value(s), mapped to colors

      kwargs     :  keyword arguments passed to ColorbarFigure.
      
    Return value:    
      fig    :  Figure instance
    
    """
    
    # set domain if not present yet:
    if domain is None : domain = [xx.min(),xx.max(),yy.min(),yy.max()]
    
    # setup figure:
    fig = ColorbarFigure( domain=domain, **kwargs )
    
    # add image:
    fig.Scatter( xx, yy, vmin=vmin, vmax=vmax, vbounds=vbounds, **scatter )

    # ok
    return fig
    
#enddef # QuickScatter


# *


def QuickDens( x, y, bins=10, normed=False, **kwargs ) :

    """
    Density plot.
    
    Arguments passed to plt.hist2d :
      x, y       =  values
      bins       =  number, array, or 2-element list of arrays
      normed     =  frequency instead of numbers
    
    Other arguments are passed to 'ColorbarFigure' :
       figsize = (width,height)
       cmap    = dict( 'colors'=['blue','white','red'], ncolor=9 )
       cbar    = dict( 'ntick'=9, extend='both', label='entity [units]' )
       
    Return value:    
      fig    :  Figure instance
    
    """
    
    # setup figure:
    fig = ColorbarFigure( **kwargs )
    
    # add image:
    fig.Density( x, y, bins=bins, normed=normed )

    # ok
    return fig
    
#enddef # QuickDens


# *


def QuickMap( cc, xx=None, yy=None, x=None, y=None, xm=None, ym=None, xxm=None, yym=None, domain=None,
               vmin=None, vmax=None, vbounds=None, bmp=dict(), cbar=dict(), **kwargs ) :

    """
    
    Plot 2D data as colored pattern on map.
    
    
    Arguments passed to 'GetGrid' :
       xm, ym     :  1D mid of cell, shapes (nx) and (ny); used to define 2D corners
       x, y       :  1D corner points, shapes (nx+1) and (ny+1); used to define 2D corner fields
       xxm, yym   :  2D mid of cell, shapes (ny,nx)
       xx, yy     :  2D corner fields, shapes (ny+1,nx+1)
    Also a domain=(left,right,bottom,top) might be used to setup the map;
    could be part of the optional 'bmp' argument.
    
    Arguments passed to 'ColorbarFigure.PColor':
       cc         :  data values, shape (ny,nx)
       vmin,vmax  :  data range for coloring
       vbounds    :  value bounds for colorbar (defines interval with same color)

    Other arguments are passed to 'ColorbarMapFigure' :
       figsize = (width,height)
       bmp     = dict( resolution='c', countries=True, 
                         domain=(west,east,south,north), raster=[30.0,30.0],
                         title='displayed on top' )
       cmap    = dict( 'colors'=['blue','white','red'], ncolor=9 )
       cbar    = dict( 'ntick'=9, extend='both', label='entity [units]' )
                 False   # to skip colorbar
                 
    Return value:    
      fig    :  Figure instance
    
    """
    
    # domain defined as part of bmp ?
    if (domain is None) and ('domain' in bmp.keys()) : domain = bmp['domain']
    
    # corner points:
    xx,yy = GetGrid( cc.shape, xx=xx, yy=yy, x=x, y=y, xm=xm, ym=ym, xxm=xxm, yym=yym, domain=domain )
    
    # set domain if not present yet:
    if domain is None : domain = [xx.min(),xx.max(),yy.min(),yy.max()]
    
    # replace or add domain in bmp argument:
    bmp['domain'] = domain
    
    # auto extend?
    if 'extend' not in cbar.keys() :
        # extend below minimum?
        extend_min = (vmin is not None) and (cc.min() < vmin)
        # extend above maximum?
        extend_max = (vmax is not None) and (cc.max() > vmax)
        # set keyword:
        if extend_min and (not extend_max) :
            cbar['extend'] = 'min'
        elif extend_min and extend_max :
            cbar['extend'] = 'both'
        elif (not extend_min) and extend_max :
            cbar['extend'] = 'max'
        else :
            cbar['extend'] = 'neither'
        #endif
    #endif

    # setup figure:
    fig = ColorbarMapFigure( bmp=bmp, cbar=cbar, **kwargs )
    
    # add image:
    fig.PColor( xx, yy, cc, vmin=vmin, vmax=vmax, vbounds=vbounds )

    # ok
    return fig
    
#enddef



########################################################################
###
### test
###
########################################################################

if __name__ == "__main__" :

    # dummy map:
    nx,ny = 4,3
    xb = numpy.arange(nx+1)
    yb = numpy.arange(ny+1)
    xxb,yyb = numpy.meshgrid(xb,yb)
    cc = numpy.random.randn(ny,nx)

    # empty figure:
    #fig = Figure()
    
    # figure with global map:
    #fig = MapFigure()
    
    ## figure with colorbar:
    #fig = ColorbarFigure()
    #fig.PColor( xxb, yyb, cc )
    
    ## figure with colorbar:
    #fig = ColorbarMapFigure()
    #fig.PColor( xxb, yyb, cc )
    
    #QuickPat( cc )
    QuickMap( cc )

#endif


########################################################################
###
### end
###
########################################################################

