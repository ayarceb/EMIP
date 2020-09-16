
#-------------------------------------------------
# help
#-------------------------------------------------

"""
********************
``emip_satx`` module
********************

The :py:mod:`emip_satx` module provides classes
for processing satellite extract data.


Class hierchy
=============

The classes and are defined according to the following hierchy:


* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

    * :py:class:`.EmipTask`
    
      * :py:class:`.EmipSatxRegrid`
      * :py:class:`.EmipSatxCatalogue`


Classes
=======


"""


#-------------------------------------------------
# modules
#-------------------------------------------------

# tools:
import emip


#-------------------------------------------------
# OMI catalogue task
#-------------------------------------------------

class EmipSatxCatalogue( emip.EmipTask ) :

    """
    EMIP task to create catalogue of satellite figures 
    from extract files.

    .. figure:: figs/omi-catalogue.png
       :scale: 50 %
       :align: center
       :alt: OMI image catalogue

       *Example of image catalogue produced from converted OMI files.*

    In the settings, specify a time range for which images should be created::

      <rcbase>.timerange.start  :  2012-01-01 00:00
      <rcbase>.timerange.end    :  2012-12-31 23:59

    The images as well as an html index are written to a single directory
    specified with::

      <rcbase>.dir    :  /data/OMI-selection/catalgoue

    The location of the converted OMI files could be specified with
    time templates and a filename filter::

      ! converted OMI files, absolute path or relative to catalogue:
      <rcbase>.input.filenames        :  ../%Y/%m/OMI-Aura_NO2_*.nc

    Specify a list of variables to be plotted; usually this is the tropospheric
    vertical column density that is the main product, but also variables might be
    of interest::

      <rcbase>.vars                   :  vcd_trop

    Per variable the maximum value for the color bar could be specified;
    if not defined, the color bar is simply stretched to the maximum value
    present in the data:

      <rcbase>.var.vcd_trop.vmax      :  20.0

    Specify the domain of the map, projection is regular longitude/latitude:

      ! map domain (west east south north):
      <rcbase>.domain       :  -30 45 35 75

    Enable the following flag to re-create existing files,
    by default only non-existing files are created::

      <rcbase>.renew                  :  False

    When finished, the :py:mod:`catalogue` module is used create index pages.
    The url that should be loaded in a browser is shown::

      Point your browser to :
        file:///data/OMI-selection/catalgoue/index.html

    """
    
    def __init__( self, rcfile, rcbase='', env={}, indent='' ) :
    
        """
        Create catalogue of OMI figures.
        """
                        
        # modules:
        import os
        import datetime
        import fnmatch
        import numpy
        import matplotlib.pyplot as plt
        
        # tools:
        import satx
        import go
        import catalogue

        # init base object:
        emip.EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # time range:
        t1 = self.GetSetting( 'timerange.start', totype='datetime' )
        t2 = self.GetSetting( 'timerange.end'  , totype='datetime' )
        
        # target directory:
        workdir = self.GetSetting( 'dir' )

        # info ..
        self.logger.info( indent+'work directory:' )
        self.logger.info( indent+'  %s' % workdir )
        # create?
        if not os.path.isdir(workdir) : os.makedirs( workdir )
        # store current:
        owd = os.getcwd()
        # change to:
        os.chdir( workdir )
        
        # template for input files:
        template = self.GetSetting( 'input.filenames' )
        # split:
        inputdir = os.path.dirname( template )
        filename = os.path.basename( template )

        # init list with input filenames:
        filenames = []
        # content of input directories is listed,
        # name of directory changes in time:
        indir_prev = None
        # loop over time range:
        t = t1
        while t <= t2 :
            # resolve time templates:
            indir = t.strftime(inputdir)
            # new?
            if indir != indir_prev :
                # check ...
                if not os.path.isdir( indir ) :
                    self.logger.error( 'input directory not found:' )
                    self.logger.error( '  %s' % indir )
                    raise Exception
                #endif
                # list:
                fnames = os.listdir( indir )
                # sort (in place)
                fnames.sort()
                # add full paths:
                for fname in fnames :
                    # filter ...
                    if not fnmatch.fnmatch( fname, filename ) : continue
                    # store:
                    filenames.append( os.path.join(indir,fname) )
                #endfor
                # reset:
                indir_prev = indir
            #endif
            # next:
            t = t + datetime.timedelta(1)
        #endwhile
        
        # variables to be plotted:
        varnames = self.GetSetting( 'vars' ).split()

        # renew existing figures?
        renew = self.GetSetting( 'renew', totype='bool' )
        
        # plot domain:
        domain = list(map( float, self.GetSetting('domain').split() ))
        
        # mid of domain, used for local time:
        domain_lon_mid = 0.5*(domain[0]+domain[1])
        # time zone shift based on bands of 360/24=15 deg
        tzone_dhour = int(numpy.round( domain_lon_mid / 15.0 ))
        # increment:
        tzone_dt = datetime.timedelta(0,3600*tzone_dhour)
        # annote:
        if tzone_dhour < 0 :
            tzone_label = 'UTC%i' % tzone_dhour
        elif tzone_dhour > 0 :
            tzone_label = 'UTC+%i' % tzone_dhour
        else :
            tzone_label = 'UTC'
        #endif
        
        # figure size:
        figsize = eval( self.GetSetting('figsize') )
        
        # collections for index:
        dates = []
        orbits = []
        
        # info ..
        self.logger.info( indent+'loop over orbit files ...' )
        # loop:
        for filename in filenames :
            # info ..
            self.logger.info( indent+'  %s ...' % filename )
            
            # basenames:
            #    OMI-Aura_NO2_20120101_o39701
            #    OMI-Aura_NO2_20120101_1200
            bname = os.path.basename(filename).replace('.nc','')
            # split:
            prefix,comp,date,orbit = bname.rsplit('_',3)
            
            # filter on time ..
            orb_t = datetime.datetime.strptime( date, '%Y%m%d' )
            if (orb_t < t1) or (orb_t > t2) :
                self.logger.info( indent+'    outside time range, skip ..' )
                continue
            #endif
            
            # store:
            if date not in dates : dates.append( date )
            if orbit not in orbits : orbits.append( orbit )
            
            # no orbit file read yet ..
            orb = None
            
            # loop over variables ...
            for varname in varnames :
            
                # target file:
                figfile = '%s__%s.png' % (bname,varname)
                # renew?
                if (not os.path.isfile(figfile)) or renew :
                    # info ..
                    self.logger.info( indent+'    create %s ...' % figfile )
                    
                    # read?
                    if orb is None :
                        # info ...
                        self.logger.info( indent+'      read ...' )
                        # init storage and read:
                        orb = satx.SatXFile( filename=filename )
                    #endif
                    
                    # time average:
                    taver = orb.GetTimeAverage()
                    # display as local standard time:
                    taver_lst = taver + tzone_dt
                    # annote:
                    title = '%s (%s) %s' % (taver_lst.strftime('%Y-%m-%d %H:%M'),tzone_label,orbit)
                    
                    # style:
                    vkey = 'var.%s' % varname
                    vunits = self.GetSetting( vkey+'.units' , default='None' )
                    vmin   = eval( self.GetSetting( vkey+'.vmin'  , default='None' ) )
                    vmax   = eval( self.GetSetting( vkey+'.vmax'  , default='None' ) )
                    colors = eval( self.GetSetting( vkey+'.colors', default='None' ) )
                    color  = self.GetSetting( vkey+'.color', default='#1f77b4' )
                    
                    # dimensions:
                    dims = orb.GetDims( varname )
                    # switch:
                    #~ maps:
                    if dims == ('pixel',) :

                        # extract corner grids and values:
                        xx,yy,values,units = orb.GetTrack( varname )
                    
                        # convert if necessary:
                        values,units = self._ConvertUnits( values, units, vunits )

                        # annote:
                        label = '%s [%s]' % (varname,units)

                        # create map figure:
                        fig = go.plot2.QuickMap( values, xx=xx, yy=yy, vmin=vmin, vmax=vmax,
                                                  cmap=dict(colors=colors),
                                                  bmp=dict(resolution='i',countries=True,domain=domain,title=title),
                                                  cbar=dict(label=label),
                                                  figsize=figsize )
                        # save:
                        fig.Export( figfile )

                        # close:
                        fig.Close()
                        
                    #~ profiles:
                    elif dims in [('pixel','layer'), ('pixel','layer_interface')] :
                    
                        # extract:
                        values = orb.var[varname]['data']
                        units  = orb.var[varname]['attrs']['units']
                    
                        # convert if necessary:
                        values,units = self._ConvertUnits( values, units, vunits )
                        # annote:
                        label = '%s [%s]' % (varname,units)
                        
                        # stats:
                        mu    = values.mean(axis=0)
                        sigma = values.std (axis=0)
                        mins  = values.min (axis=0)
                        maxs  = values.max (axis=0)
                        
                        # levels:
                        nlev = values.shape[1]
                        levels = numpy.arange(1,nlev+1)
                        ylim = [0.5,nlev+0.5]
                        
                        # new figure:
                        fig = plt.figure()
                        ax = fig.add_axes([0.1,0.1,0.8,0.8])
                        # plot mean:
                        ax.plot( mu, levels, linestyle='-' , color=color )
                        # mean +/- sigma
                        ax.plot( numpy.maximum(mu-sigma,mins), levels, linestyle='--', color=color )
                        ax.plot( numpy.minimum(mu+sigma,maxs), levels, linestyle='--', color=color )
                        # extrema:
                        ax.plot( mins, levels, linestyle='None', marker='.', color=color )
                        ax.plot( maxs, levels, linestyle='None', marker='.', color=color )
                        
                        # zero axis?
                        xlim = ax.get_xlim()
                        if (xlim[0] < 0.0) and (xlim[1] > 0.0) :
                            ax.plot( [0,0], ylim, 'k-' )
                        #endif

                        # x-axis:
                        ax.set_xlim([vmin,vmax])
                        ax.set_xlabel(label)
                        
                        # y-axis
                        ax.set_ylim(ylim)
                        ax.set_ylabel('level')
                        
                        # save:
                        fig.savefig( figfile )

                        # close:
                        plt.close(fig)

                #endif
                
            #endfor # variables
            
            ## testing ...
            #break
            
        #endfor # orbit files
        
        # create index:
        cata = catalogue.CatalogueCreator( '%s_%s_<date>_<orbit>__<varname>.png' % (prefix,comp), 
                                             title='OMI catalgoue',
                                             logger=self.logger, verbose=False )
        # define levels:
        cata.AddLevel( 'date'   , dates   , 'ul', True )
        cata.AddLevel( 'orbit'  , orbits  , 'tr', False )
        cata.AddLevel( 'varname', varnames, 'td', False )
        # html image template
        cata.AddTemplate( '*.png' , '<img src="%(filename)s" width="400" border="0">' )
        # create:
        cata.Create( basename='index', silent=False )
        
        # back:
        os.chdir( owd )
    
    #enddef __init__
    
    # *
    
    def _ConvertUnits( self, values, units, vunits ) :
    
        """
        Convert from 'units' to 'vunits' if necessary.
        """

        # no target units defined?
        if (vunits is None) or (vunits == 'None') :

            # unchanged:
            vunits = units
        
        # different units?
        elif units != vunits :

            # which conversion?
            conversion = '%s -> %s' % (units,vunits)
            # switch:
            if conversion == 'Pa -> hPa' :
                values = values/1e2
            else :
                self.logger.error( 'unsupported conversion "%s"' % conversion )
                raise Exception
            #endif

        #endif
        
        # ok
        return values,vunits
        
    #enddef _ConvertUnits

#endclass EmipSatxCatalogue


#-------------------------------------------------
# SatX regrid task
#-------------------------------------------------

class EmipSatxRegrid( emip.EmipTask ) :

    """
    EMIP task to regrid converted OMI files to a regular grid.
    Input files are expected to be produced by the
    :py:class:`.EmipSatxConvert` class.
    For each input file, an output file with similar format is created
    which has 'pixels' with a footprint equal to a grid cell.

    The remapping is done by distributing the footprint polygon over the grid cells.
    Each grid cell is therefore filled with a weighed sum of contributions from pixels
    that (partly) overlap the cell:
    
    .. math::
        y ~=~ \sum\limits_{i=1}^{npix} w_i\ x_i
        
    The weights are relative to the area covered by a pixel, thus the more area of a cell
    is covered by a pixel the more weight that pixel has.
    
    Variables that start with '``sigma_``' are assumed to be error estimates.
    The errors are assumed to be uncorrrelated between pixels, and therefore
    the combined error could be computed as weighted sum over variances:
    
    .. math::
        \sigma_y ~=~ \sqrt{ \sum\limits_{i=1}^{npix} w_i\ \sigma_{x,i}^2 }
    
    The configuration is read from the '``rcfile``',
    the keywords start with the '``rcbase``' which is for example '``emip.omi.regrid``'.
    
    The first setting that is read is the time range over which 
    files should be converted::

      <rcbase>.timerange.start        :  2012-01-01 00:00
      <rcbase>.timerange.end          :  2012-12-31 23:59

    A loop over this time range will be performed with steps per day.
    For each day, the content of an input directory is scanned for files.
    Specify the directory using time templates::

      <rcbase>.input.dir       :  /data/OMI-selection/%Y/%m

    Filenames are specified by a filter for the date and orbit number::

      <rcbase>.input.filenames :  OMI-Aura_NO2_%.nc
      
    The output directory is specfied including time templates;
    output files will have the same name as the input file::

      <rcbase>.output.dir       :  /data/OMI-gridded/%Y/%m

    The target filename will have a timestamp instead of an orbit number,
    which is easier for the assimilation system to search if any observations
    are present. The time stamp that is used will be the whole hour that
    is most nearby the average of pixel times on the track::
    
      <rcbase>.output.filename    :  OMI-Aura_NO2_%Y%m%d_%H%M.nc

    The following flag controls whether existing output files are renewed,
    otherwise it is assumed that the target file is present already and
    regridding is simply skipped::

      <rcbase>.output.renew       :  False
    
    If existing files should not be overwritten, enable the following flag::
    
      <rcbase>.output.overwrite   :  True
      
    The grid onto which the OMI data should be regridded should be defined
    by a lon/lat domain for the edges and a resolution::

      <rcbase>.grid.west      :  -30.0
      <rcbase>.grid.east      :   40.0
      <rcbase>.grid.south     :   30.0
      <rcbase>.grid.north     :   75.0

      <rcbase>.grid.dlon      :  0.2
      <rcbase>.grid.dlat      :  0.1
    
    Since computing the source pixels and weights is expensive, 
    all info on the mapping could be saved to (intermediate) files.
    Next time the regridding is run, for example to add extra variables,
    the weights could be loaded directly instead of recomputing them.
    The following flags control whether the weights sould be loaded,
    if they should be stored when (re)computed, and where to store them::

        <rcbase>.mapping.load     :  True
        <rcbase>.mapping.store    :  True
        <rcbase>.mapping.dir      :  /data/OMI-gridded__mapping/%Y/%m

    """
    
    def __init__( self, rcfile, rcbase='', env={}, indent='' ) :
    
        """
        Regrid satellite data.
        """
                        
        # modules:
        import os
        import datetime
        import fnmatch
        import numpy
        import netCDF4
        
        # tools:
        import satx
        import grid
        import go
        
        ## testing ..
        #import importlib
        #importlib.reload(grid.cg)
        #importlib.reload(go)

        # init base object:
        emip.EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # info ..
        self.logger.info( indent+'' )
        self.logger.info( indent+' ** Regrid satellite files **' )
        self.logger.info( indent+'' )

        # time range:
        t1 = self.GetSetting( 'timerange.start', totype='datetime' )
        t2 = self.GetSetting( 'timerange.end'  , totype='datetime' )
        
        # target domain:
        west  = self.GetSetting( 'grid.west' , totype='float' )
        east  = self.GetSetting( 'grid.east' , totype='float' )
        south = self.GetSetting( 'grid.south', totype='float' )
        north = self.GetSetting( 'grid.north', totype='float' )
        # target resolution:
        dlon = self.GetSetting( 'grid.dlon', totype='float' )
        dlat = self.GetSetting( 'grid.dlat', totype='float' )
        # shape:
        nlon = int( (east -west )/dlon )
        nlat = int( (north-south)/dlat )
        # info ...
        self.logger.info( indent+'define grid ...' )
        # define grid, eventually pre-compute polygons:
        cg = grid.cg.CarthesianGrid(  west=west , dlon=dlon, nlon=nlon, 
                                     south=south, dlat=dlat, nlat=nlat,
                                     polygons=False )
        # corner arrays:
        cxx,cyy = cg.GetCorners()
        
        # template for input files:
        input_dir       = self.GetSetting( 'input.dir' )
        input_filenames = self.GetSetting( 'input.filenames' )
        
        # template for output directory:
        output_dir      = self.GetSetting( 'output.dir' )
        # target filanme:
        output_filename = self.GetSetting( 'output.filename' )
        
        # renew target files?
        output_renew = self.GetSetting( 'output.renew', totype='bool' )
        # overwrite existing files?
        output_overwrite = self.GetSetting( 'output.overwrite', totype='bool' )

        # store/restore mapping?
        mapping_store = self.GetSetting( 'mapping.store', totype='bool' )
        mapping_load  = self.GetSetting( 'mapping.load', totype='bool' )
        # template for output directory:
        mapping_dir   = self.GetSetting( 'mapping.dir' )
        
        # scanned input directories:
        indirs = []
        
        # info ...
        self.logger.info( indent+'time loop ...' )
        # loop over time range:
        t = t1
        dt = datetime.timedelta(1)
        while t <= t2 :
        
            # resolve time templates:
            indir  = t.strftime(input_dir)
            outdir = t.strftime(output_dir)
            mapdir = t.strftime(mapping_dir)
            
            # already scanned?
            if indir in indirs :
                # next:
                t = t + dt
                continue
            #endif
            # store current:
            indirs.append( indir )

            # info ...
            self.logger.info( indent+'  scan %s ...' % indir )
            # list:
            fnames = os.listdir( indir )
            # sort (in place)
            fnames.sort()
            
            # loop over input files:
            for fname in fnames :
                # filter ...
                if not fnmatch.fnmatch( fname, input_filenames ) : continue
                
                # info ...
                self.logger.info( indent+'  %s ...' % fname )

                # info ...
                self.logger.info( indent+'    read time values...' )
                # full path:
                filename = os.path.join(indir,fname)
                # init storage and read:
                orb = satx.SatXFile( filename=filename, varnames=['time'] )
                # average time:
                taver = orb.GetTimeAverage()
                # round to nearby hour:
                hour = numpy.round( taver.hour + taver.minute/60.0 )
                tstamp = taver.replace(hour=0,minute=0,second=0) + datetime.timedelta(0,hour*3600)

                # target file:
                outfile = os.path.join( outdir, tstamp.strftime(output_filename) )
                
                # create?
                if (not os.path.isfile(outfile)) or output_renew :

                    # info ...
                    self.logger.info( indent+'    create %s ...' % os.path.basename(outfile) )
                
                    # check ...
                    if os.path.isfile(outfile) and (not output_overwrite) :
                        self.logger.error( 'output file already present:' )
                        self.logger.error( '  %s' % outfile )
                        raise Exception
                    #endif

                    # * read

                    # info ...
                    self.logger.info( indent+'      read all ...' )
                    # full path:
                    filename = os.path.join(indir,fname)
                    # init storage and read:
                    orb = satx.SatXFile( filename=filename )

                    # * mapping

                    # file with mapping:
                    mapping_file = os.path.join( mapdir, fname.replace('.nc','_mapping.nc') )

                    # load mapping?
                    if mapping_load and os.path.isfile(mapping_file) :

                        # info ...
                        self.logger.info( indent+'      read mapping ...' )
                        # read:
                        mapper = go.mapping.Mapping_1D_to_2D( mapping_file )

                    else :

                        # init mapping object:
                        mapper = go.mapping.Mapping_1D_to_2D( shp=(nlat,nlon) )

                        # inquuire:
                        npix = orb.GetDimension( 'pixel' )
                        
                        # info ...
                        self.logger.info( indent+'    process %i pixels ...' % npix )                        
                        # loop over pixels
                        for ipix in range(npix) :
                        
                            # testing ...
                            #if ipix < 543 : continue
                            #if ipix < 546 : continue

                            ## info ...
                            #self.logger.info( indent+'      pixel %i / %i ...' % (ipix+1,npix) )

                            # get polygon of footprint:
                            pgpix = orb.GetPixelFootprint( ipix )

                            # indices of covered grid cells and covered area:
                            pgpix_area,jj,ii,aa,units = cg.PolygonCoverage( pgpix, debug=False )
                            

                            ## info ..
                            #self.logger.info( indent+'        assigned to %i cells ...' % len(aa) )
                            # add contibution to mapper:
                            mapper.Append( ipix, pgpix_area, jj, ii, aa, units )
                            
                            ## testing ..
                            #raise Exception

                            ## testing ...
                            #if ipix == 10 : break
                            #if ipix == 559 : raise Exception
                        #endfor

                        # store?
                        if mapping_store :
                            # attributes:
                            attrs = {}
                            attrs['source_file'] = filename
                            # write:
                            mapper.Write( mapping_file, attrs=attrs )
                        #endif

                    #endif # load mapping or calculate new

                    # * apply mapping
                    
                    # any source values to be mapped?
                    if len(mapper) > 0 :
                        
                        # info ...
                        self.logger.info( indent+'    create gridded file ...' )

                        # init storage for gridded orbit:
                        orbg = satx.SatXFile()

                        # add target grid as 'track':
                        orbg.AppendTrack( 'track_longitude', cg.xxm,
                                            ('track_image','track_pixel'),
                                            attrs={ 'units' : 'degrees_east', 'long_name' : 'longitudes' } )
                        orbg.AppendTrack( 'track_latitude', cg.yym, 
                                            ('track_image','track_pixel'),
                                            attrs={ 'units' : 'degrees_north', 'long_name' : 'latitudes' } )
                        orbg.AppendTrack( 'track_corner_longitudes', cxx,
                                            ('track_image','track_pixel','corner'),
                                            attrs={ 'units' : 'degrees_east', 'long_name' : 'corner longitudes' } )
                        orbg.AppendTrack( 'track_corner_latitudes', cyy, 
                                            ('track_image','track_pixel','corner'),
                                            attrs={ 'units' : 'degrees_north', 'long_name' : 'corner latitudes' } )

                        # index arrays for target:
                        jj,ii = mapper.GetIndices()
                        # add pixel coordinates:
                        orbg.AppendTrack( 'longitude', cg.xxm[jj,ii], ('pixel',),
                                            attrs={ 'units' : 'degrees_east', 'long_name' : 'longitudes' } )
                        orbg.AppendTrack( 'latitude', cg.yym[jj,ii], ('pixel',),
                                            attrs={ 'units' : 'degrees_north', 'long_name' : 'latitudes' } )
                        orbg.AppendTrack( 'corner_longitudes', cxx[jj,ii,:], ('pixel','corner'),
                                            attrs={ 'units' : 'degrees_east', 'long_name' : 'corner longitudes' } )
                        orbg.AppendTrack( 'corner_latitudes', cyy[jj,ii,:], ('pixel','corner'),
                                            attrs={ 'units' : 'degrees_north', 'long_name' : 'corner latitudes' } )

                        # info for compression coordinate:
                        vname = 'pixel'
                        indx,attrs = mapper.GetCompression( vname, ('track_image','track_pixel') )
                        # extend
                        # add:
                        orbg.AddVariable( vname, indx, (vname,), attrs )

                        # loop over variables:
                        for vname in orb.var.keys() :

                            # some already done ...
                            if vname in ['longitude','latitude'] : continue
                            if vname in ['corner_longitudes','corner_latitudes'] : continue
                            if vname in ['track_longitude','track_latitude'] : continue
                            if vname in ['track_corner_longitudes','track_corner_latitudes'] : continue
                            if vname in ['pixel'] : continue

                            # info ..
                            self.logger.info( indent+'    regrid "%s" ..' % vname )
			  
                            # current:
                            var = orb.var[vname]

                            # convert?
                            if vname == 'time' :

                                # convert:
                                values = netCDF4.date2num( var['data'], var['attrs']['units'] )
                                # apply mapping:
                                values = mapper.Apply( values )
                                # convert:
                                values = netCDF4.num2date( values, var['attrs']['units'] )

                            # errors:
                            elif vname.startswith('sigma_') :

                                # apply mapping to variance:
                                values = numpy.sqrt( mapper.Apply( var['data']**2 ) )

                            # pixel variables:
                            elif var['dims'] in [ ('pixel',), ('pixel','layer'), ('pixel','layer_interface') ] :

                                # apply mapping:
                                values = mapper.Apply( var['data'] )

                            # level definitions:
                            elif var['dims'] in [ ('layer',), ('layer_interface',) ] :

                                # keep:
                                values = var['data']

                            else :

                                # unknown ..
                                self.logger.error( 'could not remap "%s" with dims "%s"' % (vname,var['dims']) )
                                raise Exception

                            #endif

                            # add variable:
                            orbg.AddVariable( vname, values, var['dims'], var['attrs']  )

                        #endfor # variables

                        # global attributes of input:

                        attrs = orb.attrs
                        # adhoc files had copies of all previously processed histories ...
                        attrs['history'] = attrs['history'].split('\n')[0]

                        # write, copy input attributes and extend history:
                        orbg.Write( outfile, attrs=attrs, history=['averaged %s to regular grid' % fname] )

                    else :
                                           
                        # info ...
                        self.logger.info( indent+'    no source values covering target grid ...' )
                        
                    #endif  # source values to be mapped
                    
                #endif  # create gridded file

                ## testing ...
                #break

            #endfor # input files
            
            ## testing ...
            #break

            # next:
            t = t + dt

        #endwhile  # time loop
        
    #enddef __init__
    
#endclass EmipSatxRegrid



#-------------------------------------------------
# end
#-------------------------------------------------




