
#-------------------------------------------------
# help
#-------------------------------------------------

"""
***********************
``emip_tropomi`` module
***********************

The :py:mod:`emip_tropomi` module provides classes
for processing TROPOMI data.


Class hierchy
=============

The classes and are defined according to the following hierchy:


* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

    * :py:class:`.EmipTask`

      * :py:class:`.EmipTropomiConvert`

  * :py:class:`.SatXFile`

    * :py:class:`.TropomiExtract`


Classes
=======


"""


#-------------------------------------------------
# modules
#-------------------------------------------------

# tools:
import emip


#-------------------------------------------------
# OMI extract file
#-------------------------------------------------

# modules:
import satx

# testing ...
import importlib
importlib.reload(satx)


class TropomiExtract( satx.SatXFile ) :

    """
    Storage for converted OMI file.
    """
    
    def AddSelection( self, tropdata, selected, rcf, rckey, orbnr, indent='' ) :
    
        """
        Add selected OMI to satellite extract file.
        
        Arguments:
        
        * ``tropdata``    :  :py:class:`.OMI_HDF5_File` object
        * ``selected``   :  boolean array as provided by :py:meth:`OMI_HDF5_File.SelectPixels` method
        * ``rcf``        :  :py:class:`.RcFile` instance with settings
        * ``rckey``      :  base of rcfile keys
        * ``orbnr``      :  orbit number to identify original track
            
        The first setting that is read is a list with variable names to be 
        created in the target file::

          <rcbase>.output.vars    :  longitude corner_longitudes \
                                     latitude corner_latitudes \
                                     vcd_trop  ...

        For each variable settings should be specified that describe 
        how to obtain the values and the target units (for automatic conversion if possible).

        For most variables it is sufficient to provide only the name of the original
        variable from which the data should be read::

          <rcbase>.output.var.longitude.from    :   Geolocation Fields/Longitude
          <rcbase>.output.var.longitude.units   :   degrees_east

        For some variables some special processing needs to be done.
        For these variables a key '``special``' is used which will enable the 
        correct code for conversion. For example, the following setting will
        ensure that the variable '``image_number``' will be filled with the scan
        number within the track::

          <rcbase>.output.var.image_number.special            :   scan_number
          <rcbase>.output.var.image_number.units              :   1

        If new variables require special processing, just insert a new '``special``' keyword
        and wait for the conversion class to complain about an unsupported value.

        """
        
        # modules:
        import os
        import numpy
        import netCDF4
        
        # info ..
        self.logger.info( indent+'add TROPOMI pixels to extract file ...' )
        
        # check ...
        if tropdata.ntime != 1 :
            self.logger.error( 'only single time record supported, found %i' % tropdata.ntime )
            self.logger.error( 'input file: %s' % tropdata.filename )
            raise Exception
        #endif
        # time record:
        itime = 0
        
        # copy dimension:
        #self.nlayer = tropdata.nlayer
        nlayer = tropdata.nlayer
        
        # selection indices:
        iit,iip = numpy.where( selected[itime,:,:] )
        # count:
        npix = len(iit)
        
        # target variables:
        varnames = rcf.get( rckey+'.output.vars' ).split()
        
        # loop over target variables:
        for varname in varnames :
        
            # info ...
            self.logger.info( indent+'  add variable "%s" ...' % varname )
        
            # target units:
            vunits = rcf.get( '%s.output.var.%s.units' % (rckey,varname) )
            
            # specials ?
            special = rcf.get( '%s.output.var.%s.special' % (rckey,varname), default='None' )
            # switch:
            if special == 'None' :
            
                # name of original field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # get variable dictionairy:
                #   odat['data'][:]
                #   odat['attrs'] = { 'units' : '1', ..}
                odat = tropdata.GetVar( opath )

                # extact values:
                values = odat['data']
                
                # check ...
                if 'units' not in odat['attrs'].keys() :
                    self.logger.error( 'no units in "%s" attributes' % opath )
                    raise Exception
                #endif
                # copy:
                ounits    = odat['attrs']['units']
                
                # copy:
                olongname = odat['attrs']['long_name']
                
                # dimension names:
                odimnames = odat['dimnames']
                
                # copy units?
                if vunits == '__native__' : vunits = ounits
                
                # conversion ?
                if ounits != vunits :
                    # required conversion:
                    conversion = '%s -> %s' % (ounits,vunits)
                    # switch:
                    if conversion in ['mol m-2 -> 1e15 mlc/cm2'] :
                        # check ...
                        aname = 'multiplication_factor_to_convert_to_molecules_percm2'
                        if aname not in odat['attrs'].keys() :
                            self.logger.error( 'attribute "%s" not found for conversion "%s"' % (aname,conversion) )
                            self.logger.error( 'variable : %s' % opath )
                            self.logger.error( 'filename : %s' % odat.filename )
                            raise Exception
                        #endif
                        # extract:
                        factor = odat['attrs'][aname]
                        # apply:
                        values = values * factor * 1.0e-15
                    #if conversion in ['deg -> degrees_east','deg -> degrees_north',
                    #                  'Degrees -> degrees_east', 'Degrees -> degrees_north',
                    #                  'NoUnits -> 1'] :
                    #    # no conversion needed
                    #    pass
                    #elif conversion in ["molecules cm^-2 -> 1e15 cm**-2",
                    #                    "molecules/cm^2 -> 1e15 cm**-2" ] :
                    #    # apply factor:
                    #    values = values * 1.0e-15
                    #elif conversion in ["hPa -> Pa"] :
                    #    # apply factor:
                    #    values = values * 1.0e2
                    #elif conversion in ["NoUnits -> %"] :
                    #    # assume original is in [0,1], convert to % :
                    #    values = values * 100.0
                    else :
                        self.logger.error( 'unsupported conversion "%s" requested ...' % conversion )
                        raise Exception
                    #endif
                #endif  # units differ
                
                # extract 1D array of pixels from 2D track:
                if odimnames in [('time','scanline','ground_pixel')] :

                    values = values[itime,iit,iip]
                    dimnames = ('pixel',)

                elif odimnames in [('time','scanline','ground_pixel','corner')] :

                    values = values[itime,iit,iip,:]
                    dimnames = ('pixel','corner')

                else :

                    # probably no pixels but layer definition etc,
                    # just copy dimension names and keep values:
                    dimnames = odimnames
                    
                #endif
                
                # attributes:
                attrs = {}
                attrs['units'] = vunits
                attrs['long_name'] = olongname

                # pixels?
                if dimnames[0] == 'pixel' :
                    # extend variable with latest pixels:
                    self.AppendPixels( varname, values, dimnames, attrs )
                else :
                    # new variable (or overwrite existing):
                    self.AddVariable( varname, values, dimnames, attrs )
                #endif

            # time per pixel
            elif special == 'time' :
            
                # reference time variable:
                tvar = tropdata.GetVar( 'PRODUCT/time' )
                tref = netCDF4.num2date( tvar['data'][itime], tvar['attrs']['units'] )
            
                # name of time delta field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # get variable dictionairy:
                #   odat['data'][:]
                #   odat['attrs'] = { 'units' : '1', ..}
                odat = tropdata.GetVar( opath )
                
                # step units:
                dtunits = odat['attrs']['units']
                # convert:
                if dtunits.startswith('milliseconds') :
                    # extract per pixel, convert to seconds:
                    values = odat['data'][itime,iit] * 1.0e-3
                else :
                    self.logger.error( 'unsupported "%s" units "%s"' % (opath,dtunits) )
                    raise Exception
                #endif
                
                # unit description:
                tunits = 'seconds since %s' % tref.strftime('%Y-%m-%d %H:%M:%S')
                # convert:
                values = netCDF4.num2date( values, tunits )
            
                # attributes:
                attrs = {}
                attrs['units'    ] = tunits
                attrs['long_name'] = 'time'
                
                # add:
                self.AppendPixels( varname, values, ('pixel',), attrs )

            # full track:
            elif special in ['track_longitude','track_corner_longitudes',\
                             'track_latitude' ,'track_corner_latitudes'] :
            
                # name of original field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # extract:
                odat = tropdata.GetVar( opath )

                # extact values:
                values = odat['data']
                
                # some differences between center and corner arrays:
                if 'corner' in special :
                    ## adhoc: swap from (corner,scan,pix) to (scan,pix,corner):
                    #values = values.transpose((1,2,0))
                    # set dimension names:
                    dimnames = ('track_image','track_pixel','corner')
                else :
                    # set dimension names:
                    dimnames = ('track_image','track_pixel')
                #endif
                
                # keep full swath, but select only range of scans with selected pixels:
                jj, = numpy.where( selected[itime,:,:].sum(axis=1) > 0 )
                j0 = jj.min()
                j1 = jj.max()
                
                # selection:
                vals = values[itime,j0:j1+1,:]
                
                # original units:
                ounits = odat['attrs']['units']
                # convert?
                if ounits != vunits :
                    self.logger.error( 'input variable "%s" units "%s" do no match with target variable "%s" units "%s"' % (oname,ounits,varname,vunits) )
                    raise Exception
                #endif
            
                # attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = odat['attrs']['long_name']
                
                # add:
                self.AppendTrack( varname, vals, dimnames, attrs )

                # store track info:
                self.trackinfo = {}
                self.trackinfo['iscan0'] = j0
                self.trackinfo['nscan' ] = vals.shape[0]
                self.trackinfo['npix'  ] = vals.shape[1]

            # full track:
            elif special == 'compress="track_image track_pixel"' :
            
                # extract attribute:
                clist = special.lstrip('compress="').rstrip('"')
            
                # create description:
                Y,X = clist.split()
                description = "original zero-based indices: `Y` = `P`/len(`X`), `X` = `P` mod len(`X`)"
                description = description.replace('X',X)
                description = description.replace('Y',Y)
                description = description.replace('P',varname)

                # attributes:
                attrs = {}
                attrs['units'      ] = vunits
                attrs['long_name'  ] = 'compressed pixel coordinate within track'
                attrs['compress'   ] = clist
                attrs['description'] = description
            
                # numbered index:
                values = numpy.array( (iit - self.trackinfo['iscan0']) * self.trackinfo['npix'] + iip, dtype='i4' )
                
                ## testing ...
                #yy = numpy.floor( values / self.track['npix'] )
                #xx = numpy.floor( values % self.track['npix'] )
                #for k in range(len(iit)) :
                #    print( 'xxx %6i  [%4i,%2i]  [%4i,%2i]' % (k,iit[k]-self.track['iscan0'],iip[k],yy[k],xx[k]) )
                ##endfor
                
                # add:
                self.AppendPixels( varname, values, ('pixel',), attrs )
            
                
            # pressure
            elif special == 'hybounds_to_pressure' :
            
                ## set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'half level pressures'
            
                # original fields:
                odat = {}
                for key in ['sp','hyam','hybm'] :
                    # full name:
                    opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                    # extract:
                    odat[key] = tropdata.GetVar( opath )
                #endfor
                
                # extract pixels:
                sp = odat['sp']['data'].squeeze()[iit,iip]
                
                # half level pressures:
                phalf = numpy.zeros((npix,tropdata.nlayer+1),float)
                # surface:
                iz = 0
                phalf[:,iz] = odat['hyam']['data'][iz,0] + odat['hybm']['data'][iz,0] * sp 
                # loop over layers:
                for iz in range(1,tropdata.nlayer+1) :
                    # top pressure:
                    phalf[:,iz] = odat['hyam']['data'][iz-1,1] + odat['hybm']['data'][iz-1,1] * sp
                #endfor
                
                # append:
                self.AppendPixels( varname, phalf, ('pixel','layer_interface'), attrs )

            # tropospheric kernel
            elif special == 'kernel_trop' :
            
                # original fields:
                odat = {}
                for key in ['avk','amf','amft'] :
                    # full name:
                    opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                    # extract:
                    odat[key] = tropdata.GetVar( opath )
                #endfor
                
                # set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'averaging kernel'
                
                # extract subset:
                for key in ['avk'] :
                    # dimension names
                    odimnames = odat[key]['dimnames']
                    #odimnames = odat['dimnames']
                    # location of layer index might be different;
                    # target should (levels,pixels)
                    if odimnames == ('time','scanline','ground_pixel','layer') :
                        odat[key]['data'] = (odat[key]['data'][itime,iit,iip,:]).T
                    else:
                        self.logger.error( 'unsupported dimnames "%s" for special "%s"' % (odimnames,special) )
                        raise Exception
                    #endif
                #endfor
                for key in ['amf','amft'] :
                    odat[key]['data'] = odat[key]['data'][itime,iit,iip]
                #endfor

                # convert, result should have shape (pixels,levels) :
                avk = numpy.zeros( (len(iit),tropdata.nlayer), float )
                for k in range(tropdata.nlayer) :
                    avk[:,k] = odat['avk']['data'][k,:] * odat['amf']['data'] / odat['amft']['data']
                #endfor
                
                # is troposphere layer defined ?
                key = 'troplayer'
                opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                # defined ?
                if len(opath) > 0 :
                    # split:
                    ofield,oname = opath.split('/')
                    # extract data:
                    odat[key] = tropdata.GetVar( opath )
                    # loop over pixels:
                    for ipix in range(npix) :
                        # layer index (1-based?):
                        trpl = int(odat[key]['data'][itime,iit[ipix],iip[ipix]])
                        # check ..
                        if (trpl < 1) or (trpl > tropdata.nlayer) :
                            self.logger.error( 'found tropo layer %i outside range 1,..,%i' % (trpl,tropdata.nlayer) )
                            raise Exception
                        #endif
                        # reset to zero above layer holding tropopause:
                        avk[ipix,trpl:] = 0.0
                    #endfor
                #endif
                
                # append:
                self.AppendPixels( varname, avk, ('pixel','layer'), attrs )

            # hybride coeff:
            elif special == 'hyblevel' :
            
                # index values:
                values = numpy.arange(nlayer,dtype='i4') + 1
                # attributes:
                attrs = {}
                attrs['units'        ] = vunits
                attrs['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
                attrs['formula'      ] = 'p(n,k,i) = ap(k) + b(k)*ps(n,i)'
                attrs['formula_terms'] = 'ap: hyam b: hybm ps: surface_pressure'
                
                # (re)create variable:
                self.AddVariable( varname, values, ('layer',), attrs )

            else :
            
                self.logger.error( 'unsupported special "%s"' % special )
                raise Exception
            
            #endif
        
        #endfor # target variables                

    #enddef AddSelection

#endclass TropomiExtract


#-------------------------------------------------
# OMI convert task
#-------------------------------------------------

class EmipTropomiConvert( emip.EmipTask ) :

    """
    EMIP task to extract pixels from original TROPOMI files
    and save them in a decent format.

    The arguments are passed to the base class to read the file
    with settings. 

    The first setting that is read is the time range over which 
    files should be converted::

      <rcbase>.timerange.start        :  2012-01-01 00:00
      <rcbase>.timerange.end          :  2012-12-31 23:59

    A loop over this time range will be performed with steps per day.
    For each day, the content of an archive directory is scanned;
    specify the directory using time templates::

      <rcbase>.files.dir  :  /data/TEMIS/airpollution/no2col/data/tropomi/%Y/%m

    The original data files are named::

        S5P_OFFL_L2__NO2____20180701T024100_20180701T042230_03699_01_010002_20180707T040951.nc
        |            |      |               |               |     |  |      |
        product      comp   starttime       endtime         orbit v1 v2     productiontime

    The first elements the filenames should be specified to allow
    filtering on the correct files::

      <rcbase>.files.product                  :  S5P_OFFL_L2

    Specifiy the component part of the filename with::
    
      <rcbase>.component        :  NO2

    The TROPOMI files are in NetCDF4 format.
    An example of the header is available from:
    
    * `doc/samples/S5P_OFFL_L2__NO2____20180701T005930_20180701T024100_03698_01_010002_20180707T022838.txt <../../samples/S5P_OFFL_L2__NO2____20180701T005930_20180701T024100_03698_01_010002_20180707T022838.txt>`_

    The :py:meth:`.TROPOMI_NC4_File` class is used to read the file;
    see the class documentation for the most relevant content.
    
    xxxxxxxxxxxxxxxxx
    
    Provide the name of the variable that contains the actual product:

      <rcbase>.data.column      :  TroposphericVerticalColumn

    The :py:meth:`.OMI_HDF5_File.ConvertTime` method of the object
    with the the OMI data is called to create decent time object.

    The :py:meth:`.OMI_HDF5_File.SelectPixels` method is called
    to create a pixel selection mask.

    The converted data is stored in an :py:class:`.OmiExtract` object
    which will take care of writing to a netCDF file.
    The :py:meth:`.OmiExtract.AddSelection` method is used to process
    the selected pixels, e.g. selecting the variables, apply conversions, etc.
    
    The converted data is written to a file specified by directory
    and filename templates.

    The output directory could include templates for time values::

      <rcbase>.output.dir    :  /data/OMI-selection/%Y/%m

    Also the filenames could include time values,
    and in addition a template '``%{orbit}``' in which the orbit id
    from the OMI file is inserted::

      <rcbase>.output.filename    :  OMI-Aura_NO2_%Y%m%d_%{orbit}.nc

    Specify global attributes and their value with::

      <rcbase>.output.attrs  :  format Conventions author institution email

      <rcbase>.output.attr.format         :  1.0
      <rcbase>.output.attr.Conventions    :  CF-1.6
      <rcbase>.output.attr.author         :  Arjo Segers
      <rcbase>.output.attr.institution    :  MetNorway, Oslo, Norway
      <rcbase>.output.attr.email          :  Arjo.Segers@met.no

    """
    
    def __init__( self, rcfile, rcbase='', env={}, indent='' ) :
    
        """
        Extract TROPOMI observation and write to NetCDF file.
        
        """
                        
        # modules:
        import os
        import datetime
        import fnmatch
        import tarfile
        import numpy

        # tools:
        import TROPOMI
        
        ## testing ...
        #import importlib
        #importlib.reload( TROPOMI )

        # init base object:
        emip.EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # time range:
        t1 = self.GetSetting( 'timerange.start', totype='datetime' )
        t2 = self.GetSetting( 'timerange.end'  , totype='datetime' )
        
        # info ...
        tfmt = '%Y-%m-%d %H:%M'
        self.logger.info( indent+'timerange: [%s,%s]' % (t1.strftime(tfmt),t2.strftime(tfmt)) )

        # input directory:
        inputdir = self.GetSetting( 'files.dir' )

        # filename filters:
        fnfilters = self.GetSetting( 'files.filters' ).split()

        # target component:
        component = self.GetSetting( 'component' )

        # renew existing files?
        renew = self.GetSetting( 'renew', totype='bool' )

        # info ...
        self.logger.info( indent+'collect input file names ...' )
        # init list with input filenames:
        filenames = []
        # content of input directories is listed,
        # name of directory changes in time:
        indir_prev = None
        # loop over time range:
        t = t1
        while t <= t2 :
        
            # info ...
            self.logger.info( indent+'  %s ...' % t.strftime('%Y-%m-%d') )
            
            # resolve time templates:
            indir = t.strftime(inputdir)
            # new?
            if indir != indir_prev :
                # present ?
                if os.path.isdir( indir ) :
                    # list:
                    fnames = os.listdir( indir )
                    # sort (in place)
                    fnames.sort()
                    # loop over files:
                    for fname in fnames :
                        # full path:
                        filename = os.path.join(indir,fname)
                        # tar file?
                        if fname.endswith('.tar') :
                            # list file content:
                            with tarfile.open(filename,'r') as tar :
                                xnames = tar.getnames()
                            #endwith
                            # loop:
                            for xname in xnames :
                                # add tarfile and content (to be extracted):
                                filenames.append( filename+'%'+xname )
                            #endfor
                        else :
                            # add full path:
                            filenames.append( filename )
                        #endif # tarfile
                    #endfor
                    # reset:
                    indir_prev = indir
                #endif # present
            #endif # new dir
            # next:
            t = t + datetime.timedelta(1)
        #endwhile

        # info ...
        self.logger.info( indent+'input dir(s)    : %s' % inputdir )
        self.logger.info( indent+'number of files : %i' % len(filenames) )

        # target file:
        xomi_dir               = self.GetSetting( 'output.dir' )
        xomi_filename_template = self.GetSetting( 'output.filename' )
        
        # not yet ...
        xomi = None

        # attributes
        attrs = {}
        attrnames = self.GetSetting( 'output.attrs' ).split()
        for key in attrnames :
            value = self.GetSetting( 'output.attr.%s' % (key) )
            #if key.startswith('Window_') :
            #    value = float(value)
            #elif key.startswith('Grid_dimension_') or \
            #     key.startswith('Number_of') :
            #    value = int(value)
            ##endif
            attrs[key] = value
        #endfor

        ## max layers:
        #nlayer = self.GetSetting( 'leip.prod.xomi.max_nlayer_trop', 'int' )

        # info ..
        self.logger.info( indent+'loop over files ...' )

        # wait ?
        wait_for = self.GetSetting( 'wait.for' )
        # defined ?
        wait = len(wait_for) > 0

        # skip some ?
        skiplist = self.GetSetting( 'skip' ).split()
        
        # loop:
        for filename in filenames :

            # basename:
            fname = os.path.basename( filename )
            
            # split in "<tarfile>%<fname>" if necessary:
            if '%' in fname :
                tarf,fname = filename.split('%',1)
            else :
                tarf = None
            #endif

            # filter:
            found = False
            for fnfilter in fnfilters :
                found = found or fnmatch.fnmatch( fname, fnfilter )
                if found : break
            #endfor
            if not found : continue

            # wait ?
            if wait :
                # now match ?
                if fname == wait_for :
                    wait = False
                else :
                    continue
                #endif
            #endif

            # info ...
            self.logger.info( indent+'  file %s ...' % fname )

            # adhoc skip:
            if fname in skiplist :
                # info ...
                self.logger.warning( indent+'    file is in skip list; continue ...' )
                # next:
                continue
            #endif
            
            # filenames:
            #    S5P_OFFL_L2__NO2____20180701T005930_20180701T024100_03698_01_010002_20180707T022838.nc
            bname = os.path.basename(filename).replace('.nc','')
            # split:
            product_comp,timesetc = bname.split('____')
            starttime,endtime,orbit,key1,key2,prodtime = timesetc.split('_')
            
            # extract date:
            t0 = datetime.datetime.strptime(starttime,'%Y%m%dT%H%M%S')            
            # filter ...
            if (t0 < t1) or (t0 > t2) :
                # info ...
                self.logger.warning( '    file outside time range; continue ...' )
                # skip:
                continue
            #endif

            # convert orbit number:
            orbnr = int(orbit)

            # extract from tarfile first?
            if tarf is not None :
                # info ...
                self.logger.info( indent+'    extract from tarfile ...' )
                # extract:
                with tarfile.open(tarf,'r') as tar :
                    tar.extract( fname )
                #endwith
                # info ...
                self.logger.info( indent+'    read ...' )
                # read file:
                tropfile = TROPOMI.TROPOMI_NC4_File( fname )
                # remove temporary extract:
                os.remove( fname )
            else :
                # info ...
                self.logger.info( indent+'    read ...' )
                # read file:
                tropfile = TROPOMI.TROPOMI_NC4_File( filename )
            #endif

            # apply selections, return bool mask and list of history lines:
            selected,history = tropfile.SelectPixels( self.rcf, self.rcbase )
            # count:
            nselected = selected.sum()

            # info ...
            self.logger.info( indent+'    selected %i of %i' % (nselected,selected.size) )
            # any ?
            if nselected > 0 :

                # init:
                xomi = TropomiExtract()
                # add:
                xomi.AddSelection( tropfile, selected, self.rcf, self.rcbase, orbnr )

                # update history:
                history.append( 'added %i pixels from %s' % (nselected,os.path.basename(fname)) )
                
                # time average:
                taver = xomi.GetTimeAverage()
                # round to nearby hour:
                hour = numpy.round( taver.hour + taver.minute/60.0 )
                tstamp = taver.replace(hour=0,minute=0,second=0) + datetime.timedelta(0,hour*3600)

                # target dir:
                outdir = taver.strftime(xomi_dir)
                # target file:
                xomi_filename = xomi_filename_template
                xomi_filename = xomi_filename.replace( '%{orbit}', orbit )
                xomi_filename = tstamp.strftime(xomi_filename)
                # combine:
                xomi_filename = os.path.join( outdir, xomi_filename )

                # info ...
                self.logger.info( indent+'    create extract file %s ...' % os.path.basename(xomi_filename) )

                # keep existing?
                if os.path.isfile(xomi_filename) and (not renew) :
                    # info ..
                    self.logger.info( indent+'    output file already present; continue ...' )
                    # next:
                    continue
                #endif

                # write:
                xomi.Write( filename=xomi_filename, attrs=attrs, history=history )

            #endif  # any selected
            
            ## testing ...
            #break

        #endfor  # OMI files

        ## show:
        #print ''
        #mdf.show( xomi_filename )

    #enddef __init__

#endclass EmipTropomiConvert


#-------------------------------------------------
# end
#-------------------------------------------------

