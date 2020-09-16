
#-------------------------------------------------
# help
#-------------------------------------------------

"""
*******************
``emip_omi`` module
*******************

The :py:mod:`emip_omi` module provides classes
for processing OMP data.


Class hierchy
=============

The classes and are defined according to the following hierchy:


* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

    * :py:class:`.EmipTask`

      * :py:class:`.EmipOmiConvert`

  * :py:class:`.SatXFile`

    * :py:class:`.OmiExtract`


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


class OmiExtract( satx.SatXFile ) :

    """
    Storage for converted OMI file.
    """
    
    def AddSelection( self, omidata, selected, rcf, rckey, orbnr, indent='' ) :
    
        """
        Add selected OMI to satellite extract file.
        
        Arguments:
        
        * ``omidata``    :  :py:class:`.OMI_HDF5_File` object
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
        import numpy
        
        # info ..
        self.logger.info( indent+'add OMI pixels to extract file ...' )
        
        # copy dimension:
        #self.nlayer = omidata.nlayer
        nlayer = omidata.nlayer
        
        # attributes to be copied:
        self.track_attrnames = rcf.get( rckey+'.output.track_attrs' ).split()
        # existing names:
        trattrs = omidata.swaths[omidata.component]['attrs']
        # copy:
        attrs = {}
        for key in self.track_attrnames :
            # original name:
            okey = rcf.get( rckey+'.output_track_attr.%s' % key )
            # check ..
            if okey not in trattrs.keys() :
                self.logger.error( 'key "%s" not found in track attributes:' % okey )
                for trkey in trattrs.keys() :
                    self.logger.error( '  %s = %s' % (trkey,trattrs[trkey]) )
                #endfor
                raise Exception
            #endif
            # get:
            avalue = omidata.swaths[omidata.component]['attrs'][okey]
            if type(avalue) == numpy.ndarray : avalue = avalue[0]
            if type(avalue) == numpy.bytes_ : avalue = avalue.decode('utf-8')
            # copy:
            attrs[key] = avalue
        #endfor
        # store:
        self.track_attrs.append( attrs )
        
        # selection indices:
        iit,iip = numpy.where( selected )
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
            #if special in ['None','time'] :
            if special == 'None' :
            
                # name of original field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # split:
                ofield,oname = opath.split('/')
                # extract:
                odat = omidata.swaths[omidata.component][ofield][oname]

                # extact values:
                values = odat['data']
                # original units etc:
                ounits    = odat['units']
                olongname = odat['long_name']
                odimnames = odat['dimnames']
                
                # copy units?
                if vunits == '__native__' : vunits = ounits
                
                # adhoc fix: pressure units say 'hPa', but seems 'Pa' ..
                if (oname == 'TM4SurfacePressure') and (ounits == 'hPa') and (values.max() > 1100.0) :
                    # reset units:
                    ounits = 'Pa'
                #endif
                
                # conversion ?
                if ounits != vunits :
                    # required conversion:
                    conversion = '%s -> %s' % (ounits,vunits)
                    # switch:
                    if conversion in ['deg -> degrees_east','deg -> degrees_north',
                                      'Degrees -> degrees_east', 'Degrees -> degrees_north',
                                      'NoUnits -> 1'] :
                        # no conversion needed
                        pass
                    elif conversion in ["molecules cm^-2 -> 1e15 cm**-2",
                                        "molecules/cm^2 -> 1e15 cm**-2" ] :
                        # apply factor:
                        values = values * 1.0e-15
                    elif conversion in ["hPa -> Pa"] :
                        # apply factor:
                        values = values * 1.0e2
                    elif conversion in ["NoUnits -> %"] :
                        # assume original is in [0,1], convert to % :
                        values = values * 100.0
                    else :
                        self.logger.error( 'unsupported conversion "%s" requested ...' % conversion )
                        raise Exception
                    #endif
                #endif
                
                # extract 1D array of pixels from 2D track:
                if odimnames in [('time','pixel')] :

                    values = values[iit,iip]
                    dimnames = ('pixel',)

                elif odimnames in [('corner','time','pixel')] :

                    values = values[:,iit,iip].T
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

            # full track:
            elif special in ['track_longitude','track_corner_longitudes',\
                             'track_latitude' ,'track_corner_latitudes'] :
            
                # name of original field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # split:
                ofield,oname = opath.split('/')
                # extract:
                odat = omidata.swaths[omidata.component][ofield][oname]

                # extact values:
                values = odat['data']
                
                # some differences between center and corner arrays:
                if 'corner' in special :
                    # adhoc: swap from (corner,scan,pix) to (scan,pix,corner):
                    values = values.transpose((1,2,0))
                    # set dimension names:
                    dimnames = ('track_image','track_pixel','corner')
                else :
                    # set dimension names:
                    dimnames = ('track_image','track_pixel')
                #endif
                
                # keep full swath, but select only range of scans with selected pixels:
                jj, = numpy.where( selected.sum(axis=1) > 0 )
                j0 = jj.min()
                j1 = jj.max()
                
                # selection:
                vals = values[j0:j1+1,:]
                
                # original units:
                ounits = odat['units']
                # convert?
                if ounits != vunits :
                    # accept ...
                    if ounits in ['deg','Degrees'] :
                        pass
                    else :
                        self.logger.error( 'input variable "%s" units "%s" do no match with target variable "%s" units "%s"' % (oname,ounits,varname,vunits) )
                        raise Exception
                    #endif
                #endif
            
                # attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = odat['long_name']
                
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
            elif special == 'hym_to_pressure' :
            
                # set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'half level pressures'
            
                # original fields:
                odat = {}
                for key in ['sp','hyam','hybm'] :
                    # full name:
                    opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                    # split:
                    ofield,oname = opath.split('/')
                    # extract:
                    odat[key] = omidata.swaths[omidata.component][ofield][oname]
                #endfor
                
                # convert pressure units if necessary:
                for key in ['sp','hyam'] :
                    # different from target ?
                    if odat[key]['units'] != vunits :
                        # conversion:
                        conversion = '%s -> %s' % (odat[key]['units'],vunits)
                        # which ?
                        if conversion == 'hPa -> Pa' :
                            # apply factor:
                            odat[key]['data'] = odat[key]['data'] * 100.0
                        else :
                            self.logger.error( 'unsupported conversion "%s"' % conversion )
                            raise Exception
                        #endif
                    #endif
                #endfor  # pressure variables
                
                # extract pixels:
                sp = odat['sp']['data'][iit,iip]
                
                # mid level pressures:
                pmid = numpy.zeros((npix,omidata.nlayer),float)
                for iz in range(omidata.nlayer) :
                    pmid[:,iz] = odat['hyam']['data'][iz] + odat['hybm']['data'][iz] * sp
                #endfor
                
                # extrapolate to interface:
                phalf = numpy.zeros((npix,omidata.nlayer+1),float)
                for iz in range(omidata.nlayer-1,-1,-1) :
                    phalf[:,iz] = 2.0 * pmid[:,iz] - phalf[:,iz+1]
                #endfor
                
                # extend variable with latest pixels:
                self.AppendPixels( varname, phalf, ('pixel','layer_interface'), attrs )

            # tropospheric kernel
            elif special == 'kernel_trop' :
            
                # original fields:
                odat = {}
                for key in ['avk','amf','amft'] :
                    # full name:
                    opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                    # split:
                    ofield,oname = opath.split('/')
                    # extract:
                    odat[key] = omidata.swaths[omidata.component][ofield][oname]
                #endfor
                
                # set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'averaging kernel'
                
                # extract subset:
                for key in ['avk'] :
                    # dimension names
                    odimnames = odat[key]['dimnames']
                    # location of layer index might be different;
                    # target should (levels,pixels)
                    if odimnames == ('time','pixel','layer') :
                        odat[key]['data'] = (odat[key]['data'][iit,iip,:]).T
                    elif odimnames == ('layer','time','pixel') :
                        odat[key]['data'] = odat[key]['data'][:,iit,iip]
                    else :
                        self.logger.error( 'unsupported dimnames "%s" for special "%s"' % (odimnames,special) )
                        raise Exception
                    #endif
                #endfor
                for key in ['amf','amft'] :
                    odat[key]['data'] = odat[key]['data'][iit,iip]
                #endfor
                
                # convert, result should have shape (pixels,levels) :
                avk = numpy.zeros( (len(iit),omidata.nlayer), float )
                for k in range(omidata.nlayer) :
                    avk[:,k] = odat['avk']['data'][k,:] * odat['amf']['data'] / odat['amft']['data']
                #endfor
                
                # is troposphere layer defined ?
                key = 'troplayer'
                opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
                # defined ?
                if len(opath) > 0 :
                    # split:
                    ofield,oname = opath.split('/')
                    # extract original swath:
                    troplayer = omidata.swaths[omidata.component][ofield][oname]['data']
                    # check ..
                    if troplayer.max() > omidata.nlayer :
                        self.logger.error( 'found tropo layer %i above specified maximum %i' % (troplayer.max(),omidata.nlayer) )
                        raise Exception
                    #endif
                    # trancate to zero above the tropopause layer:
                    for ipix in range(npix) :
                        trpl = int(troplayer[iit[ipix],iip[ipix]])
                        avk[ipix,trpl:] = 0.0
                    #endfor
                #endif
                
                # add:
                self.AppendPixels( varname, avk, ('pixel','layer'), attrs )

            # hybride coeff:
            elif special == 'hym_to_hyi' :
            
                # copy:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'half level hybride coefficients'

                # name of original field:
                opath = rcf.get( '%s.output.var.%s.from' % (rckey,varname) )
                # split:
                ofield,oname = opath.split('/')
                # extract:
                odat = omidata.swaths[omidata.component][ofield][oname]
                # different from target ?
                if odat['units'] != vunits :
                    # conversion:
                    conversion = '%s -> %s' % (odat['units'],vunits)
                    # which ?
                    if conversion == 'NoUnits -> 1' :
                        # unit free
                        pass
                    elif conversion == 'hPa -> Pa' :
                        # apply factor:
                        odat['data'] = odat['data'] * 100.0
                    else :
                        self.logger.error( 'unsupported conversion "%s"' % conversion )
                        raise Exception
                    #endif
                #endif
                # convert full arrays:
                nf = len(odat['data'])
                ndat = numpy.zeros((nf+1),float)
                for k in range(nf-1,-1,-1) :
                    ndat[k] = 2*odat['data'][k] - ndat[k+1]
                #endfor
                # round:
                if varname == 'hyai' : ndat = ndat.clip(0.0)
                if varname == 'hybi' : ndat = ndat.clip(0.0,1.0)
                
                # fill:
                self.AddVariable( varname, ndat[0:omidata.nlayer+1], ('layer_interface',), attrs )

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

            # orbit number
            elif special == 'orbit_number' :
            
                # attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'orbit number'
                
                # fill array with orbit numbers:
                nt = len(iit)
                vv = numpy.ones((nt),dtype='i4') * orbnr
                
                # add:
                self.AppendPixels( varname, vv, ('pixel',), attrs )

            ## index for layers
            #elif special == 'index' :
            #  
            #    # new?
            #    if varname not in self.var.keys() :
            #       #init:
            #       self.var[varname] = {}
            #       #copy:
            #       self.var[varname]['units'    ] = vunits
            #       self.var[varname]['long_name'] = 'index'
            #       self.var[varname]['special'  ] = special
            #       #needs to be filled only once:
            #       self.var[varname]['data'     ] = numpy.arange(omidata.nlayer) + 1
            #    #endif

            # orbit number
            elif special == 'scan_number' :
            
                # set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'scan number'
                
                # add time indices (scan lines), convert to 1-based:
                self.AppendPixels( varname, numpy.array(iit+1,dtype='i4'), ('pixel',), attrs )

            # pixel within scan line
            elif special == 'pixel_number' :
            
                # set attributes:
                attrs = {}
                attrs['units'    ] = vunits
                attrs['long_name'] = 'pixel number'
                
                # add pixel indices, convert to 1-based:
                self.AppendPixels( varname, numpy.array(iip+1,dtype='i4'), ('pixel',), attrs )

            ## gridded average
            #elif special == 'gridaver' :
            #
            #    # new ?
            #    if varname not in self.var.keys() :
            #        # init:
            #        self.var[varname] = {}
            #        # domain in [west,east,south,north]
            #        line = rcf.get( '%s.output.var.%s.domain'   % (rckey,varname) )
            #        self.var[varname]['domain'  ] = list(map(float,line.split()))
            #        # gridsize in [nlon,nlat]
            #        line = rcf.get( '%s.output.var.%s.gridsize' % (rckey,varname) )
            #        self.var[varname]['gridsize'] = list(map(int,line.split()))
            #        # store grid definition i fnot done yet:
            #        if self.grid == None :
            #            # extract:
            #            west,east,south,north = self.var[varname]['domain']
            #            nlon,nlat = self.var[varname]['gridsize']
            #            # fill:
            #            self.grid = { 'nlon' : nlon, 'nlat' : nlat, 
            #                          'lons' : west  + numpy.arange(nlon)*(east -west )/(nlon-1.0),
            #                          'lats' : south + numpy.arange(nlat)*(north-south)/(nlat-1.0) }
            #        #endif
            #        # copy:
            #        self.var[varname]['units'    ] = vunits
            #        self.var[varname]['long_name'] = 'pixel number'
            #        self.var[varname]['special'  ] = special
            #        # init zero sum and counter:
            #        nlon,nlat = self.var[varname]['gridsize']
            #        self.var[varname]['data'  ] = numpy.zeros((nlat,nlon),float)
            #        self.var[varname]['weight'] = numpy.zeros((nlat,nlon),float)
            #    #endif
            #
            #    # extract data field:
            #    odat = {}
            #    for key in ['from','clons','clats','vza'] :
            #        # full name:
            #        opath = rcf.get( '%s.output.var.%s.%s' % (rckey,varname,key) )
            #        # split:
            #        ofield,oname = opath.split('/')
            #        # extract:
            #        odat[key] = omidata.swaths[omidata.component][ofield][oname]
            #    #endfor
            #
            #    # different from target ?
            #    if odat['from']['units'] != vunits :
            #        # conversion:
            #        conversion = '%s -> %s' % (odat['from']['units'],vunits)
            #        # which ?
            #        if conversion in ["molecules cm^-2 -> 1e15 cm**-2",
            #                           "molecules/cm^2 -> 1e15 cm**-2" ] :
            #            # apply factor:
            #            odat['from']['data'] = odat['from']['data'] * 1.0e-15
            #        else :
            #            self.logger.error( 'unsupported conversion "%s"' % conversion )
            #            raise Exception
            #        #endif
            #    #endif
            #    
            #    # loop over pixels:
            #    for ipix in range(npix) :
            #    
            #        # original indices:
            #        it = iit[ipix]
            #        ip = iip[ipix]
            #        
            #        # observation weight:
            #        vza = odat['vza']['data'][it,ip]
            #        obsweight = numpy.cos( numpy.radians(vza) )**2 
            #        # check ...
            #        if obsweight < 1.e-6 :
            #            self.logger.error( 'observation weight zero or negative:' )
            #            self.logger.error( '  pixel   :  %i %i' % (imag,ipix) )
            #            self.logger.error( '  vza     :  %f  deg' % vza )
            #            self.logger.error( '  weight  :  %f' % obsweight )
            #            raise Exception
            #        #endif
            #
            #        # indices and weights:
            #        inds = self.GridListRegion( odat['clons']['data'][:,it,ip], 
            #                                    odat['clats']['data'][:,it,ip],  
            #                                    self.var[varname]['gridsize'],
            #                                    self.var[varname]['domain'] )
            #        # loop over cells indices:
            #        for i,j in inds :
            #             # add contributions:
            #             self.var[varname]['data'  ][j,i] += obsweight * odat['from']['data'][it,ip]
            #             self.var[varname]['weight'][j,i] += obsweight
            #        #endfor
            #
            #        ## testing ...
            #        #break
            #
            #    #endfor  # pixels
                
            else :
            
                self.logger.error( 'unsupported special "%s"' % special )
                raise Exception
            
            #endif
        
        #endfor # target variables                

    #enddef AddSelection

#endclass OmiExtract


#-------------------------------------------------
# OMI convert task
#-------------------------------------------------

class EmipOmiConvert( emip.EmipTask ) :

    """
    EMIP task to extract pixels from original OMI files
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

      <rcbase>.files.dir  :  /data/TEMIS/airpollution/no2col/data/omi/data_v2/%Y

    The original data files are named::

        OMI-Aura _ L2-OMDOMINO _ 2003m0601t0105 - o00014 _ v001 - 2013m0320t110311.he5
        |          |             |                |        |      |
        project    product       starttime        orbit  version  production

    The first elements the filenames should be specified to allow
    filtering on the correct files::

      <rcbase>.files.project                  :  OMI-Aura
      <rcbase>.files.product                  :  L2-OMDOMINO

    If a file following this specification is found it is read into an
    :py:class:`.OMI_HDF5_File` object.
    The OMI files are in HDFEOS format and are assumed to have
    the following structure of groups and variables::

      /HDFEOS
             /SWATHS
                    /<swathgroup>
                                /Data Fields
                                                   /<component>RetrievalVerticalColumn
                                                      .title
                                                      .units
                                                   /...
                                /Geolocation Fields
                                                   /Latitude
                                                   /Longitude
                                                   /...

    Information on what can be expected in these files should be provided
    in the settings::

      <rcbase>.component        :  NO2
      <rcbase>.data.swathgroup  :  DominoNO2

    Also provide the name of the variable that contains the actual product:

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
        Extract OMI observation and write to NetCDF file.
        
        """
                        
        # modules:
        import os
        import datetime

        # tools:
        import OMI
        
        ## testing ...
        #import importlib
        #importlib.reload( OMI )

        # init base object:
        emip.EmipTask.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # time range:
        t1 = self.GetSetting( 'timerange.start', totype='datetime' )
        t2 = self.GetSetting( 'timerange.end'  , totype='datetime' )

        # input directory:
        inputdir = self.GetSetting( 'files.dir' )

        # filename prefix:
        project = self.GetSetting( 'files.project' )
        product = self.GetSetting( 'files.product' )

        # full prefix:
        prefix = '%s_%s' % (project,product)

        # target component:
        component = self.GetSetting( 'component' )

        # name of swathgroup:
        swathgroup = self.GetSetting( 'data.swathgroup' )
        # name of column product:
        columnname = self.GetSetting( 'data.column' )
        
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
        
            ## info ...
            #self.logger.info( indent+'%s ...' % t.strftime('%Y-%m-%d') )
            
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
                    # add full paths:
                    for fname in fnames :
                        filenames.append( os.path.join(indir,fname) )
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

            # filter:
            if not fname.startswith(prefix) : continue
            if not fname.endswith('.he5') : continue

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

            # remove extension:
            bname,ext = os.path.splitext(fname)
            # split:
            proj,prod,time_and_orbit,production = bname.split('_')
            # split:
            time,orbit = time_and_orbit.split('-')
            # extract date:
            t0 = datetime.datetime.strptime(time,'%Ym%m%dt%H%M')
            
            # filter ...
            if (t0 < t1) or (t0 > t2) :
                # info ...
                self.logger.warning( '    file outside time range; continue ...' )
                # skip:
                continue
            #endif

            # extract orbit number:
            orbnr = int(orbit[1:])

            # target dir:
            outdir = t0.strftime(xomi_dir)
            # target file:
            xomi_filename = xomi_filename_template
            xomi_filename = xomi_filename.replace( '%{orbit}', orbit )
            xomi_filename = t0.strftime(xomi_filename)
            # combine:
            xomi_filename = os.path.join( outdir, xomi_filename )
            
            # keep existing?
            if os.path.isfile(xomi_filename) and (not renew) :
                # info ..
                self.logger.info( indent+'    output file already present; continue ...' )
                # next:
                continue
            #endif

            # adhoc skip:
            if fname in skiplist :
                # info ...
                self.logger.warning( '    file is in skip list; continue ...' )
                # next:
                continue
            #endif

            # read:
            try :
                omidata = OMI.OMI_HDF5_File( filename, swathgroup, component, verbose=False )
            except OMI.OmiLevelException :
                self.logger.warning( '    problematic levels in "%s", skip ...' % fname )
                continue
            except IOError :
                self.logger.warning( '    file "%s" seems corruped, skip ...' % fname )
                continue
            except :
                raise
            #endtry

            # convert time field:
            omidata.ConvertTime()

            ## show map ...
            #if False :
            #    lons = omidata.swaths[omidata.component]['Geolocation Fields']['Longitude']['data']
            #    lats = omidata.swaths[omidata.component]['Geolocation Fields']['Latitude' ]['data']
            #    #vv = omidata.swaths[omidata.component]['Data Fields']['HCHORetrievalVerticalColumn']['data']
            #    plt.figure()
            #    plt.plot( lons, lats, 'r.' )
            #    plt.plot( [-15,35,35,-15,-15], [35,35,70,70,35], 'k-' )
            ##endif

            # apply selections, return bool mask and list of history lines:
            selected,history = omidata.SelectPixels( self.rcf, self.rcbase )
            # count:
            nselected = selected.sum()

            # info ...
            self.logger.info( indent+'    selected %i of %i' % (nselected,selected.size) )
            # any ?
            if nselected > 0 :

                # info ...
                self.logger.info( indent+'    create extract file %s ...' % os.path.basename(xomi_filename) )

                # update history:
                history.append( 'added %i pixels from %s' % (nselected,os.path.basename(fname)) )

                # init:
                xomi = OmiExtract()
                # add:
                xomi.AddSelection( omidata, selected, self.rcf, self.rcbase, orbnr )
                # write:
                xomi.Write( filename=xomi_filename, attrs=attrs, history=history )

            #endif  # any selected
            
        #endfor  # OMI files

        ## show:
        #print ''
        #mdf.show( xomi_filename )

    #enddef __init__

#endclass EmipOmiConvert



#-------------------------------------------------
# end
#-------------------------------------------------

