
"""
***************
``OMI`` module
***************

Access OMI data in native HDF-EOS-5 format.


Class hierchy
=============

The classes are defined according to the following hierchy:

* :py:class:`.OMI_HDF5_File`

The following exception classes are defined:

* :py:class:`.OmiLevelException`


Classes
=======

"""

########################################################################
###
### OMI data class
###
########################################################################

class OmiLevelException( Exception ) :
    """
    Exception raised if inconsequent numbers of levels are found.
    """
#endclass OmiLevelException

# *

class OMI_HDF5_File( object ) :

    """
    Base class to read OMI HDF5 data file.

    Assumed structure::
    
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
        

    Arguments:
    
    * ``filename``    : name of input file
    * ``swatchgroup`` : name of group, e.g. "NO2'
    * ``component``   : component name, e.g. 'NO2'

    The file content is stored in attributes:
    
    * ``swaths[component]``   : dictionairy with per-component the data and geolocation fields

    """
    
    def __init__( self, filename, swathgroup, component, verbose=False ) :
    
        """
        Read OMI data file.        
        """
        
        # modules:
        import logging
        import os
        import h5py
        import numpy
        
        # setup logger:
        self.logger = logging.getLogger( 'OMI' )
        
        # store:
        self.filename = filename
        self.component = component

        # check ...
        if not os.path.isfile(filename) :
            self.logger.error( 'file not found : %s' % filename )
            raise Exception
        #endif
        
        # open:
        try :
            hid = h5py.File( filename, 'r' )
        except IOError :
            self.logger.error( 'IOError from opening file : %s' % filename )
            raise
        except :
            self.logger.error( 'could not open file : %s' % filename )
            raise
        #endtry
        
        # init result:
        self.swaths = {}
        self.swaths[component] = {}
        
        # form group name:
        compgrpname = '/HDFEOS/SWATHS/%s' % swathgroup
        
        # get group:
        compgrp = hid.get( compgrpname )
        if compgrp == None :
            self.logger.error( 'group "%s" not found in %s' % (compgrpname,filename) )
            raise Exception
        #endif
        
        # copy attributes:
        self.swaths[component]['attrs'] = {}
        for key in compgrp.attrs.keys() :
            self.swaths[component]['attrs'][key] = compgrp.attrs[key]
        #endfor
        
        # loop over groups:
        for field in compgrp.keys() :
            # info ...
            if verbose : print( '  %s' % field )
            # init result:
            flds = {}
            # get field group:
            fieldgrp = compgrp.get( field )
            # loop over values:
            for varname in fieldgrp.keys() :
                # testing ...
                #if varname != 'Pressure' : continue
                # info ...
                if verbose : print( '    %s' % varname )
                # init result:
                fld = {}
                # handle to variable:
                varid = fieldgrp.get(varname)
                # info ...
                if verbose :
                    print( '      data type : ', varid.dtype )
                    print( '      shape     : ', varid.shape )
                    print( '      attributes:' )
                    for key in varid.attrs.keys() :
                        print( '        %s = %s' % (key,varid.attrs[key]) )
                    #endfor  # attributes
                #endfor
                # data; geo-stationair files have (sometimes?) only a single time value:
                if len(varid.shape) == 0 :
                    values = numpy.array( varid.value )
                else :
                    values = varid[:]
                #endif
                # mask for the unfilled values:
                _FillValue = varid.attrs['_FillValue']
                mask = values == _FillValue
                # adhoc fix: also mask the 'nan' values:
                if numpy.any(numpy.isnan(values)) :
                    # info ...
                    self.logger.warning( '    found NaN values in "%s" in "%s"' % (varname,filename) )
                    # mask also these values:
                    mask = mask | numpy.isnan(values)
                #endif
                # created masked array:
                values = numpy.ma.array( data=values, mask=mask )
                # unpack:
                add_offset = varid.attrs['Offset']
                scale_factor = varid.attrs['ScaleFactor']
                values = add_offset + scale_factor * values
                # store:
                fld['data'      ] = values
                # extract attribute:
                avalue = varid.attrs['Title']
                if type(avalue) == numpy.ndarray : avalue = avalue[0]
                fld['long_name' ] = avalue
                # extract attribute:
                avalue = varid.attrs['Units']
                if type(avalue) == numpy.ndarray : avalue = avalue[0]
                fld['units'] = avalue.decode('latin1')
                # extract attribute:
                avalue = varid.attrs['_FillValue']
                if type(avalue) == numpy.ndarray : avalue = avalue[0]
                fld['_FillValue'] = avalue
                # fill units from long_name if possible:
                if fld['units'] in ['NoUnits'] :
                    # replace if possible:
                    if varname == 'CloudRadianceFraction' :
                        # problem: longname says '%', but sometimes
                        # files are found with with 0-1 (ISOTROP!)
                        if numpy.nanmax(fld['data']) <= 1.0 :
                            fld['units'] = '1'
                        else : 
                            fld['units'] = '%'
                        #endif
                    elif ('AirMassFactor'   in varname) or \
                         ('AveragingKernel' in varname) or \
                         ('CloudFraction'   in varname) or \
                         ('Id'              in varname) or \
                         ('Flag'            in varname) or \
                         ('Albedo'          in varname) or \
                         ('SurfaceAlbedo'   in varname) or \
                         ('PressurelevelB'  in varname) or \
                         ('TropoPauseLevel' in varname) or \
                         ('volumetric mixing ratio' in fld['long_name']) :
                        fld['units'] = '1'
                    else :
                        self.logger.error( 'could not replace units "%s" for variable "%s" based on on long_name "%s"' \
                                % (fld['units'],varname,fld['long_name']) )
                        self.logger.error( 'data range : %f - %f' % (numpy.nanmin(fld['data']),numpy.nanmax(fld['data'])) )
                        raise Exception
                    #endif
                #endif
                # store:
                flds[varname] = fld
            #endfor  # variables
            # store:
            self.swaths[component][field] = flds
        #endfor # field groups
        
        # done:
        hid.close()
        
        # extract time dimension:
        shp = self.swaths[self.component]['Geolocation Fields']['Time']['data'].shape
        # single value for geo, time line for leo:
        if len(shp) == 0 :
            # get shape of longitude field:
            lon_shp = self.swaths[self.component]['Geolocation Fields']['Longitude']['data'].shape
            # assume leading is time:
            self.ntime = lon_shp[0]
            # extend time array:
            time_value = self.swaths[self.component]['Geolocation Fields']['Time']['data']
            self.swaths[self.component]['Geolocation Fields']['Time']['data'] = numpy.ones((self.ntime),time_value.dtype) + time_value
        elif len(shp) == 1 :
            self.ntime = shp[0]
        else :
            self.logger.error( 'expected single value time (GEO) or 1D (LE), found shape : %s' % str(shp) )
            raise Exception
        #endif

        # extract scan dimension from longitude array:
        shp = self.swaths[self.component]['Geolocation Fields']['Longitude']['data'].shape
        if len(shp) != 2 :
            self.logger.error( 'expected 2D longitude array, found shape : %s' % str(shp) )
            raise Exception
        #endif
        if shp[0] == self.ntime :
            self.np = shp[1]
        elif shp[1] == self.ntime :
            self.np = shp[0]
        else :
            self.logger.error( 'expected (time,np) or (np,time) shape for longitudes;' )
            self.logger.erorr( 'found shape %s while time is %i' % (str(shp),self.ntime) )
            raise Exception
        #endif
        
        # extract level dimension from kernel array:
        shp = self.swaths[self.component]['Data Fields']['AveragingKernel']['data'].shape
        if len(shp) != 3 :
            self.logger.error( 'expected 3D kernel array, found shape : %s' % str(shp) )
            raise Exception
        #endif
        if shp[0:2] == (self.ntime,self.np) :
            self.nlayer = shp[2]
        elif shp[1:3] == (self.ntime,self.np) :
            self.nlayer = shp[0]
        else :
            self.logger.error( 'expected (time,np,nlayer) or (nlayer,time,np) shape for kernel;' )
            self.logger.erorr( 'found shape %s while (time,np) is (%i,%i)' % (str(shp),self.ntime,self.np) )
            raise Exception
        #endif
        
        # loop over all data fields to set dimension names:
        for field in self.swaths[component].keys() :
            # skip some ...
            if field == 'attrs' : continue
            # variables in field group:
            for varname in self.swaths[component][field].keys() :
                # get shape:
                shp = self.swaths[component][field][varname]['data'].shape
                # search known combinations:
                if shp == (self.ntime,) :
                    dimnames = ('time',)
                elif shp == (self.ntime,self.np) :
                    dimnames = ('time','pixel')
                elif shp == (4,self.ntime,self.np) :
                    dimnames = ('corner','time','pixel')
                elif shp == (self.ntime,self.np,self.nlayer) :
                    dimnames = ('time','pixel','layer')
                elif shp == (self.ntime,self.np,self.nlayer+1) :
                    dimnames = ('time','pixel','layer_interface')
                elif shp == (self.nlayer,self.ntime,self.np) :
                    dimnames = ('layer','time','pixel')
                elif shp == (self.nlayer,) :
                    dimnames = ('layer',)
                else :
                    self.logger.error( 'unknown shape %s for variable "%s";' % (str(shp),varname) )
                    self.logger.error( 'dimensions time %i, pixel %i, layer %i' % (self.ntime,self.np,self.nlayer) )
                    raise OmiLevelException
                #endif
                # store:
                self.swaths[component][field][varname]['dimnames'] = dimnames
            #endfor  # variables
        #endfor # field groups
        
        #
        # Adhox fix: corners have strange order, convert to mathematical default of counter-clock wise:
        #                                    
        #            _ - o 0                         _ - o 1
        #      2 o -      \         --->       2 o -      \     
        #         \    _ - o 1                    \    _ - o 0
        #        3 o -                           3 o -
        #
        # loop over corner variables:
        for vname in ['LongitudeCornerpoints','LatitudeCornerpoints'] :
            # copy original:
            values = self.swaths[self.component]['Geolocation Fields'][vname]['data'].copy()
            # swap, dimension order is (corner,scan,pix)
            self.swaths[self.component]['Geolocation Fields'][vname]['data'][0,:,:] = values[1,:,:]
            self.swaths[self.component]['Geolocation Fields'][vname]['data'][1,:,:] = values[0,:,:]
        #endfor
        
    #enddef __init__
    
    # *
    
    def ConvertTime( self ) :
    
        """
        Convert time values (one per scanline) to datetime structures (2D, for each pixel)
        
        New variable is added::
        
          self.swaths[component]['Geolocation Fields']['Datetime']
        """
    
        # modules:
        import logging
        import datetime
        import netCDF4
        import numpy
        
        #
        # Original 'Time' units and description:
        #
        #   title = "Time at Start of Scan (s, UNIX time / POSIX time), number of seconds that have elapsed since midnight Coordinated Universal Time (UTC), 1 January 1970."
        #   units = "s"
        #
        # Create new field 'Datetime' field with units:
        #   units = "Seconds since 1970-01-01 00:00'
        #
        # select:
        varid = self.swaths[self.component]['Geolocation Fields']['Time']
        # values:
        tvalues   = varid['data']
        # extract description:
        long_name = varid['long_name'].decode('latin-1')
        # check ...
        key = 'Time at Start of Scan (s, UNIX time / POSIX time), number of seconds that have elapsed since midnight Coordinated Universal Time (UTC),'
        if long_name.startswith(key) :
            # remove leading description:
            time0 = long_name.replace(key,'').replace('.','').strip()
            # extract datetime object:
            t0 = datetime.datetime.strptime(time0,'%d %B %Y')
            # convert:
            var = {}
            var['units'    ] = t0.strftime('seconds since %Y-%m-%d %H:%M:%H')
            var['long_name'] = long_name
            if 'mask' in dir(tvalues) :
                values1d = netCDF4.num2date( tvalues.data, var['units'] )
            else :
                values1d = netCDF4.num2date( tvalues     , var['units'] )
            #endif
        # alternative:
        #   "Time at Start of Scan (s, TAI93)"
        elif 'TAI' in long_name :
            # find start:
            i0 = long_name.index('TAI')
            # extract:
            year = int(long_name[i0+3:].replace(')',''))
            # convert to 4-digits if necessary:
            if year < 100 :
                if year > 50 :
                    year = 1900 + year
                else :
                    year = 2000 + year
                #endif
            #endif
            # reference time:
            t0 = datetime.datetime(year,1,1,0,0,0)
            # convert:
            var = {}
            var['units'    ] = t0.strftime('seconds since %Y-%m-%d %H:%M:%H')
            var['long_name'] = long_name
            values1d         = netCDF4.num2date( tvalues, var['units'] )
        else :
            self.logger.error( 'could not convert time units "%s"' % long_name )
            self.logger.error( 'first value : %f' % tvalues[0] )
            raise Exception
        #endif
        
        # expand to 2D:
        var['data'] = numpy.zeros( (self.ntime,self.np), values1d.dtype )
        for ip in range(self.np) :
            var['data'][:,ip] = values1d
        #endfor
        
        # set dim names:
        var['dimnames'] = ('time','pixel')
        
        # store:
        self.swaths[self.component]['Geolocation Fields']['Datetime'] = var

    #enddef ConvertTime
    
    # *
    
    def SelectPixels( self, rcf, rckey ) :
    
        """
        Apply filters specified in rcfile.
        
        Arguments:
        
        * ``rcf``       :  :py:class:`.RcFile` object with settings
        * ``rckey``     :  basename for rcfile keys, e.g. "omi" for the example below
            
        Return values:
        
        * ``selected`` : boolean array with shape of track (ntime,npix) which
          is True if a pixel passed all checks;
        * ``history``  : list of character string describing operations.
          
        Example configuration::
        
            ! Specifiy a list of filter names.
            ! For each name, specify the variable which values are used for testing,
            ! the type test, and eventually some thresholds or other settings for this type.
            ! The units of the thresholds should match with the units in the variables,
            ! the expected units have to be defined too.
            ! The examples below show possible types and their settings.

            ! filters:
            omi.filters                        :  lons lats albedo valid

            ! select range of values:
            omi.filter.lons.var                :  Geolocation Fields/Longitude
            omi.filter.lons.units              :  degrees
            omi.filter.lons.type               :  minmax
            omi.filter.lons.minmax             :  -15.0 35.0

            ! select above a minimum:
            omi.filter.lats.var                :  Geolocation Fields/Latitude
            omi.filter.lats.units              :  degrees
            omi.filter.lats.type               :  minmax
            omi.filter.lats.minmax             :  35.0 70.0

            ! select below a maximum:
            omi.filter.albedo.var              :  Data Fields/SurfaceAlbedo
            omi.filter.albedo.units            :  1
            omi.filter.albedo.type             :  max
            omi.filter.albedo.max              :  0.3

            ! select only values with data (no "_FillValue"):
            omi.filter.valid.var               :  Data Fields/NO2RetrievalTroposphericVerticalColumn
            omi.filter.valid.type              :  valid
        
        """
        
        # modules:
        import logging
        import numpy
   
        # init selection, by default accept all:
        selected = numpy.ones( (self.ntime,self.np), 'bool' )
        
        # init history:
        history = []

        # which filters ?
        filternames = rcf.get( '%s.filters' % rckey ).split()
        # loop over filters:
        for filtername in filternames :
            # filter key in rcfiles:
            fkey = '%s.filter.%s' % (rckey,filtername)
            # which variable ?
            varpath = rcf.get(fkey+'.var')
            fields,varname = varpath.split('/')
            # extract values:
            values = self.swaths[self.component][fields][varname]['data']
            # extract units:
            vunits = self.swaths[self.component][fields][varname]['units']
            # which type ?
            filtertype = rcf.get(fkey+'.type')
            # switch:
            #~ minimum value:
            if filtertype == 'equal' :
                # read value:
                svalue = rcf.get(fkey+'.value')
                fvalue = float(svalue)
                # apply filter:
                selected = selected & (values == fvalue)
                # update history:
                history.append( "selected pixels with '%s' == %s" % (varname,svalue) )
            #~ minimum value:
            elif filtertype == 'min' :
                # read threshold:
                smin = rcf.get(fkey+'.min')
                fmin = float(smin)
                # check units ...
                funits = rcf.get(fkey+'.units')
                if funits != vunits :
                    self.logger.error( 'variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
                    raise Exception
                #endif
                # apply filter:
                selected = selected & (values >= fmin)
                # update history:
                history.append( "selected pixels with '%s' >= %s" % (varname,smin) )
            #~ maximum value:
            elif filtertype == 'max' :
                # read threshold:
                smax = rcf.get(fkey+'.max')
                fmax = float(smax)
                # check units ...
                funits = rcf.get(fkey+'.units')
                if funits != vunits :
                    self.logger.error( 'variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
                    self.logger.error( 'value range : %f - %f' % (numpy.nanmin(values),numpy.nanmax(values)) )
                    raise Exception
                #endif
                # apply filter:
                selected = selected & (values <= fmax)
                # update history:
                history.append( "selected pixels with '%s' <= %s" % (varname,smax) )
            #~ range of valid values:
            elif filtertype == 'minmax' :
                # read thresholds:
                smin,smax = rcf.get(fkey+'.minmax').split()
                fmin = float(smin)
                fmax = float(smax)
                # check units ...
                funits = rcf.get(fkey+'.units')
                if funits != vunits :
                    self.logger.error( 'variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
                    raise Exception
                #endif
                # apply filter:
                selected = selected & (values >= fmin) & (values <= fmax)
                # update history:
                history.append( "selected pixels with '%s' in [%s,%s]" % (varname,smin,smax) )
            #~ only unmasked values:
            elif filtername == 'valid' :
                # only accept non-masked fields:
                selected = selected & (~ values.mask)
                # update history:
                history.append( "selected pixels with valid '%s' values" % varname )
            #~ not yet ...
            else :
                self.logger.error( 'unsupported filter name : %s' % filtername )
                raise Exception
            #endif
        #endfor
        
        # if ncecessary, convert from a masked to a normal array:
        if 'mask' in dir(selected) : selected = selected.filled(fill_value=False)

        # ok
        return selected,history
    
    #enddef SelectPixels
        
        

#endclass OMI_HDF5_File


########################################################################
###
### test
###
########################################################################

# runned as main program ?
if __name__ == "__main__" :

    # test file:
    fname = '/modas/scratch02/projects/ISOTROP/synth-sat-obs/sample/L2-files/S5_LEO/ISOTROP-S5_L2-L2-NO2_2003m0601t0105-o00014_v001-2013m0320t110311.he5'

    # read:
    omi = OMI_HDF5_File( fname, 'NO2', verbose=True )


#endif  # main program


########################################################################
###
### end
###
########################################################################


