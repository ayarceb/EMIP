
"""
******************
``TROPOMI`` module
******************

Access TROPOMI data in NetCDF4 as available on TEMIS site.


Class hierchy
=============

The classes are defined according to the following hierchy:

* :py:class:`.TROPOMI_NC4_File`


Classes
=======

"""

########################################################################
###
### OMI data class
###
########################################################################

class TROPOMI_NC4_File( object ) :

    """
    Base class to store content of TROPOMI NetCDF4 file.
    
    Most relevant variables::

        PRODUCT/
            dimensions:
                scanline = 3245 ;
                ground_pixel = 450 ;
                corner = 4 ;
                time = 1 ;
                layer = 34 ;
            variables:
                int time(time) ;
                float latitude (time, scanline, ground_pixel) ;
  		                units = "degrees_east" ;
                float longitude(time, scanline, ground_pixel) ;
  		                units = "degrees_north" ;
                float nitrogendioxide_tropospheric_column(time, scanline, ground_pixel) ;
  		                units = "mol m-2" ;
  		                multiplication_factor_to_convert_to_molecules_percm2 = 6.02214e+19f ;
  	            float nitrogendioxide_tropospheric_column_precision(time, scanline, ground_pixel) ;
  		                units = "mol m-2" ;
  		                multiplication_factor_to_convert_to_molecules_percm2 = 6.022141e+19f ;
          	    float averaging_kernel(time, scanline, ground_pixel, layer) ;
  		                units = "1" ;
  	            ubyte qa_value(time, scanline, ground_pixel) ;
  		                units = "1" ;
  		                long_name = "data quality value" ;
  		                comment = "A continuous quality descriptor, varying between 0 (no data) and 1 (full quality data). Recommend to ignore data with qa_value < 0.5" ;
            SUPPORT_DATA/
                GEOLOCATIONS/
                    variables:
                        float solar_zenith_angle(time, scanline, ground_pixel) ;
      		                    units = "degree" ;
      	                float latitude_bounds   (time, scanline, ground_pixel, corner) ;
  		                        units = "degrees_east" ;
      	                float longitude_bounds  (time, scanline, ground_pixel, corner) ;
  		                        units = "degrees_north" ;
                INPUT_DATA/
                    variables:
      	                float surface_altitude(time, scanline, ground_pixel) ;
      		                    units = "m" ;
      	                float surface_pressure(time, scanline, ground_pixel) ;
      		                    units = "Pa" ;
      	                float surface_albedo(time, scanline, ground_pixel) ;
      		                    units = "1" ;
      	                float cloud_fraction_crb(time, scanline, ground_pixel) ;
      		                    units = "1" ;

    Arguments:
    
    * ``filename``    : name of input file

    """
    
    def __init__( self, filename, verbose=False ) :
    
        """
        Read OMI data file.        
        """
        
        # modules:
        import os
        import netCDF4
        import numpy
        
        # store:
        self.filename = filename

        # check ...
        if not os.path.isfile(filename) :
            print( 'ERROR - file not found : %s' % filename )
            raise Exception
        #endif
        
        # open:
        ncid = netCDF4.Dataset( filename, 'r' )
        
        # relevant stuff in this group:
        group = 'PRODUCT'
        
        # copy some dimensions for convenience:
        gid = ncid.groups[group]
        self.ntime = len(gid.dimensions['time'])
        self.nscan = len(gid.dimensions['scanline'])
        self.ngpix = len(gid.dimensions['ground_pixel'])
        self.nlayer= len(gid.dimensions['layer'])

        # init storage:
        self.groups = {}
        ## info ...
        #print( 'xxx group "%s"' % group )
        # recursive fill:
        self.groups[group] = self._ReadGroup( ncid.groups[group], indent='  ' )
        
        # close:
        ncid.close()
        
    #enddef __init__
    
    # *
    
    def _ReadGroup( self, gid, indent='' ) :
    
        """
        Return group content in dictionairy.
        """
        
        # init result:
        data = {}
        # loop over variables:
        for varname in gid.variables.keys() :
           
            ## info ...
            #print( indent+'xxx var "%s"' % varname )

            # create storage:
            data[varname] = {}
            # access:
            varid = gid.variables[varname]
            
            # copy dimension names:
            data[varname]['dimnames'] = varid.dimensions
            
            # read data:
            data[varname]['data'] = varid[:]
            
            # copy attributes:
            data[varname]['attrs'] = {}
            for key in varid.ncattrs() :
                data[varname]['attrs'][key] = varid.getncattr(key)
            #endfor
            
            # add some missing units ...
            if 'units' not in data[varname]['attrs'].keys() :
                # assume something ...
                if 'longitude' in varname :
                    data[varname]['attrs']['units'] = 'degrees_east'
                elif 'latitude' in varname :
                    data[varname]['attrs']['units'] = 'degrees_north'
                #else :
                #    print( 'ERROR - could not guess units for variable "%s"' % varname )
                #    raise Exception
                #endif
            #endif # no units
            
            # add missing long name ...
            if 'long_name' not in data[varname]['attrs'].keys() :
                data[varname]['attrs']['long_name'] = varname.replace('_',' ')
            #endif
            
        #endfor # variables
        
        # loop over sub-groups:
        for group in gid.groups.keys() :
            ## info ...
            #print( indent+'xxx group "%s"' % group )
            # recursive call:
            data[group] = self._ReadGroup( gid.groups[group], indent=indent+'  ' )
        #endfor

        # ok
        return data

    #enddef _ReadGroup
    
    # *
    
    def GetVar( self, path ) :
    
        """
        Return variable dictionairy for path::
        
           var['data'][:]   # numpy (masked) array
           var['attrs'] = { 'units' : '1', ... }
        """
        
        # recursive call:
        return self._GetVar( self.groups, path, '' )
        
    #enddef GetVar
    
    # *
    
    def _GetVar( self, groups, path, root ) :
    
        """
        Extract variable data.
        """
        
        # split:
        if '/' in path :
            # split:
            group,rpath = path.split('/',1)
            # check ...
            if group not in groups.keys() :
                print( 'ERROR - group "%s" not found on root "%s"' % (group,root) )
                print( 'ERROR - filename: %s' % self.filename )
                raise Exception
            #endif
            # recursive call:
            var = self._GetVar( groups[group], rpath, root+'/'+group )
        else :
            # check ...
            if path not in groups.keys() :
                print( 'ERROR - variable "%s" not found on root "%s"' % (path,root) )
                print( 'ERROR - filename: %s' % self.filename )
                raise Exception
            #endif
            # copy:
            var = groups[path]
        #endif
        
        # ok
        return var
        
    #enddef _GetVar

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
            omi.filter.lons.var                :  PRODUCT/longitude
            omi.filter.lons.units              :  degrees_east
            omi.filter.lons.type               :  minmax
            omi.filter.lons.minmax             :  -15.0 35.0

            ! select above a minimum:
            omi.filter.lats.var                :  PRODUCT/latitude
            omi.filter.lats.units              :  degrees_north
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
        import numpy
   
        # init selection, by default accept all:
        selected = numpy.ones( (self.nscan,self.ngpix), 'bool' )
        
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
            group,varname = varpath.split('/')
            # extract values:
            values = self.groups[group][varname]['data']
            # extract units:
            vunits = self.groups[group][varname]['attrs']['units']
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
                    print( 'ERROR - variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
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
                    print( 'ERROR - variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
                    print( 'ERROR - value range : %f - %f' % (numpy.nanmin(values),numpy.nanmax(values)) )
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
                    print( 'ERROR - variable "%s" has units "%s" while threshold has units "%s"' % (varname,vunits,funits) )
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
                print( 'ERROR - unsupported filter name : %s' % filtername )
                raise Exception
            #endif
        #endfor # filters
        
        # if ncecessary, convert from a masked to a normal array:
        if 'mask' in dir(selected) : selected = selected.filled(fill_value=False)

        # ok
        return selected,history
    
    #enddef SelectPixels
        
        

#endclass TROPOMI_NC4_File


########################################################################
###
### end
###
########################################################################


