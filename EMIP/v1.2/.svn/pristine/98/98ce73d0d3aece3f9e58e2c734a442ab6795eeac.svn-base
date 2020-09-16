########################################################################
###
### mapping tools
###
########################################################################

class Mapping_1D_to_2D( object ) :
    
    """
    Storage for mapping arrays.
    
    If filename is provided, read content from file.
    Otherwise provide ``shp=(ny,nx)`` to initialize as empty.
    
    
    """
    
    def __init__( self, filename=None, shp=None ) :
    
        """
        Initialize mapping for source shape.
        """
        
        # modules:
        import numpy
        
        # read?
        if filename is not None :
        
            # read:
            self.Read( filename )
            
        elif shp is not None :
        
            # store:
            self.ny,self.nx = shp

            # storage for target index lists:
            self.ii = numpy.array([],dtype='i4')
            self.jj = numpy.array([],dtype='i4')

            # storage for source index lists:
            self.pp = numpy.array([],dtype='i4')
            # storage for source weights:
            self.pw = numpy.array([],dtype='f4')

            # storage for mapping weights:
            self.ww = numpy.array([],dtype='f4')
            self.units = '1'

            # init weights sum:
            self.wsum = numpy.ma.array( data=numpy.zeros((self.ny,self.nx),dtype='f4'),
                                        mask=numpy.ones ((self.ny,self.nx),bool ) )
                                        
        else :
            print( 'ERROR - either provide "filename" or "shp"' )
            raise Exception
        #endif
        
    #enddef __init__
    
    # *
    
    def __len__( self ) :
    
        """
        Return number of source values.
        """
        
        # ok
        return len(self.ww)
        
    #enddef __len__
    
    # *
    
    def Append( self, ip, pw, jj, ii, ww, units ) :
    
        """
        Add distribution of source value 'ip' over target cells.
        Integer lists jj and ii define target cells, 
        float list ww the weights of the source values in the cell.
        """
        
        # modules:
        import numpy
        
        # add:
        self.pp = numpy.append( self.pp, numpy.array( [ip] * len(ww), dtype='i4' ) )
        self.pw = numpy.append( self.pw, numpy.array( [pw] * len(ww), dtype='f4' ) )
        self.ii = numpy.append( self.ii, numpy.array( ii, dtype='i4' ) )
        self.jj = numpy.append( self.jj, numpy.array( jj, dtype='i4' ) )
        self.ww = numpy.append( self.ww, numpy.array( ww, dtype='f4' ) )
        self.units = units
        
        # loop over contributions:
        for j,i,w in zip(jj,ii,ww) :
            # update weight sum in cells:
            self.wsum.data[j,i] += w
            # reset flags:
            self.wsum.mask[j,i] = False
        #endfor
        
    #enddef Append
    
    # *
    
    def Apply( self, source ) :
    
        """
        Apply mapping to source array.
        First dimension is supposed to be the source index.
        
        Return values:
        
        * ``values`` : weighted combinations of source values
        """
        
        # modules:
        import numpy
        
        
        # cells with contribution:
        jj,ii = self.GetIndices()

        # rank:
        rank = len(source.shape)
        # switch:
        if rank == 1 :
        
            # target field:
            field = numpy.zeros((self.ny,self.nx),float)
            # loop over source values:
            for k in range(len(self.ww)) :
                # add:
                field[self.jj[k],self.ii[k]] += self.ww[k] * source[self.pp[k]]
            #endfor

            # extract and normalize:
            values = field[jj,ii] / self.wsum[jj,ii]
        
        elif rank == 2 :
        
            # number of levels:
            nz = source.shape[1]
            
            # target field:
            field = numpy.zeros((self.ny,self.nx,nz),float)
            # loop over source values:
            for k in range(len(self.ww)) :
                # add:
                field[self.jj[k],self.ii[k],:] += self.ww[k] * source[self.pp[k],:]
            #endfor

            # extract and normalize:
            values = numpy.zeros((len(ii),nz),float)
            for iz in range(nz) :
                values[:,iz] = field[jj,ii,iz] / self.wsum[jj,ii]
            #endfor
            
        else :
        
            print( 'ERROR - unsupported source shape "%s"' % str(source.shape) )
            raise Exception
        
        #endif

        # ok
        return values
        
    #enddef Apply
    
    # *
    
    def GetIndices( self ) :
    
        """
        Index arrays in target grid where weight sum exceeds zero.
        
        Return values:

        * ``jj`, ``ii`` : indices on target grid
        """
        
        # modules:
        import numpy
        
        # ok
        return numpy.where( self.wsum > 0.0 )
        
    #enddef GetIndices
    
    # *
    
    def GetCompression( self, name, dimnames ) :
    
        """
        Return compression coordinate and related attributes.
        
        Arguments:
        
        * ``name`` : coordinate name, for example ``pixel``
        * ``dimanames`` : names of y and x coordinates, e.g. ``('lat','lon')``
        
        Return values:
        
        * ``index`` : 1D index of cells with non-zero weights
        * ``attrs`` : dictionairy with attributes
        
        """
        
        # modules:
        import numpy
        
        # attribute:
        compress = '%s %s' % dimnames
        # attributet:
        description = "original zero-based indices: `Y` = `P`/len(`X`), `X` = `P` mod len(`X`)"
        description = description.replace('X',dimnames[1])
        description = description.replace('Y',dimnames[0])
        description = description.replace('P',name)
        # collect:
        attrs = {}
        attrs['units'] = '1'
        attrs['long_name'] = 'compressed coordinate'
        attrs['compress'] = compress
        attrs['description'] = description

        # non-zero weights:
        jj,ii = self.GetIndices()
        # numbered index: ...
        indices = numpy.array( jj * self.nx + ii, dtype='i4' )
        
        # ok
        return indices,attrs
        
    #enddef GetCompression
                
    
    # *
    
    def Write( self, filename, attrs={} ) :
    
        """
        Write mapping weights to netcdf file.
        """
        
        # modules:
        import os
        import netCDF4
        
        # target directory:
        dirname = os.path.dirname(filename)
        # create if necessary:
        if not os.path.isdir(dirname) : os.makedirs( dirname )
        
        # total number of mappings:
        n = len(self.ww)
        
        # create:
        ncid = netCDF4.Dataset( filename, 'w' )
        
        # attributes:
        for key in attrs.keys() : ncid.setncattr( key, attrs[key] )

        # add dimensions:
        ncid.createDimension( 'x', self.nx )
        ncid.createDimension( 'y', self.ny )
        ncid.createDimension( 'n', n )

        # write list:
        varid = ncid.createVariable( 'ii', self.ii.dtype, ('n',) )
        varid.setncattr( 'long_name', 'x-index in target array (zero based)' )
        varid.setncattr( 'units', '1' )
        varid[:] = self.ii
        # write list:
        varid = ncid.createVariable( 'jj', self.jj.dtype, ('n',) )
        varid.setncattr( 'long_name', 'y-index in target array (zero based)' )
        varid.setncattr( 'units', '1' )
        varid[:] = self.jj
        # write list:
        varid = ncid.createVariable( 'pp', self.pp.dtype, ('n',) )
        varid.setncattr( 'long_name', 'index in source array (zero based)' )
        varid.setncattr( 'units', '1' )
        varid[:] = self.pp
        # write list:
        varid = ncid.createVariable( 'pw', self.pw.dtype, ('n',) )
        varid.setncattr( 'long_name', 'weight of source value' )
        varid.setncattr( 'units', self.units )
        varid[:] = self.pw
        # write list:
        varid = ncid.createVariable( 'ww', self.ww.dtype, ('n',) )
        varid.setncattr( 'long_name', 'weights of source values' )
        varid.setncattr( 'units', self.units )
        varid[:] = self.ww
        
        # write total weight:
        varid = ncid.createVariable( 'wsum', self.wsum.dtype, ('y','x'), fill_value=-999.9 )
        varid.setncattr( 'long_name', 'weight sums' )
        varid.setncattr( 'units', self.units )
        varid[:] = self.wsum
        
        # close:
        ncid.close()
        
    #enddef Write
    
    # *
    
    def Read( self, filename ) :
    
        """
        Read mapping weights from netcdf file.
        """
        
        # modules:
        import os
        import netCDF4
        
        # check ...
        if not os.path.isfile(filename) :
            self.logger.error( 'file not found: %s' % filename )
            raise Exception
        #endif
        
        # open:
        ncid = netCDF4.Dataset( filename, 'r' )
        
        # dims:
        self.nx = len(ncid.dimensions['x'])
        self.ny = len(ncid.dimensions['y'])
        
        # read lists:
        self.ii = ncid.variables['ii'][:]
        self.jj = ncid.variables['jj'][:]
        self.pp = ncid.variables['pp'][:]
        self.pw = ncid.variables['pw'][:]
        self.ww = ncid.variables['ww'][:]
        
        # read map with total weights:
        self.wsum = ncid.variables['wsum'][:]
        
        # close:
        ncid.close()
        
    #enddef Read
    
#endclass Mapping_1D_to_2D
        

########################################################################
###
### end
###
########################################################################

