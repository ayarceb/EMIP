"""
Tools for Carthesian-Grid .

A Carthesian-Grid is defined by two 1D axes 
for longitudes and latitudes with a regular spacing.
The assumption on regular spacing is used to speedup computations.
"""

class CarthesianGrid( object ) :

    """
    Regular lon/lat grid cells.
    Assumes constant spacing in lon and in lat direction.
    
    Set polygons flag to define attribute 'pg[j,i]' with cell polygons.
    """

    def __init__( self, west=None, dlon=None, nlon=None, lons=None,
                        south=None, dlat=None, nlat=None, lats=None,
                        center=False, polygons=False ) :

        """
        Setup lon/lat grid cell.
        By default, west and south define the edges of grid cells;
        if center is set to True they define the grid box center however.
        """

        # modules:
        import numpy
        
        # tools
        import grid_tools

        # grid points or cells ?
        if lons is None :
            if center :
                self.west = west - dlon/2.0
            else :
                self.west = west
            #endif
            self.dlon = dlon
            self.nlon = nlon
            self.lons = self.west + (numpy.arange(self.nlon)+0.5)*self.dlon
        else :
            self.lons = lons
            self.dlon = lons[1]-lons[0]
            self.nlon = len(lons)
            self.west = lons[0] - self.dlon*0.5
        #endif
        self.east = self.west + self.nlon*self.dlon

        if lats is None :
            if center :
                self.south = south - dlon/2.0
            else :
                self.south = south
            #endif
            self.dlat  = dlat
            self.nlat  = nlat
            self.lats  = self.south + (numpy.arange(self.nlat)+0.5)*self.dlat
        else :
            self.lats  = lats
            self.dlat  = lats[1]-lats[0]
            self.nlat  = len(lats)
            self.south = lats[0] - self.dlat*0.5
        #endif
        self.north = self.south + self.nlat*self.dlat

        # number of cells per degree:
        self.nlon_per_degree = int(round(1.0/self.dlon))
        self.nlat_per_degree = int(round(1.0/self.dlat))

        # store domain bounding box:
        self.domain = [self.west,self.east,self.south,self.north]

        # grid cell boundaries:
        self.blons = numpy.zeros((self.nlon+1))
        self.blons[0] = self.west
        for i in range(self.nlon) : self.blons[i+1] = self.lons[i] + 0.5*self.dlon
        # grid cell boundaries:
        self.blats = numpy.zeros((self.nlat+1))
        self.blats[0] = self.south
        for j in range(self.nlat) : self.blats[j+1] = self.lats[j] + 0.5*self.dlat

        # grid cell boundaries:
        self.lon_bounds = numpy.zeros((2,self.nlon))
        self.lon_bounds[0,:] = self.lons - 0.5*self.dlon
        self.lon_bounds[1,:] = self.lons + 0.5*self.dlon
        # grid cell boundaries:
        self.lat_bounds = numpy.zeros((2,self.nlat))
        self.lat_bounds[0,:] = self.lats - 0.5*self.dlat
        self.lat_bounds[1,:] = self.lats + 0.5*self.dlat

        # 2D corner fields:
        self.xx,self.yy = numpy.meshgrid( self.blons, self.blats )
        
        # 2D mid fields:
        self.xxm,self.yym = numpy.meshgrid( self.lons, self.lats )

        # cell area in m2:
        self.area = numpy.zeros((self.nlat,self.nlon),float)
        for j in range(self.nlat) :
            for i in range(self.nlon) :
                self.area[j,i] = grid_tools.ll_area( 
                                          self.lon_bounds[0,i],   # west
                                          self.lon_bounds[1,i],   # east
                                          self.lat_bounds[0,j],   # south
                                          self.lat_bounds[1,j]  ) # north
            #endfor
        #endfor
        
        #
        # polygons
        #
        
        # create?
        if polygons :
            # create storage:
            self.pg = []
            # loop lat rows:
            for j in range(self.nlat) :
                # initialze empty row:
                row = []
                # loop over lon cells:
                for i in range(self.nlon) :
                    # define polygon:
                    row.append( self.GetPolygon( j, i ) )
                #endfor # i
                # store:
                self.pg.append( row )
            #endfor # j
        else :
            # not defined:
            self.pg = None
        #endif # define polygons

    #enddef __init__
    
    # *
    
    def GetPolygon( self, j, i ) :
    
        """
        Return polygon for grid cell [j,i].
        """

        # tools:
        import go

        # sides:
        west  = self.lon_bounds[0,i]
        east  = self.lon_bounds[1,i]
        south = self.lat_bounds[0,j]
        north = self.lat_bounds[1,j]
        # loop over corners (counter clock wise):
        #    2 o--o 1
        #      |  | 
        #    3 o--o 0
        corners = []
        corners.append( go.vector.Vector( east, south ) )
        corners.append( go.vector.Vector( east, north ) )
        corners.append( go.vector.Vector( west, north ) )
        corners.append( go.vector.Vector( west, south ) )

        # define polygon:
        return go.vector.Polygon( corners=corners )
        
    #enddef GetPolygon
    
    # *
    
    def GetCorners( self ) :
    
        """
        Return corner arrays:
        
        * ``cxx`` : longitude corners shaped ``(nlat,nlon,4)``
        * ``cyy`` : idem for latitude
        
        Corners are ordered counter-clock wise::
        
            2 o--o 1
              |  |
            3 o--o 0
        """
        
        # modules:
        import numpy
        
        # storage:
        cxx = numpy.zeros((self.nlat,self.nlon,4),float)
        cyy = numpy.zeros((self.nlat,self.nlon,4),float)
        
        # fill longitudes:
        cxx[:,:,0] = self.xxm + 0.5 * self.dlon
        cxx[:,:,1] = self.xxm + 0.5 * self.dlon
        cxx[:,:,2] = self.xxm - 0.5 * self.dlon
        cxx[:,:,3] = self.xxm - 0.5 * self.dlon
        
        # fill latitudes:
        cyy[:,:,0] = self.yym - 0.5 * self.dlat
        cyy[:,:,1] = self.yym + 0.5 * self.dlat
        cyy[:,:,2] = self.yym + 0.5 * self.dlat
        cyy[:,:,3] = self.yym - 0.5 * self.dlat
        
        # ok:
        return cxx,cyy
        
    #enddef GetCorners
        
    
    # *

    def in_domain( self, lat, lon ) :

        """
        Returns True if point is in domain.
        """

        # check location:
        indom = (lon >= self.west ) and (lon <= self.east ) and \
                (lat >= self.south) and (lat <= self.north)

        # ok
        return indom

    #enddef

    # *

    def sample_indices( self, lat, lon, missing_value=None ) :

        """
        Returns (j,i) with cell indices to sample (lat,lon) .
        """

        # modules:
        import numpy

        # outside region:

        # i-index
        if (lon < self.west ) or (lon > self.east ) :
            i = missing_value
        else :
            i = min( int(numpy.floor( (lon - self.west)/self.dlon )), self.nlon-1 )
        #endif

        # j-index
        if (lat < self.south) or (lat > self.north) :
            j = missing_value
        else :
            j = min( int(numpy.floor( (lat - self.south)/self.dlat )), self.nlat-1 )
        #endif

        # ok
        return j,i

    #enddef sample_indices
    
    # *
    
    def interp_indices( self, lat, lon ) :
    
        """
        Return indices and weights for bi-linear interpolation.

        Return values:
          ww  : weights
          jj  : lat indices
          ii  : lon indices
          
        All values are negative if location is outside domain.
        """
        
        # init results; indices of surrounding points:
        #
        #  lat direction
        #   ^ 
        #   |  2    3
        #   |    *
        #   |  0    1
        #   +----------> lon direction
        #
        ii = [-999,-999,-999,-999]
        jj = [-999,-999,-999,-999]
        ww = [-999.9,-999.9,-999.9,-999.9]

        # indices of source cell including source point:
        j,i = self.sample_indices( lat, lon, missing_value=-999 )
        #
        if (i < 0) or (j < 0) :
        
            # outside domain ...
            pass
            #print 'ERROR location ',lat,lon, ' outside domain ', self.domain
            #raise Exception
            
        else :

            # lon direction
            if lon < self.lons[i] :
                alfa = ( lon - self.lons[i-1] ) / (self.lons[i] - self.lons[i-1] )
                ii[0] = i - 1 ; ww[0] = 1.0 - alfa
                ii[1] = i     ; ww[1] = alfa
                ii[2] = i - 1 ; ww[2] = 1.0 - alfa
                ii[3] = i     ; ww[3] = alfa
            else :
                alfa = ( lon - self.lons[i] ) / (self.lons[i+1] - self.lons[i] )
                ii[0] = i     ; ww[0] = 1.0 - alfa
                ii[1] = i + 1 ; ww[1] = alfa
                ii[2] = i     ; ww[2] = 1.0 - alfa
                ii[3] = i + 1 ; ww[3] = alfa
            #endif
            #
            # lat direction
            if lat < self.lats[j] :
                alfa = ( lat - self.lats[j-1] ) / (self.lats[j] - self.lats[j-1] )
                jj[0] = j - 1 ; ww[0] = ww[0] * (1.0 - alfa)
                jj[1] = j - 1 ; ww[1] = ww[1] * (1.0 - alfa)
                jj[2] = j     ; ww[2] = ww[2] * alfa
                jj[3] = j     ; ww[3] = ww[3] * alfa
            else :
                alfa = ( lat - self.lats[j] ) / (self.lats[j+1] - self.lats[j] )
                jj[0] = j     ; ww[0] = ww[0] * (1.0 - alfa)
                jj[1] = j     ; ww[1] = ww[1] * (1.0 - alfa)
                jj[2] = j + 1 ; ww[2] = ww[2] * alfa
                jj[3] = j + 1 ; ww[3] = ww[3] * alfa
            #endif
            
        #endif
        
        # ok
        return ww,jj,ii

    #enddef interp_indices
    
    # * overlap
    
    def SharedDomain( self, cg=None, bbox=None ) :
    
        """
        Return domain (west,east,south,north) that overlaps
        with the domain defined by either the :py:class:`CarthesianGrid` object ``cg``
        or the bounding box ``bbox=(west,easth,south,north)``.
        In case the domains do not overlap with this grid, ``None`` is returned.
        """
        
        # extract own bounding box:
        west1,east1,south1,north1 = self.domain
        # second bounding box:
        if cg is not None :
            west2,east2,south2,north2 = cg.domain
        elif bbox is not None :
            west2,east2,south2,north2 = bbox
        else :
            print( 'ERROR - provide eigher `cg` or `bbox` argument' )
            raise Exception
        #endif
        
        # no overlap?
        if (west2  > east1 ) or (east2  < west1 ) or \
           (south2 > north2) or (north2 < south1) :
            # no overlap:
            domain = None
        else :
            # shared:
            west  = max(west1,west2)
            east  = min(east1,east2)
            south = max(south1,south2)
            north = min(north1,north2)
            # fill:
            domain = (west,east,south,north)
        #endif
        
        # ok
        return domain
        
    #enddef SharedDomain
    
    # *
    
    def SharedRange( self, cg=None, bbox=None ) :
    
        """
        Returns lat/lon index ranges for overlap with square domain,
        either defined by other grid (cg) or by bounding box (west,east,south,north)::
        
           j0,j1,i0,i1 = self.SharedRange( cg=cg )
           field[j0:j1,i0:i1] = ...
        
        If the grid/box do not overlap the result will be ``None``.
        """
        
        # modules:
        import numpy
        
        # shared domain:
        domain = self.SharedDomain( cg=cg, bbox=bbox )
        # defined?
        if domain is None :
            # dummy:
            result = None
        else :
            # expand:
            west,east,south,north = domain
        
            # range of cells including west and east:
            i0 = numpy.argmin( abs( self.lons - west ) )
            i1 = numpy.argmin( abs( self.lons - east ) ) + 1

            # range of cells including south and north:
            j0 = numpy.argmin( abs( self.lats - south ) )
            j1 = numpy.argmin( abs( self.lats - north ) ) + 1
            
            # fill:
            result = j0,j1,i0,i1
            
        #endif
        
        # ok:
        return result
        
    #enddef SharedRange
    
    # *
    
    def SharedGrid( self, cg ) :
    
        """
        Return grid definition of shared domain.
        """
        
        # range:
        j0,j1,i0,i1 = self.SharedRange( cg )
        # define:
        return CarthesianGrid( lons=self.lons[i0:i1], lats=self.lats[j0:j1] )
        
    #enddef SharedGrid
    
    # *
    
    def PolygonCoverage( self, pg, debug=False ) :
    
        """
        Return indices of grid cells covered by polygon,
        and fractions of the polygon area covering each of these cells.
        
        Arguments:
        
        * `pg` : :py:class:`go.vector.Polygon` object
        
        Return values:
        
        * `pga` : polygion area
        * `ii` : list of longitude cell indices
        * `jj` : list of latitude cell indices
        * `aa` : list of cell area covered by polygon
        * `units` : area units
        
        """
        
        # modules:
        import numpy
        
        # init results:
        jj,ii,aa = [],[],[]
        area_units = 'km2'

        # area of polygon, convert from m2 to km2:
        pg_area = pg.LonLat_Area() / 1e6
        
        # index range of cells potentially covered by polygon:
        srange = self.SharedRange( bbox=pg.BoundingBox() )
        # defined?
        if srange is not None :
        
            # extract:
            j0,j1,i0,i1 = srange

            # testing ...
            if debug :
                # modules:
                import matplotlib.pyplot as plt
                # new figure:
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                # polygon in red:
                pg.PlotEdges( ax, color='red', linestyle='-' )
                # bounding box:
                west,east,south,north = pg.BoundingBox()
                ax.plot( [west,east,east,west,west], [south,south,north,north,south], 'b-' )
            #endif # debug

            # loop over cells:
            for j in range(j0,j1) :
                for i in range(i0,i1) :

                    # get polygon of this cell:
                    if self.pg is None :
                        cell = self.GetPolygon( j, i )
                    else :
                        cell = self.pg[j][i]
                    #endif

                    # testing ..
                    if debug :
                        # show cell boundaries in gray:
                        cell.PlotEdges( ax, color='0.5', linestyle='-' )
                    #endif

                    # intersection with grid cell;
                    # might fail if edge of footprint is almost parrallel to to cell:
                    try :
                        segment = cell.Intersection( pg )
                    except :
                        #cell.PlotEdges( ax, color='purple', linestyle='-' )
                        #raise
                        segment = None
                    #endtry
                    # defined?
                    if segment is not None :
                        # testing ..
                        if debug :
                            segment.PlotFill( ax, color='red', alpha=0.5  )
                        #endif
                        # area of segment covering cell, convert from m2 to km2:
                        segment_area = segment.LonLat_Area() / 1e6
                        # store:
                        jj.append( j )
                        ii.append( i )
                        aa.append( segment_area )
                    #endif # segment present

                    #break

                #endfor # i

                #break
            #endfor # j

            # testing ...
            if debug :
                # show:
                plt.show()
            #endif # debug
            
        #endif  # overlap

        # ok
        return pg_area,numpy.array(jj),numpy.array(ii),numpy.array(aa),area_units
    
    #enddef PolygonCoverage

#endclass


#
# Grid - horizontal sampling and interpolation
#


class hori_sample_points( object ) :

    """
    Object to store indices to sample at list of locations.
    """

    # *

    def __init__( self, plats, plons, lli_in ) :

        """
        Pre-compute indices to sample points defined by lists
        lats and lons (note order!).
        """

        # modules:
        import logging

        # check ...
        if len(plons) != len(plats) :
            logging.error( 'number of lons and lats different : %i %i' % (len(plons),len(plats)) )
            raise ValueError
        #endif

        # store input:
        self.nlat_in = lli_in.nlat
        self.nlon_in = lli_in.nlon
        self.np = len(plons)
        #self.plons = plons
        #self.plats = plats

        # empty indices:
        self.ii_in = []
        self.jj_in = []
        # loop over locations:
        for ip in range(self.np) :
            # get indices:
            j,i = lli_in.sample_indices( plats[ip], plons[ip] )
            # store:
            self.ii_in.append(i)
            self.jj_in.append(j)
        #endfor

    #enddef  # __init__

    # *

    def apply( self, data_in, missing_value=-999 ) :

        """
        Sample data at cell locations stored in jj and ii.
        """

        # modules:
        import logging
        import numpy

        # input rank:
        ndim = len(data_in.shape)
        # check:
        if ndim < 2 :
            logging.error( 'to apply ll sampling array should have at least rank 2' )
            raise ValueError
        #endif
        # check ...
        if (data_in.shape[ndim-1] != self.nlon_in) or \
           (data_in.shape[ndim-2] != self.nlat_in) :
            logging.error( 'shape of input (%i,%i) does not match with sampling (%i,%i)' % \
                      (data_in.shape[ndim-2],data_in.shape[ndim-1],
                       self.nlat_in,self.nlon_in) )
            raise ValueError
        #endif

        # per dimension:
        if ndim == 2 :
            # original shape:
            ny,nx = data_in.shape
            # output array:
            data = numpy.ma.array( numpy.zeros((1,self.np)), dtype=float,
                                   fill_value=missing_value,
                                   mask=True)
            # loop over target cells:
            for ip in range(self.np) :
                i_in = self.ii_in[ip]
                j_in = self.jj_in[ip]
                if (i_in is not None) and (j_in is not None) :
                    data.data[0,ip] = data_in[j_in,i_in]
                    data.mask[0,ip] = False
                #endif
            #endfor
        elif ndim == 3 :
            # original shape:
            nt,ny,nx = data_in.shape
            # output array:
            data = numpy.ma.array( numpy.zeros((nt,1,self.np)), dtype=float,
                                   fill_value=missing_value,
                                   mask=True)
            # loop over target cells:
            for ip in range(self.np) :
                i_in = self.ii_in[ip]
                j_in = self.jj_in[ip]
                if (i_in is not None) and (j_in is not None) :
                    data.data[:,0,ip] = data_in[:,j_in,i_in]
                    data.mask[:,0,ip] = False
                #endif
            #endfor
        elif ndim == 4 :
            # original shape:
            nt,nz,ny,nx = data_in.shape
            # output array:
            data = numpy.ma.array( numpy.zeros((nt,nz,1,self.np)), dtype=float,
                                   fill_value=missing_value,
                                   mask=True)
            # loop over target cells:
            for ip in range(self.np) :
                i_in = self.ii_in[ip]
                j_in = self.jj_in[ip]
                if (i_in is not None) and (j_in is not None) :
                    data.data[:,:,0,ip] = data_in[:,:,j_in,i_in]
                    data.mask[:,:,0,ip] = False
                #endif
            #endfor
        else :
            logging.error( 'rank not supported yet : %i' % ndim )
            raise ValueError
        #endif

        # ok
        return data

    #endif

#endclass   # hori_sample_points



# **********************************************************************


class hori_sample_grid( object ) :

    """
    Object to store indices to sample from grid cells into grid points.
    """

    # *

    def __init__( self, lli, lli_in ) :

        """
        Pre-compute indices to sample an array defined on cells lli_in
        into points defind on lli.
        """

        # modules:
        import numpy

        # store:
        self.nlat = lli.nlat
        self.nlon = lli.nlon
        self.nlat_in = lli_in.nlat
        self.nlon_in = lli_in.nlon
        # output:
        self.jj_in = numpy.zeros((lli.nlat,lli.nlon),int)
        self.ii_in = numpy.zeros((lli.nlat,lli.nlon),int)
        # loop over target grid points:
        for j in range(lli.nlat) :
            for i in range(lli.nlon) :
                # indices:
                j_in,i_in = lli_in.sample_indices( lli.lats[j], lli.lons[i],
                                                   missing_value=-999 )
                # store:
                self.jj_in[j,i] = j_in
                self.ii_in[j,i] = i_in
            #endfor
        #endfor

    #enddef
    
    # *
    
    def getindex( self, j, i ) :
    
        """
        Return sample indices (j_in,i_in) of cell in input grid
        that contains center of cell (j,i) in target grid.
        """
        
        return self.jj_in[j,i],self.ii_in[j,i]
        
    #enddef

    # *

    def apply( self, data_in ) :

        """
        Horizontal remapping.
        Sample data at y and x cells stored in jj and ii.
        """

        # modules:
        import logging
        import numpy

        # input rank:
        ndim = len(data_in.shape)
        # check:
        if ndim < 2 :
            logging.error( 'to apply ll sampling array should have at least rank 2' )
            raise ValueError
        #endif
        # check ...
        if (data_in.shape[ndim-1] != self.nlon_in) or \
           (data_in.shape[ndim-2] != self.nlat_in) :
            logging.error( 'shape of input (%i,%i) does not match with sampling (%i,%i)' % \
                      (data_in.shape[ndim-2],data_in.shape[ndim-1],
                       self.nlat_in,self.nlon_in) )
            raise ValueError
        #endif

        # per dimension:
        if ndim == 2 :
            # original shape:
            ny,nx = data_in.shape
            # output array:
            data = numpy.zeros((self.nlat,self.nlon),float)
            # loop over target cells:
            for j in range(self.nlat) :
                for i in range(self.nlon) :
                    i_in = self.ii_in[j,i]
                    j_in = self.jj_in[j,i]
                    if (i_in < 0) or (j_in < 0) :
                        value = None
                    else :
                        value = data_in[j_in,i_in]
                    #endif
                    data[j,i] = value
                #endfor
            #endfor
        elif ndim == 3 :
            # original shape:
            nz,ny,nx = data_in.shape
            # output array:
            data = numpy.zeros((nz,self.nlat,self.nlon),float)
            # loop over target cells:
            for j in range(self.nlat) :
                for i in range(self.nlon) :
                    i_in = self.ii_in[j,i]
                    j_in = self.jj_in[j,i]
                    if (i_in < 0) or (j_in < 0) :
                        value = None
                    else :
                        value = data_in[:,j_in,i_in]
                    #endif
                    data[:,j,i] = value
                #endfor
            #endfor
        elif ndim == 4 :
            # original shape:
            nt,nz,ny,nx = data_in.shape
            # output array:
            data = numpy.zeros((nt,nz,self.nlat,self.nlon),float)
            # loop over target cells:
            for j in range(self.nlat) :
                for i in range(self.nlon) :
                    i_in = self.ii_in[j,i]
                    j_in = self.jj_in[j,i]
                    if (i_in < 0) or (j_in < 0) :
                        value = None
                    else :
                        value = data_in[:,:,j_in,i_in]
                    #endif
                    data[:,:,j,i] = value
                #endfor
            #endfor
        else :
            logging.error( 'rank not supported yet : %i' % ndim )
            raise ValueError
        #endif

        # ok
        return data

    #enddef # apply


#endclass hori_interp_grid



# **********************************************************************


class hori_aver_grid( object ) :

    """
    Object to store indices to area aver from grid cells to grid cells.
    """

    # *

    def __init__( self, lli, lli_in ) :

        """
        Pre-compute indices to sample data defined on cells lli_in
        into cells defind on lli using aera averages.
        """

        # modules:
        import numpy
        
        # tools:
        import grid_tools
        
        ## testing ...
        #reload(grid_tools)

        # store:
        self.nlat = lli.nlat
        self.nlon = lli.nlon
        self.nlat_in = lli_in.nlat
        self.nlon_in = lli_in.nlon
        
        # mappings for coordinates:
        #  lon_mapping = [ (ii,ww), (ii,ww), .. ]
        self.lon_mapping = grid_tools.cells_mapping( lli_in.blons, lli.blons )
        self.lat_mapping = grid_tools.cells_mapping( lli_in.blats, lli.blats )

    #enddef
    
    # *
    
#    def getindex( self, j, i ) :
#    
#        """
#        Return sample indices (j_in,i_in) of cell in input grid
#        that contains center of cell (j,i) in target grid.
#        """
#        
#        return self.jj_in[j,i],self.ii_in[j,i]
#        
#    #enddef

    # *

    def apply( self, data_in ) :

        """
        Horizontal averaging of lon/lat grids.
        """

        # modules:
        import logging
        import numpy

        # input rank:
        ndim = len(data_in.shape)
        # check:
        if ndim < 2 :
            logging.error( 'to apply ll sampling array should have at least rank 2' )
            raise ValueError
        #endif
        # check ...
        if (data_in.shape[ndim-1] != self.nlon_in) or \
           (data_in.shape[ndim-2] != self.nlat_in) :
            logging.error( 'shape of input (%i,%i) does not match with sampling (%i,%i)' % \
                      (data_in.shape[ndim-2],data_in.shape[ndim-1],
                       self.nlat_in,self.nlon_in) )
            raise ValueError
        #endif

        # per dimension:
        if ndim == 2 :

            # output array:
            data = numpy.zeros((self.nlat,self.nlon),float)
            # loop over target cells:
            for j in range(self.nlat) :
                # source info:
                jj,ww = self.lat_mapping[j]
                # loop over target cells:
                for i in range(self.nlon) :
                    # source info:
                    ii,vv = self.lon_mapping[i]
                    # loop over source cells:
                    for qj in range(len(jj)) :
                        for qi in range(len(ii)) :
                            # add contribution:
                            data[j,i] += data_in[jj[qj],ii[qi]] * ww[qj] * vv[qi]
                        #endfor
                    #endfor
                #endfor # i
            #endfor # j

#        elif ndim == 3 :
#            # original shape:
#            nz,ny,nx = data_in.shape
#            # output array:
#            data = numpy.zeros((nz,self.nlat,self.nlon),float)
#            # loop over target cells:
#            for j in range(self.nlat) :
#                for i in range(self.nlon) :
#                    i_in = self.ii_in[j,i]
#                    j_in = self.jj_in[j,i]
#                    if (i_in < 0) or (j_in < 0) :
#                        value = None
#                    else :
#                        value = data_in[:,j_in,i_in]
#                    #endif
#                    data[:,j,i] = value
#                #endfor
#            #endfor
#
#        elif ndim == 4 :
#            # original shape:
#            nt,nz,ny,nx = data_in.shape
#            # output array:
#            data = numpy.zeros((nt,nz,self.nlat,self.nlon),float)
#            # loop over target cells:
#            for j in range(self.nlat) :
#                for i in range(self.nlon) :
#                    i_in = self.ii_in[j,i]
#                    j_in = self.jj_in[j,i]
#                    if (i_in < 0) or (j_in < 0) :
#                        value = None
#                    else :
#                        value = data_in[:,:,j_in,i_in]
#                    #endif
#                    data[:,:,j,i] = value
#                #endfor
#            #endfor

        else :
            logging.error( 'rank not supported yet : %i' % ndim )
            raise ValueError
        #endif

        # ok
        return data

    #enddef # apply


#endclass hori_aver_grid



# **********************************************************************


class hori_interp_grid( object ) :

    """
    Object to store indices to interpolate from cell centers to cell centers.
    """

    # *

    def __init__( self, lli, lli_in ) :

        """
        Pre-compute indices to interpolate an array defined on
        centers of cells lli_in to centers of cells defined on lli.
        """

        # modules:
        import numpy

        # store:
        self.nlat = lli.nlat
        self.nlon = lli.nlon
        self.nlat_in = lli_in.nlat
        self.nlon_in = lli_in.nlon
        self.lli = lli
        self.lli_in = lli_in

        # target point outside domain ?
        outside = numpy.zeros((self.nlat,self.nlon),bool)
        for iy in range(self.nlat) :
            for ix in range(self.nlon) :
                outside[iy,ix] = not lli_in.in_domain( lli.lats[iy], lli.lons[ix] )
            #endfor
        #endfor
        self.ii_outside = numpy.where( outside )

    #enddef __init__

    # *

    def apply( self, data_in, missing_value=None ) :

        """
        Horizontal interpolation to grid cell centers.
        """

        # modules:
        import logging
        import numpy

        # testing ...
        import go
        import matplotlib.pyplot as plt

        # input rank:
        ndim = len(data_in.shape)
        # check:
        if ndim < 2 :
            logging.error( 'to apply ll sampling array should have at least rank 2' )
            raise ValueError
        #endif
        # check ...
        if (data_in.shape[ndim-1] != self.nlon_in) or \
           (data_in.shape[ndim-2] != self.nlat_in) :
            logging.error( 'shape of input (%i,%i) does not match with sampling (%i,%i)' % \
                      (data_in.shape[ndim-2],data_in.shape[ndim-1],
                       self.nlat_in,self.nlon_in) )
            raise ValueError
        #endif

        # per dimension:
        if ndim == 3 :
            # original shape:
            nt,ny,nx = data_in.shape
            nz = 1
            # output array:
            data = numpy.zeros((nt,nz,self.nlat,self.nlon),float)
            pat_in = numpy.zeros((self.nlat_in,self.nlon),float)
            # loop over 2D fields:
            for it in range(nt) :
                for iz in range(nz) :
                    # linear interpolation, extended over boundaries:
                    for iy in range(self.nlat_in) :
                        pat_in[iy,:] = numpy.interp( self.lli.lons,
                                                 self.lli_in.lons, data_in[it,iy,:] )
                    #endfor
                    for ix in range(self.nlon) :
                        data[it,iz,:,ix] = numpy.interp( self.lli.lats,
                                                 self.lli_in.lats, pat_in[:,ix] )
                    #endfor
                    # mask:
                    if len(self.ii_outside) > 0 :
                        pat = data[it,iz,:,:]
                        pat[self.ii_outside] = missing_value
                        data[it,iz,:,:] = pat
                    #endif
                #endfor
            #endfor
        elif ndim == 4 :
            # original shape:
            nt,nz,ny,nx = data_in.shape
            # output array:
            data = numpy.zeros((nt,nz,self.nlat,self.nlon),float)
            pat_in = numpy.zeros((self.nlat_in,self.nlon),float)
            # loop over 2D fields:
            for it in range(nt) :
                for iz in range(nz) :
                    # linear interpolation, extended over boundaries:
                    for iy in range(self.nlat_in) :
                        pat_in[iy,:] = numpy.interp( self.lli.lons,
                                                 self.lli_in.lons, data_in[it,iz,iy,:] )
                    #endfor
                    for ix in range(self.nlon) :
                        data[it,iz,:,ix] = numpy.interp( self.lli.lats,
                                                 self.lli_in.lats, pat_in[:,ix] )
                    #endfor
                    # mask:
                    if len(self.ii_outside) > 0 :
                        pat = data[it,iz,:,:]
                        pat[self.ii_outside] = missing_value
                        data[it,iz,:,:] = pat
                    #endif
                #endfor
            #endfor
        else :
            logging.error( 'rank not supported yet : %i' % ndim )
            raise ValueError
        #endif

        # ok
        return data

    #enddef # apply


#endclass hori_interp_grid


