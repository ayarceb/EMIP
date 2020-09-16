"""
Grid tools
"""


# tools:
import binas


# ***


def ll_area( west, east, south, north, radius=binas.ae ) :

    """
    Given rectangle [west,east]x[south,north] in degrees, 
    compute area.
    Default radius is the radius of the earth (in m) 
    taken from the Binas module, unit of the result is then m2 .
    """

    # modules:
    import math

    # evaluate:
    A = ( math.sin(math.radians(max(north,south))) - math.sin(math.radians(min(north,south))) ) * math.radians(abs(east-west))  # rad^2
    A = A * radius**2

    # ok
    return A

#enddef  # ll_area


# ***


def ll_distance( lon1, lat1, lon2, lat2, radius=binas.ae ) :

    """
    Calculates the distance between two points given their (lat, lon) co-ordinates.
    It uses the Spherical Law Of Cosines (http://en.wikipedia.org/wiki/Spherical_law_of_cosines):
    
    cos(c) = cos(a) * cos(b) + sin(a) * sin(b) * cos(C)                        (1)
    
    In this case:
    a = lat1 in radians, b = lat2 in radians, C = (lon2 - lon1) in radians
    and because the latitude range is  [-pi/2, pi/2] instead of [0, pi]
    and the longitude range is [-pi, pi] instead of [0, 2pi]
    (1) transforms into:
    
    x = cos(c) = sin(a) * sin(b) + cos(a) * cos(b) * cos(C)
    
    Finally the distance is arccos(x)
    
    (source: http://cyberpython.wordpress.com/2010/03/31/python-calculate-the-distance-between-2-points-given-their-coordinates/)
    """

    # modules:
    import math
    
    # quick fix ..
    if (lon1 == lon2) and (lat1 == lat2) :
    
        # same location ...
        distance = 0.0
        
    else :

        delta = lon2 - lon1
        a = math.radians(lat1)
        b = math.radians(lat2)
        C = math.radians(delta)
        x = math.sin(a) * math.sin(b) + math.cos(a) * math.cos(b) * math.cos(C)

        # distance in radians:
        phi = math.acos(x)   # radians

        # in m over the globe:
        distance = phi * radius
        
    #endif
    
    # ok
    return distance

#enddef  # ll_distance

# *

def cell_mapping( bounds, xx ) :

    """
    Input is list of boundary values (0:nx) and target range x(2).
    Return indices in bounds and weights for fractions of bounds that sum up to range.
    """
    
    # number of bound cells:
    nb = len(bounds) - 1

    # find start cell:
    ib1 = None
    for ib in range(nb) :
        if (xx[0] >= bounds[ib]) and (xx[0] <= bounds[ib+1]) :
            ib1 = ib
            break
        #endif
    #endfor
    if ib1 is None :
        print( 'could not find start value %f in bounds:' % xx[0] )
        print( bounds )
        raise Exception
    #endif

    # find end cell:
    ib2 = None
    for ib in range(nb) :
        if (xx[1] >= bounds[ib]) and (xx[1] <= bounds[ib+1]) :
            ib2 = ib
            break
        #endif
    #endfor
    if ib2 is None :
        if abs(bounds[nb]-xx[1]) < 1.0e-4 :
            ib2 = nb-1
        else :
            print( 'could not find end value %f in bounds:' % xx[1] )
            print( bounds )
            raise Exception
        #endif
    #endif
    
    # init results:
    ii = []
    ww = []
    # loop over source cells:
    for ib in range(ib1,ib2+1) :
        # overlapping range:
        d1 = max( bounds[ib], xx[0] )
        d2 = min( bounds[ib+1], xx[1] )
        ## testing ...
        #print( '  x ', ib, d1, d2, (d2-d1), (bounds[ib+1]-bounds[ib]) )
        # weight:
        w = (d2-d1) / (bounds[ib+1]-bounds[ib])
        # store index:
        ii.append( ib )
        ww.append( w )
    #endfor # source cells:
    
    # ok
    return (ii,ww)
    
#enddef cell_mapping

# *

def cells_mapping( bounds_in, bounds ) :

    """
    Return list to map input cells defined by bounds_in to output cells:
      [(ii,ww), (ii,ww), ...]
    """
    
    # init result:
    mapping = []
    # loop over target cells:
    for i in range(len(bounds)-1) :
        # add mapping:
        mapping.append( cell_mapping( bounds_in, bounds[i:i+2] ) )
        ## testing ..
        #ii,ww = cell_mapping( bounds_in, bounds[i:i+2] )
        #print( 'xxx i         = ', i )
        #print( 'xxx bounds_in = ', bounds_in[0:4] )
        #print( 'xxx cell      = ', bounds[i:i+2] )
        #print( 'xxx ii        = ', ii )
        #print( 'xxx ww        = ', ww )
        #raise Exception
    #endfor
    
    ## testing: total source weights:
    #import numpy
    #s = numpy.zeros( (len(bounds_in),) )
    #for tup in mapping :
    #    ii,ww = tup
    #    for k in range(len(ii)) :
    #        s[ii[k]] += ww[k]
    #    #endfor
    ##endfor
    #print( 'xxx s = ', s )
    
    # ok
    return mapping
    
#enddef cells_mapping
    


