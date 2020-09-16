"""
Classes for vectors, lines, and polygons.
"""

class Vector( object ) :
    
    """
    2D vector object.
    """
    
    def __init__( self, x, y ) :
    
        """
        Inititalize vector using positions.
        
        Attributes defined:
          x, y    : position elements.
        """
        
        # modules:
        import numpy
      
        # store:
        self.x = numpy.double(x)
        self.y = numpy.double(y)
        
    #enddef __init__
    
    # *
    
    def __str__( self ) :
    
        """
        Return str with end point.
        """
        
        # ok:
        return '(%f,%f)' % (self.x,self.y)
        
    #enddef __str__
    
    # *
    
    def Length( self ) :
    
        """
        Return length following l2 norm.
        """
        
        # modules:
        import math
        
        # compute length:
        return math.sqrt( self.x**2 + self.y**2 )
        
    #enddef Length
    
    # *
    
    def Scale( self, factor ) :
    
        """
        Scale vector with provided factor.
        """
        
        # scale:
        self.x = self.x * factor
        self.y = self.y * factor
        
    #enddef Scale
    
    # *
    
    def Normalize( self ) :
    
        """
        Normalize vector.
        """
        
        # current length:
        length = self.Length()
        
        # check ...
        if length == 0.0 :
            print( 'ERROR - could not normalize zero-length vector' )
            raise Exception
        #endif
        
        # scale:
        self.Scale( 1.0/length )
        
    #enddef Normalize
    
    # *
    
    def Dot( self, vec ) :
    
        """
        Dot product between vectors.
        """
        
        # compute dot product:
        return self.x * vec.x + self.y * vec.y
        
    #enddef Dot
    
    # *
    
    def Angle( self, seg ) :
    
        """
        Angle in radians between vectors:
                         x.y
          cos(alpha) = -------
                       |x| |y|
        """
        
        # modules:
        import numpy
        
        # ratio:
        cos_alpha = min( max( 0.0, self.Dot( seg) / self.Length() / seg.Length() ), 1.0 )
        # invert:
        return numpy.arccos( cos_alpha )
        
    #enddef Angle
    
    # *
    
    def NormalVector( self ) :
    
        """
        Return normal vector with length 1.
        The vector is perpendicular to the line segment,
        with sign undetermined.
        """
        
        # check length, should be non-zero:
        if self.Length() == 0.0 :
            print( 'ERROR - no normal vector for length-zero vector' )
            raise Exception
        #endif
        
        # Dot product should be zero:
        #   (x,y).(p,q) = 0
        #   x*p + y*q = 0
        # First assume q==1, normalize later on:
        #   xp + y = 0
        #        p = -y/x,  q = 1
        # or if x is zero:
        #   p = 1   ,  q = 0
        if self.x == 0 :
            p = 1.0
            q = - self.x / self.y
        else :
            p = - self.y / self.x
            q = 1.0
        #endif
        
        # init result:
        nvec = Vector( p, q )
        # normalize:
        nvec.Normalize()
        
        # ok:
        return nvec
        
    #enddef NormalVector
    
    # *
    
    def NormalVectorOutward( self, p ) :
    
        """
        Return normal vector with length 1,
        perpendicular to the line segment,
        and outward from the provided point.
        """
        
        # normal vector perpendicular to line,
        # but undefined for sign:
        nvec = self.NormalVector()
        
        #    nvec |
        #      ^  |     / self
        #       \ |   /
        #        \| / 
        #   ------o----------->  self
        #       / | -
        #     /   |   -
        #         |     v p
        #   
        # for outward w.r.t. to p, the dot product should be negative:
        d = nvec.Dot( p )
        # swap if nvec and p are in same direction (positive dot product):
        if d > 0.0 :
            # swap:
            nvec.Scale( -1.0 )
        #endif

        # ok:
        return nvec
        
    #enddef NormalVectorOutward
    
    # *
    
    def Plot( self, ax, **kwargs ) :
      
        """
        Add marker for end point to plot.
        """
        
        # style:
        style = dict( marker='o' )
        style.update( kwargs )
        
        # add arrow:
        ax.plot( [self.x], [self.y], **style )
        
    #enddef Plot
    
    # *
    
    def PlotArrow( self, ax, **kwargs ) :
      
        """
        Add arrow to plot.
        """
        
        # add arrow:
        ax.arrow( 0.0, 0.0, self.x, self.y, 
                   head_width=0.1, length_includes_head=True, 
                   **kwargs )
        
    #enddef PlotArrow

#endclass Vector

# *

class LineSegment( Vector ) :

    """
    Line element defined by base and direction vectors.
    """
    
    def __init__( self, x0, y0, x1, y1 ) :
    
        """
        Define line segment by start and end point.
        """
        
        # init as direction vector:
        Vector.__init__( self, x1-x0, y1-y0 )
        
        # store base:
        self.x0 = x0
        self.y0 = y0
        # end point:
        self.x1 = x1
        self.y1 = y1
        
    #enddef __init__
    
    # *
    
    def __str__( self ) :
    
        """
        Return str with start and end point.
        """
        
        # ok:
        return '(%f,%f) - (%f,%f)' % (self.x0,self.y0,self.x1,self.y1)
        
    #enddef __str__
    
    # *
    
    def NormalLineSegment( self ) :
    
        """
        Return normal line with base at mid of current line.
        """
        
        # mid:
        mx = ( self.x0 + self.x1 )/2.0
        my = ( self.y0 + self.y1 )/2.0
        
        # normal vector:
        nvec = self.NormalVector()
        
        # define segment:
        return LineSegment( mx, my, mx+nvec.x, my+nvec.y )
        
    #enddef NormalLineSegment
    
    # *
    
    def NormalLineSegmentOutward( self, p ) :
    
        """
        Return normal line with base at mid of current line,
        and direction outward to provided point.
        """
        
        # mid:
        mx = ( self.x0 + self.x1 )/2.0
        my = ( self.y0 + self.y1 )/2.0
        
        # p relative to mid:
        pm = Vector( p.x-mx, p.y-my )
        
        # normal vector of direction (from origin), 
        # outward to p-base:
        nvec = self.NormalVectorOutward( pm )
        
        # define segment:
        return LineSegment( mx, my, mx+nvec.x, my+nvec.y )
        
    #enddef NormalLineSegmentOutward
    
    # *
        
    def IntersectionFactors( self, seg ) :
    
        """
        Return factors (fself,fseg) for self and provided segment
        that define the intersection point:
          p = self.base + fself * self.direction
          p = seg.base  + fseg  * seg.direction
        Both are None if self and segment are parallel.
        """
        
        # modules:
        import numpy
        
        # angle:
        alfa = self.Angle(seg)
        # lines never cross if parallel:
        if alfa == 0.0 :
            # no intersection:
            fself,fseg = None,None
        else :
            # solve factors for direction vectors:
            #   b1 + f1*d1 = b2 + f2*d2
            # or as 2x2 problem:
            #            [a1]
            #   [d1,-d2] [a2] = b2 - b1
            f = numpy.linalg.solve( numpy.array([[self.x,-seg.x],[self.y,-seg.y]]),
                                    numpy.array([seg.x0-self.x0,seg.y0-self.y0]) )
            # split:
            fself,fseg = f[0],f[1]
        #endif # parallel?
        
        # ok
        return fself,fseg
        
    #enddef IntersectionFactors
    
    # *
        
    def Intersection( self, seg ) :
    
        """
        Return Vector object with intersection with provided line segment,
        or None of the segments do not cross.
        """
        
        # factors for direction vectors for intersection point:
        fself,fseg = self.IntersectionFactors( seg )
        
        # not defined if parallel:
        if fself is None :
            # no intersection:
            p = None
        else :
            # both on line?
            if (fself >= 0.0) and (fself <= 1.0) and (fseg >= 0.0) and (fseg <= 1.0) :
                # evaluate for first line:
                p = Vector( self.x0+fself*self.x, self.y0+fself*self.y )
            else :
                # intersection outside at least one of the segments:
                p = None
            #endif
        #endif
        
        # ok
        return p
        
    #enddef Intersection
    
    # *
        
    def CutOff( self, line, nvec ) :
    
        """
        Cut the part of the current segment that is outward from the line,
        where outward is defined by nvec (probably a normal vector of the line).
        
        Return values:
          part    : remaining segment, or copy of line if no intersection
          cutoff  : bool
        """
        
        # dot product between direction and normal vector:
        d = self.Dot( nvec )
        # no cut-off if perpendicular (lines are parallel):
        if d == 0.0 :
    
            # no cut off:
            seg = LineSegment( self.x0, self.y0, self.x1, self.y1 )
            # set flag:
            cutoff = False

        else :
        
            # factors for direction vectors for intersection point:
            f,fline = self.IntersectionFactors( line )
            
            # no intersection?
            if (f is None) or (f < 0.0) or (f > 1.0) or (fline is None) or (fline < 0.0) or (fline > 1.0) :

                # no intersection:
                seg = LineSegment( self.x0, self.y0, self.x1, self.y1 )
                # set flag:
                cutoff = False
                
            else :

                #
                #                  |
                #           o----->|- - - - >  self
                #                  |--> nvec                0<f<1, d>0
                #
                #                  |
                #           o - - -|-------->  self
                #               <--| nvec                   0<f<1, d<0
                #
                # same direction?
                if d > 0.0 :
                    # first part:
                    seg = LineSegment( self.x0, self.y0, self.x0+f*self.x, self.y0+f*self.y )
                #
                else :
                    # second part:
                    seg = LineSegment( self.x0+f*self.x, self.y0+f*self.y, 
                                       self.x0+self.x, self.y0+self.y )
                #endif  # direction

                # set flag:
                cutoff = True

            #endif  # intersection?
            
        #endif # parallel?

        # ok
        return seg,cutoff

    #enddef CutOff
    
    # *
    
    def LonLat_LeftArea( self, xb, ae=None ) :
    
        """
        Integral from y-axis to line segment if (x,y) are interpreted 
        as (longitude,latitude) in degrees.
        Result has units (radians*ae)**2, if undefined ae is set to earth radius in m.
        
            y^
             |
           y1|  +-------------o
             |  |/////////////  line: x = a y + b
             |  |////////////
             |  |///////////
           y0|  +---------o
             |         
           --+--+---------+---+----> x
                xb        x0  x1
           
        Line segment as function of y:
          x = (dx/dy) y + x0
          x =    a    y + b
          
        Area in m2, weight with cos(y) to account for smaller cells towards poles:
               y1                     
          A  =  | (a y + b) cos(y) dy 
               y=y0
               
        Use that:

          [y sin(y) + cos(y)]' = sin(y) + y cos(y) - sin(y) = y cos(y)
        
        Thus:
               y1
           A =  | a y cos(y) + b cos(y) dy
               y=y0
                                                   y1
             = [ a {y sin(y) + cos(y)} + b sin(y) ]
                                                   y=y0
        
        """

        # modules:
        import math 
               
        # tools:
        import binas
        
        # radius:
        if ae is None : ae = binas.ae
        
        # line almost parallel to x-axis ?
        if abs(self.y) <= 0.001 * abs(self.x) :

            # no area:
            area = 0.0

        else :
        
            # line definition with base in middle, convert to radiances if necessary:
            a = self.x / self.y  # ratio is the same in rads
            xm   = math.radians( 0.5*(self.x0+self.x1) - xb )
            phim = math.radians( 0.5*(self.y0+self.y1) )
            ## testing ...
            #print ''
            #print 'inty: line ', xm, '+', a, '* ( phi - ', phim,')'

            # integral over positive range, symetric around middle;
            # this gives smaller differences for line definitions
            # that differ in direction:
            if self.y0 < self.y1 :
        
                # convert lat range to radians:
                phi0 = math.radians( self.y0 )
                phi1 = math.radians( self.y1 )
                
            else :
        
                # convert lat range to radians:
                phi0 = math.radians( self.y1 )
                phi1 = math.radians( self.y0 )
            
            #endif

            ## testing ...
            #print 'inty: range ', phi0, phi1

            # evaluate integrants:
            area0 = a * ( (phi0-phim) * math.sin(phi0) + math.cos(phi0) ) + xm * math.sin(phi0)
            area1 = a * ( (phi1-phim) * math.sin(phi1) + math.cos(phi1) ) + xm * math.sin(phi1)
            # combine:
            area = area1 - area0
            ## testing ...
            #print 'inty: area0 ', a * ( phi0 * math.sin(phi0) + math.cos(phi0) ), ' + ', b * math.sin(phi0), ' = ', area0
            #print 'inty: area1 ', a * ( phi1 * math.sin(phi1) + math.cos(phi1) ), ' + ', b * math.sin(phi1), ' = ', area1
            #print 'inty: area ', area1, ' - ', area0, ' = ', area
            
            # swap sign if necessary:
            #if phi1 < phi0 : area = -1.0 * area
            
            # convert to m2:
            area = area * ae**2
            
        #endif
        
        # ok
        return area
        
    #endef LonLat_LeftArea
        
    
    # *
    
    def Plot( self, ax, **kwargs ) :
      
        """
        Add line segment to plot.
        """
        
        # add line:
        ax.plot( (self.x0,self.x1), (self.y0,self.y1), **kwargs )
        
    #enddef Plot
    
    # *
    
    def PlotArrow( self, ax, **kwargs ) :
      
        """
        Add line segment to plot as arrow.
        """
        
        # add arrow:
        ax.arrow( self.x0, self.y0, self.x, self.y, 
                   head_width=0.1, length_includes_head=True, 
                   **kwargs )
        
    #enddef PlotArrow

#endclass LineSegment


# ***


class Polygon( object ) :

    """
    Polygon defined by corner points or edge segments.
    """
    
    def __init__( self, corners=None, edges=None ) :
    
        """
        Define polygon by list of corners or line segments.
        
        Optional arguments:
          corners   : list of Vector objects
        """
        
        # how defined?
        if corners is not None :
        
            # check ..
            if edges is not None :
                print( 'ERROR - polygon could not be defined by both corners and edges.' )
                raise Exception
            #endif
        
            # number of corners:
            self.n = len(corners)

            # check ..
            if self.n < 3 :
                print( 'ERROR - at least 3 corners needed for polygon, received %i:' % self.n )
                for i in range(self.n) :
                    print( 'ERROR -   %s' % corners[i] )
                #endfor
                raise Exception
            #endif
            
            # store:
            self.corners = corners

            # define edge segments:
            self.edges = []
            for i in range(self.n-1) :
                # add segment:
                self.edges.append( LineSegment( corners[i].x, corners[i].y, 
                                                corners[i+1].x, corners[i+1].y ) )
            #endfor
            self.edges.append( LineSegment( corners[self.n-1].x, corners[self.n-1].y, 
                                            corners[0].x, corners[0].y ) )        
        # edges:
        elif edges is not None :
        
            # check ..
            if corners is not None :
                print( 'ERROR - polygon could not be defined by both edges and corners.' )
                raise Exception
            #endif
        
            # number of corners:
            self.n = len(edges)

            # check ..
            if self.n < 3 :
                print( 'ERROR - at least 3 edges needed for polygon, received %i:' % self.n )
                for i in range(self.n) :
                    print( 'ERROR -   %s' % edges[i] )
                #endfor
                raise Exception
            #endif
            
            # init flags:
            used = []
            for i in range(self.n) : used.append( False )

            # init edge segments with first:
            self.edges = [edges[0]]
            used[0] = True
            # loop over remaining edges
            for i in range(1,self.n) :
                # init flag:
                qj = -999
                # loop over edges that are not used yet:
                for j in range(self.n) :
                    # skip if already used:
                    if used[j] : continue
                    # start connected with latest end point?
                    dist = (self.edges[-1].x1 - edges[j].x0)**2 + (self.edges[-1].y1 - edges[j].y0)**2
                    # update?
                    if (qj < 0) or (dist < qdist) :
                        # store current:
                        qj = j
                        qstart = True
                        # update minimum:
                        qdist = dist
                    #endif
                    # end connected to latest end point?
                    dist = (self.edges[-1].x1 - edges[j].x1)**2 + (self.edges[-1].y1 - edges[j].y1)**2
                    # update?
                    if dist < qdist :
                        # stoe:
                        qj = j
                        qstart = False
                        # update minimum:
                        qdist = dist
                    #endif
                #endfor
                # check ..
                if qj < 0 :
                    print( 'ERROR - could not find connection between:' )
                    print( 'ERROR -   %s' % self.edges[-1] )
                    print( 'ERROR - and edges:' )
                    for j  in range(self.n) :
                        print( 'ERROR -   %s  (used %s)' % (edges[j],used[j]) )
                    #endfor
                    raise Exception
                #endif
                # define:
                if qstart :
                    # store copy:
                    self.edges.append( LineSegment( edges[qj].x0, edges[qj].y0, edges[qj].x1, edges[qj].y1 ) )
                else :
                    # store swapped copy:
                    self.edges.append( LineSegment( edges[qj].x1, edges[qj].y1, edges[qj].x0, edges[qj].y0 ) )
                #endif
                # mark as used:
                used[qj] = True
            #endfor # target edges
            
            # corners:
            self.corners = []
            for i in range(self.n) :
                self.corners.append( Vector( self.edges[i].x0, self.edges[i].y0 ) )
            #endfor
        
        else :
            print( 'ERROR - need corners to define polygon' )
            raise Exception
        #endif
        
        # center point:
        xm = 0.0
        ym = 0.0
        for i in range(self.n) :
            xm = xm + self.corners[i].x
            ym = ym + self.corners[i].y
        #endfor
        xm = xm / self.n
        ym = ym / self.n
        # define:
        self.center = Vector( xm, ym )
        
        # outward normal vectors:
        self.normals = []
        for i in range(self.n) :
            # add normal vector directed outward from center:
            self.normals.append( self.edges[i].NormalLineSegmentOutward( self.center ) )
        #endfor
        
    #enddef __init__
    
    # *
    
    def Inside( self, p ) :
      
        """
        True if p is inside polygon.
        """
        
        # line from center to p:
        line = LineSegment( self.center.x, self.center.y, p.x, p.y )
        
        # find intersection with edge; 
        # init index, remains negative if no edge is crossed and p is inside:
        k = -999
        # loop:
        for i in range(self.n) :
            # intersection point, or None:
            c = self.edges[i].Intersection( line )
            # found ?
            if c is not None :
                # store:
                k = i
                # leave:
                break
            #endif
        #endfor
        
        # inside if no intersection was found:
        return k < 0
        
    #enddef Inside
    
    # *
    
    def InnerLineSegment( self, seg ) :
    
        """
        Return part of line segment that is inside current polygon,
        or None if completely outside.
        """
        
        # init result as copy:
        part = LineSegment( seg.x0, seg.y0, seg.x1, seg.y1 )
        # set flag:
        anycutoff = False
        # loop over edges:
        for i in range(self.n) :
            # cut of part that is outside, None if no intersection:
            part,cutoff = part.CutOff( self.edges[i], self.normals[i] )
            # update flag:
            anycutoff = anycutoff or cutoff
        #endfor
        
        # not any cutoff ? then completely inside or outside:
        if not anycutoff :
            # start point:
            p = Vector( part.x0, part.y0 )
            # is p inside?
            inside = self.Inside( p )
            # reset to empty if outisde:
            if not inside : part = None
        #endif
        
        # ok
        return part
        
    #enddef InnerLineSegment

    # *
    
    def Intersection( self, pg ) :
    
        """
        Return polygon of intersection of current and external polygon.
        """
        
        ## testing ...
        #print( '' )

        # init collection of edges:
        edges = []
        # loop over own edges:
        for i in range(self.n) :
            # part inside other polygon:
            part = pg.InnerLineSegment( self.edges[i] )
            # store if not empty:
            if (part is not None) and (part.Length() > 0.0) :
                #print( 'xxx added part of own  edge ', i+1, str(part), ' with length %e' % part.Length() )
                edges.append( part )
            #endif
        #endfor
        # loop over external edges:
        for i in range(pg.n) :
            # part inside current polygon:
            part = self.InnerLineSegment( pg.edges[i] )
            # store if not empty:
            if (part is not None) and (part.Length() > 0.0) :
                #print( 'xxx added part of poly edge ', i+1, str(part), ' with length %e' % part.Length() )
                edges.append( part )
            #endif
        #endfor
        
        # empty?
        if len(edges) == 0 :
            # empty:
            poly = None
        else :
            # define new polygon by edges:
            #for i in range(len(edges)) : print( 'xxx edge ', i, str(edges[i]), edges[i].Length() )
            poly = Polygon( edges=edges )
        #endif
        
        # ok:
        return poly
    
    #enddef Intersection
    
    # *
    
    def LonLat_Area( self, ae=None ) :
    
        """
        Polygon area if (x,y) is interpeted as (longitude,latitude).
        Result in units (radians*ae)**2 .
        """
        
        # left most value:
        xb = self.corners[0].x
        for i in range(1,self.n) :
            xb = min( xb, self.corners[i].x )
        #endfor
        
        # init sum:
        area = 0.0
        # loop over edge segments:
        for i in range(self.n) :
            # contribution:
            inty = self.edges[i].LonLat_LeftArea( xb, ae=ae )
            # sign is positive if normal vector points to right of y-axis,
            # and negative if normal vector points to left ;
            # sufficent to check x-value:
            if self.normals[i].x < 0.0 : inty = -1.0 * inty
            ## testing ...
            #print( 'area: edge %i  %s   : %f' % (i+1,self.edges[i],inty) )
            # add contribution:
            area = area + inty
        #endfor
        ## testing ...
        #print( 'area: total : %f' % area )
        
        # ok
        return area
        
    #enddef LonLat_Area
    
    # *
    
    def PlotFill( self, ax, **kwargs ) :
    
        """
        Add polygon to ax.
        """
        
        # collect corners:
        xx = []
        yy = []
        for i in range(self.n) :
            xx.append( self.corners[i].x )
            yy.append( self.corners[i].y )
        #endfor
        
        # add polygon:
        ax.fill( xx, yy, **kwargs )
        
    #enddef PlotFill
    
    # *
    
    def PlotCenter( self, ax, **kwargs ) :
    
        """
        Add center point to plot.
        """
        
        # add:
        self.center.Plot( ax, **kwargs )
        
    #enddef PlotCenter
    
    # *
    
    def PlotEdges( self, ax, **kwargs ) :
    
        """
        Add edge segments to plot.
        """
        
        # loop over segments:
        for i in range(self.n) :
            # add line:
            self.edges[i].Plot( ax, **kwargs )
        #endfor
        
    #enddef PlotEdges
    
    # *
    
    def PlotNormals( self, ax, **kwargs ) :
    
        """
        Add edge normal vectors to plot as arrows.
        """
        
        # loop over segments:
        for i in range(self.n) :
            # add line:
            self.normals[i].PlotArrow( ax, **kwargs )
        #endfor
        
    #enddef PlotNormals
    
    # *
    
    def BoundingBox( self ) :
    
        """
        Return (west,east,south,north) bounding box enclosing polygon.
        """
        
        # loop over corners:
        for i in range(self.n) :
            # first?
            if i == 0 :
                west  = self.corners[i].x
                east  = self.corners[i].x
                south = self.corners[i].y
                north = self.corners[i].y
            else :
                west  = min( west , self.corners[i].x )
                east  = max( east , self.corners[i].x )
                south = min( south, self.corners[i].y )
                north = max( north, self.corners[i].y )
            #endif
        #endfor
        
        # ok
        return west,east,south,north
    
    #enddef BoundingBox
        

#endclass Polygon

## =====================================================================
## test
## =====================================================================
#
## main program?
#if __name__ == '__main__' :
#  
#    # modules:
#    import matplotlib.pyplot as plt
#    import math
#    
#    # tools:
#    import binas
#    import grid
#    
#    # vector, lines, angles
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#
#        # new point:
#        p0 = Vector( 2, 1 )
#        # add vector:
#        p0.PlotArrow( ax, color='b', linestyle='-' )
#        
#        # define line:
#        line1  = LineSegment( 1,1, 4,5 )
#        line1p = LineSegment( 1,2, 4,6 )  # parallel to line1
#        line2  = LineSegment( 3,4.5, 5,1 )
#        line2s = LineSegment( 3,2, 3.7,1 )  # shorter than line2
#        # add lines:
#        line1.PlotArrow( ax, color='r', linestyle='-' )
#        line1p.PlotArrow( ax, color='r', linestyle='-' )
#        line2.PlotArrow( ax, color='green', linestyle='-' )
#
#        # info ...
#        print( 'angle between line1 and line1p : %f deg' % numpy.degrees(line1.Angle(line1p)) )
#        print( 'angle between line1 and line2  : %f deg' % numpy.degrees(line1.Angle(line2)) )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # plot vector, lines, angles
#    
#    # intersection
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#
#        # new point:
#        line1  = LineSegment( 1,1, 4,5 )
#        line2  = LineSegment( 3,4.5, 5,1 )
#        line2s = LineSegment( 3,2, 3.7,1 )  # shorter than line2
#
#        # add lines:
#        line1.PlotArrow( ax, color='r', linestyle='-' )
#
#        line2.PlotArrow( ax, color='green', linestyle='-' )
#        c12 = line1.Intersection( line2 )
#        if c12 is not None : c12.Plot( ax, color='green' )
#
#        line2s.PlotArrow( ax, color='green', linestyle='-' )
#        c12s = line1.Intersection( line2s )
#        if c12s is not None : c12s.Plot( ax, color='g' )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # plot vector, lines, angles
#    
#    # normal vectors
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#
#        # new point:
#        p0 = Vector( 2, 1 )
#        line1  = LineSegment( 1,1, 4,5 )
#
#        # add vector:
#        p0.PlotArrow( ax, color='b', linestyle='-' )
#        # normal vector:
#        nvec1 = p0.NormalVector()
#        nvec1.PlotArrow( ax, color='b', linestyle='--' )
#        # normal vector outward to point:
#        p1 = Vector( 0.5, 0.7 )
#        p1.Plot( ax, color='purple' )
#        nvec2 = p0.NormalVectorOutward( p1 )
#        nvec2.PlotArrow( ax, color='purple', linestyle='--' )
#
#        # add lines:
#        line1.PlotArrow( ax, color='r', linestyle='-' )
#        # normal segment:
#        nseg1 = line1.NormalLineSegment()
#        nseg1.PlotArrow( ax, color='r', linestyle='--' )
#        # normal vector outward to point:
#        p2 = Vector( 2.5, 4 )
#        p2.Plot( ax, color='orange' )
#        nseg2 = line1.NormalLineSegmentOutward( p2 )
#        nseg2.PlotArrow( ax, color='orange', linestyle='--' )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # plot normal vectors
#    
#    # cut segment
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#
#        # main line:
#        line1  = LineSegment( 1,1, 4,2.5 )
#        # cut-off line:
#        #line2  = LineSegment( -0.5,2, 1,0.2 )
#        line2  = LineSegment( 2,3.5, 4,1 )
#        #line2  = LineSegment( 4, 4, 4.8, 2.2 )
#        # "mid" point
#        #p2 = Vector(0,0)
#        p2 = Vector(5,4)
#        # outward vector:
#        nvec2   = line2.NormalLineSegmentOutward(p2)
#        res,cutoff = line1.CutOff( line2, nvec2 )
#        if not cutoff : print( 'WARNING - no cutoff' )
#
#        # add lines:
#        line1.PlotArrow( ax, color='r' )#, linestyle='-' )
#        p2.Plot( ax, color='g' )
#        line2.PlotArrow( ax, color='g' )#, linestyle='-' )
#        nvec2.PlotArrow( ax, color='g' )#, linestyle='--' )
#        if res is not None : res.PlotArrow( ax, color='purple' )#, linestyle='-', linewidth=3 )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # cut segment
#    
#    # *
#
#    # polygon
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#        
#        # define corners
#        corners = []
#        corners.append( Vector( 1, 1 ) )
#        #corners.append( Vector( 2.5, 2 ) )
#        corners.append( Vector( 3, 1.2 ) )
#        corners.append( Vector( 4, 4 ) )
#        corners.append( Vector( 2, 3.5 ) )
#        # define polygon:
#        poly1 = Polygon( corners=corners )
#        # show:
#        poly1.PlotFill( ax, color='red', alpha=0.5 )
#        poly1.PlotEdges( ax, color='red' )
#        #poly1.PlotCenter( ax, color='red' )
#        #poly1.PlotNormals( ax, color='red', linestyle='--' )
#        
#        # line crossing polygon:
#        #line = LineSegment( 0.5, 1, 5, 2 )
#        line = LineSegment( 1.5, 2, 3, 3 )
#        #line = LineSegment( 1.5, 0, 3, 1 )
#        # show:
#        line.Plot( ax, color='green', linestyle='-' )
#
#        # inner part:
#        line1 = poly1.InnerLineSegment( line )
#        # show:
#        if line1 is not None :
#            line1.Plot( ax, color='green', linestyle='-', linewidth=3 )
#        #endif
#        
#        # define corners
#        corners = []
#        x0,y0 = 2,0
#        corners.append( Vector( x0+0, y0+0.2 ) )
#        corners.append( Vector( x0+2, y0+0.6 ) )
#        corners.append( Vector( x0+1.5, y0+2 ) )
#        corners.append( Vector( x0+0.2, y0+1.8 ) )
#        # define polygon:
#        poly2 = Polygon( corners=corners )        
#        # show:
#        poly2.PlotFill( ax, color='blue', alpha=0.5 )
#        poly2.PlotEdges( ax, color='blue' )
#        
#        ## testing ...
#        #poly1.edges[0].Plot( ax, color='green', linestyle='-' )
#        #qline = poly2.InnerLineSegment( poly1.edges[0] )
#        #if qline is not None : qline.Plot( ax, color='green', linestyle='-', linewidth=3 )
#
#        # intersection:
#        pg = None
#        pg = poly1.Intersection( poly2 )
#        # add ?
#        if pg is not None :
#            pg.PlotFill( ax, color='green', alpha=0.5 )
#            pg.PlotEdges( ax, color='green' )
#        #endif
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # plot normal vectors
#    
#    # *
#    
#    # polygon area
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#        
#        # define corners
#        x0,y0 = 1,1
#        #x0,y0 = -1,-1.2
#        corners = []
#        corners.append( Vector( x0+0, y0+0 ) )
#        corners.append( Vector( x0+2, y0+0.001 ) )
#        corners.append( Vector( x0+3, y0+2.001 ) )
#        corners.append( Vector( x0+1, y0+2 ) )
#        # define polygon:
#        poly1 = Polygon( corners=corners )
#        
#        # add:
#        poly1.PlotFill( ax, color='red', alpha=0.5 )
#        poly1.PlotEdges( ax, color='red' )
#        
#        # area ...
#        ae = numpy.degrees(1.0)
#        #ae = 1.0
#        A1 = poly1.LonLat_Area( ae=ae )
#        # info ...
#        print( 'area = %f' % A1 )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # polygon area
#    
#    # *
#    
#    # polygon overlap
#    if True :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#        
#        # define corners
#        corners = []
#        corners.append( Vector( 1, 1 ) )
#        corners.append( Vector( 3, 1.2 ) )
#        corners.append( Vector( 4, 4 ) )
#        corners.append( Vector( 2, 3.5 ) )
#        # define polygon:
#        poly1 = Polygon( corners=corners )
#        
#        # define corners
#        corners = []
#        x0,y0 = 2,0
#        corners.append( Vector( x0+0, y0+0.2 ) )
#        corners.append( Vector( x0+2, y0+0.6 ) )
#        corners.append( Vector( x0+1.5, y0+2 ) )
#        corners.append( Vector( x0+0.2, y0+1.8 ) )
#        # define polygon:
#        poly2 = Polygon( corners=corners )
#        
#        # add:
#        poly1.PlotFill( ax, color='red', alpha=0.5 )
#        poly1.PlotEdges( ax, color='red' )
#        
#        # add:
#        poly2.PlotFill( ax, color='blue', alpha=0.5 )
#        poly2.PlotEdges( ax, color='blue' )
#        
#        # intersection:
#        pg = None
#        pg = poly1.Intersection( poly2 )
#        # add ?
#        if pg is not None :
#            pg.PlotFill( ax, color='green', alpha=0.5 )
#            pg.PlotEdges( ax, color='green' )
#        #endif
#        
#        # area ...
#        ae = numpy.degrees(1.0)
#        #ae = 1.0
#        A1 = poly1.LonLat_Area( ae=ae )
#        A2 = poly2.LonLat_Area( ae=ae )
#        A3 = pg.LonLat_Area( ae=ae )
#        # info ...
#        print( 'area first   = %f' % A1 )
#        print( 'area second  = %f' % A2 )
#        print( 'area overlap = %f' % A3 )
#
#        # canvas:
#        xlim  = [-1,6]
#        ylim  = [-1,6]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # polygon overlap
#    
#    # *
#    
#    # grid
#    if False :
#  
#        # new figure:
#        fig = plt.figure()
#        ax = fig.add_axes([0.1,0.1,0.8,0.8])
#        ax.set_aspect('equal')
#        
#        # corner:
#        west,south = 5.0,52.0
#        #west,south = -1.1,-2.0
#        
#        # shape:
#        nlat,nlon = 8,4
#        # grid:
#        lli = grid.cg.CarthesianGrid( west=west, dlon=0.5, nlon=nlon, 
#                                      south=south, dlat=0.25, nlat=nlat )
#                                      
#        # corner grids:
#        cxx,cyy = numpy.meshgrid(lli.blons,lli.blats)
#        # disturb?
#        cxx = cxx + numpy.random.random((nlat+1,nlon+1))*0.1
#        cyy = cyy + numpy.random.random((nlat+1,nlon+1))*0.1
#        
#        # same as polygons:
#        pols = []
#        for j in range(lli.nlat) :
#            row = []
#            for i in range(lli.nlon) :
#                # define corners:
#                corners = []
#                corners.append( Vector( cxx[j  ,i  ], cyy[j  ,i  ] ) )
#                corners.append( Vector( cxx[j  ,i+1], cyy[j  ,i+1] ) )
#                corners.append( Vector( cxx[j+1,i+1], cyy[j+1,i+1] ) )
#                corners.append( Vector( cxx[j+1,i  ], cyy[j+1,i  ] ) )
#                # add polygon:
#                row.append( Polygon( corners=corners ) )
#            #endfor # i
#            # append:
#            pols.append( row )
#        #endfor # j
#        
#        # define corners
#        corners = []
#        x0,y0 = west+0.4,south+0.3
#        corners.append( Vector( x0+0, y0+0.2 ) )
#        corners.append( Vector( x0+0.85, y0+0.3 ) )
#        corners.append( Vector( x0+0.65, y0+0.8 ) )
#        corners.append( Vector( x0-0.1, y0+0.6 ) )
#        # define polygon:
#        poly2 = Polygon( corners=corners )
#        
#        # add polygons:
#        for j in range(nlat) :
#            for i in range(nlon) :
#                # add to plot:
#                pols[j][i].PlotFill( ax, color='red', alpha=0.5 )
#                pols[j][i].PlotEdges( ax, color='red' )
#            #endfor
#        #endfor
#        
#        # add:
#        poly2.PlotFill( ax, color='blue', alpha=0.5 )
#        poly2.PlotEdges( ax, color='blue' )
#        
#        # radius to be used:
#        ae = binas.ae
#        ae = numpy.degrees(1.0)
#        # area:
#        area2 = None
#        area2 = poly2.LonLat_Area(ae=ae)
#        # init area sum:
#        area = 0.0
#        # loop over polygons:
#        for j in range(nlat) :
#            for i in range(nlon) :
#                ## testing ...
#                #if (j+1,i+1) != (3,3) : continue
#                # intersection:
#                pg = pols[j][i].Intersection( poly2 )
#                # defined ?
#                if pg is not None :
#                    # add to plot:
#                    pg.PlotFill( ax, color='green', alpha=0.5 )
#                    pg.PlotEdges( ax, color='green' )
#                    #pg.PlotNormals( ax, color='green', linestyle='--' )
#                    # overlap with poly:
#                    pg_area = pg.LonLat_Area( ae=ae )
#                    area += pg_area
#                    #print 'aaa', j+1, i+1, pg_area, area
#                #endif
#            #endfor
#        #endfor
#        # show:
#        if area2 is not None :
#            rdiff = (area-area2)/area2*100.0
#            print( 'area poly2     : %f' % area2 )
#            print( 'area overlaps  : %f' % area )
#            print( 'area rel.diff. : %8.4f %%' % rdiff )
#        #endif
#
#        # canvas:
#        xlim  = [west-0.5,west+2.5]
#        ylim  = [south-0.5,south+2.5]
#        plt.plot( [0,0], ylim, 'k-' )
#        plt.plot( xlim, [0,0], 'k-' )
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)
#
#    #endif # polygon area
#        
##endif  # main
#
## =====================================================================
## end
## =====================================================================
#
#
