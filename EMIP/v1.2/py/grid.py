
"""
Grid tools.
"""

#
# TO BE RESTRUCTURED
#
# New class hierchy:
#
#   Grid               # empty base class
#
#     IrregGrid        # list of points, e.g. stations
#  
#     RegGrid               # 2D raster of grid points
#       RegRectGrid         # 2D rectangular raster of grid points
#       RegGridCells        # 2D raster of grid cells
#         RegRectGridCells  # 2D rectangular raster of grid cells
#

import grid_cg      as cg
#import grid_icg     as icg
#import grid_ncg     as ncg
#
#import grid_tools   as tools
#import grid_vsamp   as vsamp
#import grid_slopes  as slopes
#from grid_cg import CarthesianGrid as llcells
