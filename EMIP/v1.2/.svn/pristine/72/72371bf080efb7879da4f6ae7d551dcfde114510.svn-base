#! /usr/bin/env python


########################################################################
###
### help
###
########################################################################


"""
GO - General Objects

Logging
    Messages are passed to the logging system.

"""


########################################################################
###
### logger
###
########################################################################

# modules:
import sys
import logging

# setup logging if not done yet:
logging.basicConfig( stream=sys.stdout, level=logging.INFO,
                       format='[%(levelname)-8s] %(message)s' )

# get root logger instance:
logger = logging.getLogger()


########################################################################
###
### sub modules
###
########################################################################


#import go_records as records
#import go_array   as array
import go_plot2   as plot2
#import go_match   as match
#import go_math    as math
import go_vector    as vector
import go_mapping   as mapping
        
# testing ..
import importlib
importlib.reload(vector)
importlib.reload(mapping)


########################################################################
###
### end
###
########################################################################
