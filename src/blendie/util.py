# Created by Shoma Yamanouchi (shoma.yamanouchi@ec.gc.ca)

######## Requirements ########
# Python 3+ (3.9+ recommended) 
# numpy
# datetime
# metpy
# herbie
# grib2io
##############################

# IMPORTANT NOTICE: This package will ONLY work in North America #


import logging
from decimal import Decimal


NUMBER_TYPES = (int, float, Decimal)

__version__ = "0.0.2"
__version_info__ = (0, 0, 2)

logger = logging.getLogger('blendie')


def get_version():
    return __version__
