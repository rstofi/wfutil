
"""
wfutil
========
Python package of utility functions used to work with:
    time arrays of (1xN)
    frequency arrays (1xM)
    data arrays in the time-frequency space (NxM or MxN)

developed for handling waterfall and flag plots from interferometric data

Author(s): Kristof Rozgonyi (@rstofi)
"""
from . import core, time, msutil, visutil, miscutil

#from .version import version as __version__

__name__ = "wfutil"
__author__ = ["Kristof Rozgonyi (@rstofi)"]

__all__ = ['core', 'time', 'msutil', 'visutil', 'miscutil']
