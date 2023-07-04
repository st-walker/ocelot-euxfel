"""Top-level package for oxfel."""

__author__ = """Stuart Walker"""
__email__ = "stuart.walker@desy.de"
__version__ = "0.1.0"


from .accelerator import lattice, sections, wakes
from . import fel_track
from .predefined import cat_to_i1d, cat_to_b2d
from .longlist import XFELLongList
from .optics import FIXED_MATCH_POINTS, get_match_point_optics
