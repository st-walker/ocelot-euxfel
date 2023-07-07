"""Top-level package for oxfel."""

__author__ = """Stuart Walker"""
__email__ = "stuart.walker@desy.de"
__version__ = "0.1.0"


from .accelerator import lattice, sections, wakes
from . import fel_track
from .predefined import cat_to_i1d_design, cat_to_i1d_real, cat_to_b2d_design, cat_to_b2d_real
from .longlist import XFELLongList
from .optics import FIXED_MATCH_POINTS, get_match_point_optics
