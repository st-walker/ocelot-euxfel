"""Top-level package for oxfel."""

__author__ = """Stuart Walker"""
__email__ = "stuart.walker@desy.de"
__version__ = "0.2.0"


from .accelerator import lattice# , sections
from .fel_track import EuXFELSimConfig, PeakEnergyTwissCalculator, PeakCurrentTwissCalculator
from .predefined import cat_to_i1d #, cat_to_b1d, cat_to_b2d, cat_to_tld, all_models
from .longlist import XFELComponentList
from .optics import FIXED_MATCH_POINTS, get_match_point_optics, START_SIM


from .secmatch import match_injector
from .astra import load_reference_0320_100k_distribution
