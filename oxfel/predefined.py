import toml
from pathlib import Path
import pickle

from oxfel.accelerator.lattice import i1
from oxfel.accelerator.sections import (
    G1,
    A1,
    AH1,
    LH,
    I1D_Screen,
    DL,
    BC0,
    L1,
    BC1,
    L2,
    BC2,
    B2D,
)
from oxfel.fel_track import SectionedFEL
from oxfel.conversion import make_felconfig_design_to_real

I1D_SECTIONS = [G1, A1, AH1, LH, I1D_Screen]
B2D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D]
TWISS0 = i1.tws0

CONVERSION_CONFIG = Path("/Users/stuartwalker/repos/oxfel/oxfel/accelerator/lattice/real-matched-conf.pcl")

def _make_real_felconfig():
    with CONVERSION_CONFIG.open("rb") as f:
        return pickle.load(f)

    # dconf = toml.load(CONVERSION_CONFIG)
    # return make_felconfig_design_to_real(dconf["extras"]["real"])

def _init_i1d_sections():
    return [x() for x in I1D_SECTIONS]

def _init_b2d_sections():
    return [x() for x in B2D_SECTIONS]

def cat_to_i1d_design():
    return SectionedFEL(_init_i1d_sections(), twiss0=TWISS0)

def cat_to_i1d_real():
    return SectionedFEL(_init_i1d_sections(), twiss0=TWISS0, felconfig=_make_real_felconfig())

def cat_to_b2d_design():
    return SectionedFEL(_init_b2d_sections(), twiss0=TWISS0)

def cat_to_b2d_real():
    return SectionedFEL(_init_b2d_sections(), twiss0=TWISS0, felconfig=_make_real_felconfig())
