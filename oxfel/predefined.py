from ocelot.cpbd.beam import Twiss

from oxfel.accelerator.lattice import i1
from oxfel.accelerator.sections import G1, A1, AH1, LH, I1D_Screen
from oxfel.fel_track import SectionedFEL


def cat_to_i1d():
    return SectionedFEL([G1, A1, AH1, LH, I1D_Screen], twiss0=i1.tws)


