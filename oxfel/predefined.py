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


def cat_to_i1d():
    return SectionedFEL([G1, A1, AH1, LH, I1D_Screen], twiss0=i1.tws0)


def cat_to_b2d():
    return SectionedFEL(
        [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D], twiss0=i1.tws0
    )


# def print_fixed_match_points_twiss(twiss_df):
#     for match_point in FIXED_MATCH_POINTS:
#         row = twiss_df[twiss_df[match_point]]
