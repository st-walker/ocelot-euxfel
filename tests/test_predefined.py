import pytest

import numpy as np

from oxfel.predefined import cat_to_i1d
from oxfel.optics import get_match_point_optics, get_match_point



@pytest.fixture
def i1d_real():
    return cat_to_i1d(model_type="real")

@pytest.fixture
def i1d_design():
    return cat_to_i1d(model_type="design")


def test_i1d_match_optics_design(i1d_real):
    """Build the real model, track, and check the optics at the match point are correct."""
    twiss, mlat = i1d_real.machine_twiss()
    df = get_match_point(twiss)
    # Should be perfectly matched.
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_x"].item(), 1.0)
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_y"].item(), 1.0)
    

def test_i1d_match_optics_real(i1d_design):
    """Build the design model, track, and check the optics at the match point are correct."""
    twiss, mlat = i1d_design.machine_twiss()
    df = get_match_point(twiss)
    # Should be perfectly matched.
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_x"].item(), 1.0)
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_y"].item(), 1.0)
    
