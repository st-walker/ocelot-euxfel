import pytest
import numpy as np

from oxfel.predefined import cat_to_i1d
from oxfel.accelerator.lattice import i1, i1d
from oxfel.fel_track import Linac, EuXFELSimConfig

from oxfel.optics import get_match_point_constraints, get_match_point, INJECTOR_MATCHING_QUAD_NAMES

@pytest.fixture
def i1d_real():
    return cat_to_i1d(model_type="real")



def test_linac_match_optics(i1d_real):
    twiss0 = i1d_real.twiss0
    twiss0.beta_x *= 1.05
    twiss0.beta_y *= 1.05

    i1d_real.twiss0 = twiss0
    twiss, _ = i1d_real.machine_twiss()
    df = get_match_point(twiss)

    # Assert beam is not matched (because we changed the inputs a bit)
    assert not np.isclose(df.loc["MATCH.52.I1"]["bmag_x"].item(), 1.0)
    assert not np.isclose(df.loc["MATCH.52.I1"]["bmag_y"].item(), 1.0)

    constraints = get_match_point_constraints()# ["MATCH.52.I1"]
    # constraints = {"MATCH.52.I1": get_match_point_constraints()["MATCH.52.I1"]}
    # Now we match
    felconf = i1d_real.match_twiss(i1d_real.twiss0,
                                   elements=INJECTOR_MATCHING_QUAD_NAMES,
                                   constraints=constraints,
                                   stop="MATCH.52.I1")

    twiss_matched, _ = i1d_real.machine_twiss(felconfig=felconf)
    df = get_match_point(twiss_matched)

    # Should be matched now at the match point.
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_x"].item(), 1.0)
    assert np.isclose(df.loc["MATCH.52.I1"]["bmag_y"].item(), 1.0)



# def test_get_sequence(i1d_real):
#     sequence = i1d_real.get_sequence()
#     assert np.isclose(sequence.length(), i1d_real.length(), atol=1e-3)
#     assert sequence[0].id == i1d_real.sections[0].sequence[0].id
#     assert sequence[-1].id == i1d_real.sections[-1].sequence[-1].id

# def test_get_sequence_with_stop(i1d_real):
#     cathode_to_point = i1d_real.get_sequence(stop=i1.start_sim.id)

#     assert np.isclose(cathode_to_point.length(), 3.2, atol=1e-3)
#     assert cathode_to_point[0].id == "STSEC.23.I1"
#     assert cathode_to_point[-1].id == "start_ocelot"

# def test_to_navigator(i1d_real):
#     navi = i1d_real.to_navigator()

#     sequence = navi.lat.sequence

#     assert np.isclose(i1d_real.length(), navi.lat.totalLen)
#     assert sequence[0].id == i1d_real.sections[0].sequence[0].id
#     assert sequence[-1].id == i1d_real.sections[-1].sequence[-1].id

#     nprocs = sum([len(x.physics_processes_array) for x in i1d_real.sections])

#     assert len(navi.process_table.proc_list) == nprocs

# def test_to_navigator_with_start_stop(i1d_real):
#     navi = i1d_real.to_navigator(start="start_ocelot", stop="MATCH.37.I1")

#     sequence = navi.lat.sequence

#     assert np.isclose(11.2488, navi.lat.totalLen)
#     assert sequence[0].id == "start_ocelot"
#     assert sequence[-1].id == "MATCH.37.I1"

#     # nprocs = sum([len(x.physics_processes_array) for x in i1d_real.sections])

#     assert len(navi.process_table.proc_list) == 4
