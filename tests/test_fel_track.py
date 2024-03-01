import pytest
import numpy as np

from oxfel.predefined import cat_to_i1d
from oxfel.accelerator.lattice import i1, i1d
from oxfel.fel_track import Linac, EuXFELSimConfig



@pytest.fixture
def linac():
    return cat_to_i1d(model_type="real")


def test_linac_length(linac):
    pass


def test_get_sequence(linac):
    sequence = linac.get_sequence()
    assert np.isclose(sequence.length(), linac.length(), atol=1e-3)
    assert sequence[0].id == linac.sections[0].sequence[0].id
    assert sequence[-1].id == linac.sections[-1].sequence[-1].id

def test_get_sequence_with_stop(linac):
    cathode_to_point = linac.get_sequence(stop=i1.start_sim.id)

    assert np.isclose(cathode_to_point.length(), 3.2, atol=1e-3)
    assert cathode_to_point[0].id == "STSEC.23.I1"
    assert cathode_to_point[-1].id == "start_ocelot"

def test_to_navigator(linac):
    navi = linac.to_navigator()

    sequence = navi.lat.sequence

    from IPython import embed; embed()
    assert np.isclose(linac.length(), navi.lat.totalLen)
    assert sequence[0].id == linac.sections[0].sequence[0].id
    assert sequence[-1].id == linac.sections[-1].sequence[-1].id

    nprocs = sum([len(x.physics_processes_array) for x in linac.sections])

    assert len(navi.process_table.proc_list) == nprocs

def test_to_navigator_with_start_stop(linac):
    navi = linac.to_navigator(start="start_ocelot", stop="MATCH.37.I1")

    sequence = navi.lat.sequence

    assert np.isclose(11.2488, navi.lat.totalLen)
    assert sequence[0].id == "start_ocelot"
    assert sequence[-1].id == "MATCH.37.I1"

    # nprocs = sum([len(x.physics_processes_array) for x in linac.sections])

    assert len(navi.process_table.proc_list) == 4
