import pytest
import numpy as np


from oxfel.accelerator.lattice import i1, i1d
from oxfel.accelerator.sections import G1, A1, AH1
from oxfel.fel_track import Linac, EuXFELSimConfig

@pytest.fixture
def fel():
    return Linac([G1, A1, AH1], twiss0=i1.tws)

def test_get_sequence(fel):
    sequence = fel.get_sequence()
    assert np.isclose(sequence.length(), fel.length(), atol=1e-3)
    assert sequence[0].id == fel.sections[0].sequence[0].id
    assert sequence[-1].id == fel.sections[-1].sequence[-1].id

def test_get_sequence_with_stop(fel):
    cathode_to_point = fel.get_sequence(stop=i1.start_sim.id)

    assert np.isclose(cathode_to_point.length(), 3.2, atol=1e-3)
    assert cathode_to_point[0].id == "STSEC.23.I1"
    assert cathode_to_point[-1].id == "START_SIM"

def test_to_navigator(fel):
    navi = fel.to_navigator()

    sequence = navi.lat.sequence

    assert np.isclose(fel.length(), navi.lat.totalLen)
    assert sequence[0].id == fel.sections[0].sequence[0].id
    assert sequence[-1].id == fel.sections[-1].sequence[-1].id

    nprocs = sum([len(x.physics_processes_array) for x in fel.sections])

    assert len(navi.process_table.proc_list) == nprocs

def test_to_navigator_with_start_stop(fel):
    navi = fel.to_navigator(start="START_SIM", stop="MATCH.37.I1")

    sequence = navi.lat.sequence

    assert np.isclose(11.2488, navi.lat.totalLen)
    assert sequence[0].id == "START_SIM"
    assert sequence[-1].id == "MATCH.37.I1"

    # nprocs = sum([len(x.physics_processes_array) for x in fel.sections])

    assert len(navi.process_table.proc_list) == 4
