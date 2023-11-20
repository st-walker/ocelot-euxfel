import os
from importlib_resources import files
from collections import defaultdict
from pathlib import Path
from collections.abc import MutableSequence
import logging

import yaml
from ocelot.cpbd.csr import CSR
from ocelot.cpbd.physics_proc import SmoothBeam, PhysProc
from ocelot.cpbd.wake3D import Wake, WakeTable
from ocelot.cpbd.sc import SpaceCharge

from .wakes import find_wake_file
from .sequence import MachineSequence


DEFAULT_PROCESSES = files("oxfel.processes") / "processes.yaml"
LOG = logging.getLogger(__name__)


class RaisingProc(PhysProc):
    def __init__(self, exception_message):
        super().__init__()
        self.exception_message = exception_message

    def apply(self, parray, dz):
        raise DontTrackHereError(self.exception_message)


class PlacedPhysics:
    def __init__(self, process, start_name, stop_name, description="", active=True):
        self.start_name = start_name
        self.stop_name = stop_name
        self.process = process
        self.description = description
        self.active = active

    def __repr__(self) -> str:
        proc = self.process
        start = self.start_name
        stop = self.stop_name
        desc = self.description
        return f'<{type(proc).__name__} from "{start}" to "{stop}". Description: "{desc}">'


# class PrecalculatingCSR..

class PlacedElementPhysics(PlacedPhysics):
    def __init__(self, proc, start_name, stop_name, description="", active=True):
        self.proc = proc
        self.start_name = str(start_name)
        self.stop_name = str(stop_name)
        self.active = active

    def modify_element(self, element):
        pass

    def modify_sequence(self, element):
        pass


class LaserHeaterProcess(PlacedElementPhysics):
    def __init__(self, proc, start, stop, undulator_name, description=""):
        pass


class PhysicsList(MutableSequence):
    def __init__(self, processes=()):
        self._processes = list(processes)

    def __getitem__(self, i):
        return self._processes[i]

    def __setitem__(self, i, value):
        self._processes[i] = value

    def __delitem__(self, i):
        del self._processes[i]

    def insert(self, i, value):
        self._processes.insert(i, value)

    def __len__(self):
        return len(self._processes)
    

    def __repr__(self):
        string = "\n".join(repr(x) for x in self)
        return f"<{type(self).__name__}: {string}>"


def load_physics_lists_from_file(fname: os.PathLike, process_builder=None) -> PhysicsList:
    with open(fname, "r") as f:
        config = yaml.safe_load(f)

    if process_builder is None:
        process_builder = ProcessBuilder()

    result = defaultdict(PhysicsList)
    for section_name in config:
        physics_definitions = config[section_name]["physics"]
        for pdict in physics_definitions:
            procname = pdict["process"]

            process = process_builder.make_process(procname, pdict["kwargs"],
                                                   # start=pdict["start"],
                                                   # stop=pdict["stop"]
                                                   )


            placed_proc = PlacedPhysics(process,
                                        pdict["start"],
                                        pdict["stop"],
                                        description=pdict["description"])

            result[section_name].append(placed_proc)

    return result


class ProcessBuilder:
    DISPATCH = {"Wake": "make_wake",
                "SmoothBeam": "make_smooth_beam",
                "CSR": "make_csr",
                "SpaceCharge": "make_space_charge",
                "RaisingProc": "make_raising_proc"}

    def make_process(self, procname, kwargs):
        method_name = self.DISPATCH[procname]
        method = getattr(self, method_name)
        return  method(**kwargs)

    def make_wake(self, filename, *, factor, step=None, w_sampling=None, filter_order=None):
        wake = Wake()
        
        wake.wake_table = WakeTable(find_wake_file(filename))
        wake.factor = factor
        if step is not None:
            wake.step = step

        if w_sampling is not None:
            wake.w_sampling = w_sampling
        if filter_order is not None:
            wake.filter_order = filter_order

        return wake

    def make_space_charge(self, *, step, nmesh_xyz=None, random_mesh=None):
        sc = SpaceCharge()
        sc.step = step
        if nmesh_xyz is not None:
            sc.nmesh_xyz = nmesh_xyz
        if random_mesh is not None:
            sc.random_mesh = random_mesh
        return sc

    def make_beam_transform(self, *, betx, bety, alfx, alfy, tr_slice):
        tws = Twiss()
        tws.beta_x = betx
        tws.beta_y = bety
        tws.alpha_x = alfx
        tws.alpha_y = alfy
        tws.gamma_x = (1 + alfx**2) / betx
        tr = BeamTransform(tws=tws)
        tr.slice = tr_slice
        return tr

    def make_csr(self, *, sigma_min, traj_step, apply_step, n_bin=None, step=None):
        csr = CSR()
        csr.sigma_min = sigma_min
        csr.traj_step = traj_step
        csr.apply_step = apply_step
        if n_bin is not None:
            csr.n_bin = n_bin
        if step is not None:
            csr.step = step
        return csr

    def make_raising_proc(self, exception_message):
        return RaisingProc(exception_message=exception_message)

    def make_laser_modulator(self, dE=300.):
        lh = LaserModulator()
        lh.dE = dE
        lh.sigma_l = 300
        lh.sigma_x = 300e-6
        lh.sigma_y = 300e-6
        lh.z_waist = None
        return lh

    def make_smooth_beam(self, mslice):
        smproc = SmoothBeam()
        smproc.mslice = mslice
        return smproc


def load_default_physics_lists() -> PhysicsList:
    return load_physics_lists_from_file(DEFAULT_PROCESSES)
