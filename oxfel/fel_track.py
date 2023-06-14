"""
Section class for s2e tracking.
S.Tomin. XFEL/DESY. 2017
"""
from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Optional, Union, Type, Iterable
from pathlib import Path
import logging
from copy import deepcopy
import textwrap
import re
import pandas as pd

import numpy as np

from ocelot.cpbd.csr import CSR
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.physics_proc import PhysProc, SaveBeam
from ocelot.cpbd.io import save_particle_array, load_particle_array
import ocelot.cpbd.track as track
from ocelot.cpbd.magnetic_lattice import (
    MagneticLattice,
    insert_markers_by_predicate,
    flatten,
)
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.beam import twiss_iterable_to_df
from ocelot.cpbd.optics import Twiss, twiss as oce_calc_twiss
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.beam import optics_from_moments, moments_from_parray
from ocelot.cpbd.match import match, match_beam


LOG = logging.getLogger(__name__)

ElementAccessType = Optional[Union[int, str, OpticElement]]


class TwissMismatchError(RuntimeError):
    pass


class MachineSequence(Sequence):
    def __init__(self, sequence: list[OpticElement], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._sequence = list(flatten(sequence))

    def __getitem__(self, key: ElementAccessType) -> Union[OpticElement, MachineSequence]:
        if isinstance(key, slice):
            start, step, stop = key.start, key.step, key.stop
            if step is not None:
                raise TypeError("Slice step is not supported")

            if start is None:
                start = 0
            if stop is None:
                stop = -1

            start = self._normalise_key(start)
            stop = self._normalise_key(stop)
            return MachineSequence(self._sequence[start:stop])
        else:
            index = self._normalise_key(key)
            return self._sequence[index]

    def closed_interval(
        self, start: ElementAccessType = None, stop: ElementAccessType = None
    ) -> MachineSequence:
        """Get interval including the stop element at the end (unlike __getitem__)"""
        if start is None:
            start = 0
        istart = self._normalise_key(start)
        if stop is None:
            stop = -1
        istop = self._normalise_key(stop)
        # Closed interval so add 1:
        istop += 1
        return self[istart:istop]

    def __iter__(self):
        yield from iter(self._sequence)

    def _normalise_key(self, key: ElementAccessType) -> int:
        # If already an int then just return it
        if isinstance(key, int):
            if key >= 0:
                return key
            # Convert negative index into a positive one.
            return key % len(self)

        # If a name of an element find the index
        if isinstance(key, str):
            try:
                return self.names().index(key)
            except:
                import ipdb

                ipdb.set_trace()

        # if an element instance then get the key
        try:
            return self._sequence.index(key)
        except (TypeError, ValueError):
            pass

        import ipdb

        ipdb.set_trace()
        raise ValueError(f"Unable to normalise key: {key}")

    def __len__(self) -> int:
        return len(self._sequence)

    def __add__(self, other: Iterable[OpticElement]) -> MachineSequence:
        return type(self)(list(self) + list(flatten(other)))

    def __str__(self) -> str:
        strs = "\n".join([repr(s) for s in self._sequence])
        return f"<{type(self).__name__}:\n{strs}>"

    def __contains__(self, key: ElementAccessType) -> bool:
        try:
            return key in self._sequence or key.id in self.names()
        except AttributeError:
            return key in self.names()

    def names(self) -> list[str]:
        return [x.id for x in self]

    def element_attributes(self, key, property_name: str) -> np.array:
        indices = [i for (i, ele) in enumerate(self) if (ele.id == key) or (ele is key)]
        if not indices:
            raise KeyError
        result = []
        for i in indices:
            result.append(getattr(self[i], property_name))
        return np.squeeze(np.array(result))

    def length(self) -> float:
        return sum(x.l for x in self)

    def element_s(self, key: Union[str, OpticElement]) -> float:
        indices = [i for (i, ele) in enumerate(self) if (ele.id == key) or (ele is key)]

        if not indices:
            raise KeyError

        result = []
        s = 0
        for i, element in enumerate(self):
            length = element.l
            s += length
            if i in indices:
                result.append(s)
        return np.squeeze(np.array(result))

    def __repr__(self) -> str:
        return f"<MachineSequence: {repr(self._sequence)}>"


class FELSection:
    def __init__(self, outdir: os.Pathlike, *args, **kwargs):
        self.outdir = Path(".")
        self.physics_processes_array = []  # list of physics procs with start/stops

        self.particle_dir = self.outdir / "particles"
        self.tws_dir = self.outdir / "tws_dir"
        self.method = {"global": SecondTM}

    @property
    def outdir_particles(self) -> Path:
        return self.outdir / "particles"

    @property
    def outdir_twiss(self) -> Path:
        return self.outdir / "twiss"

    @property
    def output_parray1_file(self) -> Path:
        return self.outdir_particles / f"section_{type(self).__name__}.npz"

    def add_physics_process(self, physics_process, start, stop):
        if not isinstance(start, str) and start is not None:
            raise TypeError(f"{start} should be a str , an element name, or None")
        if not isinstance(stop, str) and stop is not None:
            raise TypeError(f"{stop} should be a str, an element name, or None")

        self.physics_processes_array.append([physics_process, start, stop])

    def save_beam_file(self, parray: ParticleArray) -> None:
        self.output_parray1_file.parent.mkdir(parents=True, exist_ok=True)
        LOG.debug(f"Writing to: {self.output_parray1_file}")
        save_particle_array(self.output_parray1_file, parray)

    # def track(self, parray0, start=None, stop=None, felconfig=None):
    #     navi = self.get_navigator(start=start, stop=stop, felconfig=felconfig)
    #     _, parray1 = track.track(self.lattice, parray0, navi, print_progress=True, overwrite_progress=False)
    #     self.save_beam_file(parray1)
    #     self.post_track(parray1)
    #     return parray1

    def calculate_twiss(
        self,
        twiss0: Twiss,
        idstart: ElementAccessType = None,
        idstop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        sequence = self.get_sequence(
            idstart=idstart, idstop=idstop, felconfig=felconfig
        )
        mlat = MagneticLattice(sequence)
        all_twiss = oce_calc_twiss(mlat, twiss0, return_df=True)
        return all_twiss, mlat

    def get_sequence(self, idstart=None, idstop=None, felconfig=None):
        if felconfig is not None:
            return felconfig.modify_sequence(self.sequence)
        else:
            return deepcopy(self.sequence)

    @classmethod
    @property
    def name(cls):
        return cls.__name__


class SectionedFEL:
    """
    High level class to work with SectionTrack()
    """

    def __init__(
        self, sections: list[FELSection], twiss0: Twiss, outdir: os.PathLike = "./"
    ):
        """

        :param sequence: list of SectionTrack()
        """
        self.outdir = Path(outdir)
        self.sections = [sec(self.outdir) for sec in sections]
        self.twiss0 = twiss0

        self._check_sections_for_duplicate_start_stops()

    def _check_sections_for_duplicate_start_stops(self) -> None:
        sections = self.sections
        previous_section = sections[0]

        for section in sections[1:]:
            if section.sequence[0] == previous_section.sequence[-1]:
                raise ValueError(
                    "Overlapping sequence: {previous_section}: {section}, {previous_stop}"
                )

    def _calculate_twiss_between_two_points(
        self, twiss0: Twiss, start=None, stop=None, felconfig=None
    ):
        # self._check_twiss_position_matches(twiss0, start)
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        mlat = MagneticLattice(sequence)
        full_twiss = oce_calc_twiss(mlat, twiss0, return_df=True)
        return full_twiss, mlat

    def machine_twiss(
        self,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=self.twiss0, stop=stop, felconfig=felconfig
        )
        return twiss, mlat

    def calculate_twiss(
        self, twiss0, start, stop=None, felconfig=None
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        LOG.debug(f'Calculating linear optics between "{start}" and "{stop}"')
        # Final stretch, between the matching point and the provided "stop"
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=twiss0, start=start, stop=stop, felconfig=felconfig
        )

        return twiss, mlat

    def track(
        self,
        parray0,
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
        dumps=None,
    ) -> ParticleArray:
        parray1, _ = self._piecewise_matched_tracking(
            parray0,
            calculate_optics=False,
            start=start,
            stop=stop,
            felconfig=felconfig,
            dumps=dumps,
        )

        return parray1

    def track_optics(
            self,
            parray0: ParticleArray,
            start: ElementAccessType = None,
            stop: ElementAccessType = None,
            felconfig: Optional[FELSimulationConfig] = None,
            dumps=None,
            physics=False,
            opticscls=None,
            outdir: os.PathLike = None,
    ) -> tuple[ParticleArray, pd.DataFrame]:
        parray1, twiss = self._piecewise_matched_tracking(
            parray0,
            start=start,
            calculate_optics=True,
            stop=stop,
            physics=physics,
            felconfig=felconfig,
            opticscls=opticscls,
            outdir=outdir,
            dumps=dumps,
        )
        return parray1, twiss

    def _track_between_two_points(
        self,
        parray0: ParticleArray,
        *,
        calculate_optics: bool,
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        outdir: os.PathLike = None,
        dumps=None,
            physics=False,
        opticscls=None,
        felconfig: Optional[FELSimulationConfig] = None,
    ) -> tuple[ParticleArray, pd.DataFrame]:
        if outdir is None:
            outdir = Path(self.outdir)

        navi = self.to_navigator(start=start, stop=stop, felconfig=felconfig, physics=physics)
        if calculate_optics:
            calc_procs = _add_markers_to_navi_for_optics(navi, opticscls)
        if dumps:
            _add_save_beams(navi, dumps, outdir)

        LOG.debug(f'Starting tracking between "{start}" and "{stop}"')
        _, parray1 = track.track(
            navi.lat,
            parray0,
            navi=navi,
            calc_tws=False,
            print_progress=True,
            overwrite_progress=False,
        )
        if calculate_optics:
            twiss = _extract_twiss_paramters_from_twiss_processes(calc_procs)
        else:
            twiss = []

        return parray1, twiss

    def _piecewise_matched_tracking(
        self,
        parray0: ParticleArray,
        calculate_optics: bool,
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
            physics: bool = True,
        outdir: os.PathLike = None,
        opticscls=None,
        dumps=None,
    ):
        LOG.debug("Beginning tracking")
        parray0 = parray0.copy()
        tracked_twiss = []
        # Track to the matching points.  Which are also apparently assumed to be in order...
        parray1, section_twiss = self._track_between_two_points(
            parray0,
            calculate_optics=calculate_optics,
            start=start,
            stop=stop,
            outdir=outdir,
            physics=physics,
            dumps=dumps,
            opticscls=opticscls,
            felconfig=felconfig,
        )
        if calculate_optics:
            tracked_twiss.append(section_twiss)
            if tracked_twiss:
                full_twiss_df = pd.concat(tracked_twiss)
        else:
            full_twiss_df = pd.DataFrame()

        return parray1, full_twiss_df

    def to_navigator(
        self,
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
            physics: bool = True,
    ) -> Navigator:
        if felconfig is None:
            felconfig = FELSimulationConfig()

        # Make the navigator
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        mlat = MagneticLattice(sequence, method={"global": SecondTM})
        navi = Navigator(mlat, unit_step=felconfig.unit_step)

        sequence_before = self._get_sequence_before_start(start)
        sequence_after = self._get_sequence_after_stop(stop)

        if not physics:
            return navi

        # Now attach physics to the navigator we made above.
        LOG.info(f'Adding procs 2 Navi:  Sequence start = "{start}", stop = "{stop}"')
        for section in self.sections:
            for (
                process,
                proc_start_name,
                proc_stop_name,
            ) in section.physics_processes_array:
                assert isinstance(proc_start_name, str) or proc_start_name is None
                assert isinstance(proc_stop_name, str) or proc_stop_name is None

                # Special values meaning at the start/end of the section.
                if proc_start_name is None:
                    proc_start = None
                if proc_stop_name is None:
                    proc_stop = None

                # If process stops before the sequene even begins, then skip.
                if proc_stop_name in sequence_before:
                    LOG.debug(
                        f"From Section {type(self).__name__} - Skipping {process}: start={proc_start_name}, stop={proc_stop_name}. Reason: Process precedes sequence"
                    )
                    continue

                # If process doesn't even start until after the sequence ends, then skip
                if proc_start_name in sequence_after:
                    LOG.debug(
                        f"From Section {type(self).__name__} - Skipping {process}: start={proc_start_name}, stop={proc_stop_name}. Reason: Process follows sequence"
                    )
                    continue

                # If the process starts before, then
                if proc_start_name in sequence_before:
                    proc_start = sequence[0]
                elif proc_start_name is not None:
                    proc_start = sequence[proc_start_name]

                if proc_stop_name in sequence_after:
                    proc_stop = sequence[-1]
                elif proc_stop_name is not None:
                    proc_stop = sequence[proc_stop_name]

                # Set energy for CSR processes incase it needs it.
                if isinstance(process, CSR):
                    self._set_csr_process_energy(
                        process, proc_start_name, proc_stop_name, felconfig=felconfig
                    )

                try:
                    navi.add_physics_proc(process, proc_start, proc_stop)
                except:
                    import ipdb

                    ipdb.set_trace()

                print(
                    f"From Section {type(section).__name__} - Adding {process}: start={proc_start_name}, stop={proc_stop_name}"
                )

        return navi

    def _get_sequence_before_start(self, start: ElementAccessType) -> MachineSequence:
        if start is None:
            return MachineSequence([])
        else:
            return self.get_sequence(stop=start)

    def _get_sequence_after_stop(self, stop: ElementAccessType) -> MachineSequence:
        if stop is None:
            return MachineSequence([])
        else:
            return self.get_sequence(start=stop)

    def _set_csr_process_energy(
        self,
        process: PhysProc,
        start_name: ElementAccessType,
        stop_name: ElementAccessType,
        felconfig: FELSimulationConfig,
    ) -> None:
        machine_twiss, _ = self.machine_twiss(felconfig=felconfig)
        start_energy = machine_twiss.query(f"id == '{start_name}'").iloc[0].E
        stop_energy = machine_twiss.query(f"id == '{stop_name}'").iloc[0].E
        assert start_energy == stop_energy
        process.energy = start_energy

    def get_sequence(
        self,
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: FELSimulationConfig = None,
    ) -> MachineSequence:
        result = []

        found_start = False
        found_stop = False

        full_sequence = start is None and stop is None
        section_sequences = [x.get_sequence() for x in self.sections]

        # There is a difference in behaviour between stop is None and
        # stop is the last element of the sequence.  only in the latter
        # there shouldn't be...
        xstart = start
        xstop = stop

        if start is None:
            start = section_sequences[0][0].id
        if stop is None:
            stop = section_sequences[-1][-1].id

        for sequence in section_sequences:
            this_section_start = None
            this_section_stop = None
            if start in sequence:
                this_section_start = start
                found_start = True
            if stop in sequence:
                this_section_stop = stop
                found_stop = True

            if not found_start and not full_sequence:
                continue

            if this_section_stop is not None:
                result.extend(
                    sequence.closed_interval(this_section_start, this_section_stop)
                )
            else:
                result.extend(sequence.closed_interval(this_section_start, stop=None))
            if found_stop and not full_sequence:
                break

        if start is not None and not found_start:
            raise ValueError(f"Start position: {start} not found in sequence")
        if stop is not None and not found_stop:
            raise ValueError(f"Stop position {stop} not found in sequence")

        if felconfig is not None:
            return MachineSequence(felconfig.controller.modify_sequence(result))
        return MachineSequence(deepcopy(result))

    def length(self) -> float:
        x = 0
        for section in self.sections:
            x += section.sequence.length()
        return x

    @property
    def section_names(self) -> list[str]:
        return [sec.name for sec in self.sections]

    def __getitem__(self, key: str) -> FELSection:
        return dict(zip(self.section_names, self.sections))[key]

    def get_element(self, name: str) -> OpticElement:
        for section in self.sections:
            for element in section.sequence:
                if element.id == name:
                    return element
        raise KeyError(f"{name} not found")

    def get_element_end_s(self, name: str) -> float:
        s = 0
        for element in self.get_sequence():
            s += element.l
            if element.id == name:
                return s
        raise KeyError(f"{name} not found")

    def get_element_start_s(self, name: str) -> float:
        s = 0
        for element in self.get_sequence():
            if element.id == name:
                return s
            s += element.l
        raise KeyError(f"{name} not found")

    def match_beam(self, parray0: ParticleArray, *, 
                   elements: list[str],
                   constraints: dict[str: dict[str, float]],
                   start: ElementAccessType = None,
                   stop: ElementAccessType = None,
                   felconfig: Optional[FELSimulationConfig] = None,
                   **match_beam_kwargs) -> FELSimulationConfig:

        # Make navigator from given start, stop and felconfig
        navi = self.to_navigator(start=start, stop=stop, felconfig=felconfig)

        # Update match_beam_kwargs with function **kwargs.
        match_beam_kwargs = {"verbose": True, "max_iter": 500, "method": "simplex", "min_i5": False,
                             "bounds": None, "vary_bend_angle": False} | match_beam_kwargs

        # Turn ElementAccessType into Element instances, if they are
        # not already, for elements and constraints arguments.
        
        # Turn element names into element instances in the the Navigator
        elements = _get_element_instances_from_mlat(navi.lat, elements)
        # Magic one liner that takes element names as keys and turns them into Element instances
        constraints = dict(zip(_get_element_instances_from_mlat(navi.lat, constraints.keys()), constraints.values()))
        # Do the matching

        res = match_beam(navi.lat, constraints, elements, parray0, navi, **match_beam_kwargs)

        # Apply the result to the simulation config.
        if felconfig is None:
            new_conf = FELSimulationConfig()
        else:
            new_conf = deepcopy(felconfig)

        new_conf.controller.components = {quad_name: {"k1": strength} for (quad_name, strength)
                                          in zip([x.id for x in elements], res)}

        return new_conf

    def match(self, twiss0: Twiss, *, 
              elements: list[str],
              constraints: dict[str: dict[str, float]],
              start: ElementAccessType = None,
              stop: ElementAccessType = None,
              felconfig: Optional[FELSimulationConfig] = None,
              **match_kwargs) -> FELSimulationConfig:

        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        lat_matching = MagneticLattice(sequence)

        elements = _get_element_instances_from_mlat(lat_matching, elements)
        # Magic one liner that takes element names as keys and turns them into Element instances
        constraints = dict(zip(_get_element_instances_from_mlat(lat_matching, constraints.keys()), constraints.values()))

        res = match(lat_matching, constraints, elements, twiss0, **match_kwargs)

        # Apply the result to the simulation config.
        if felconfig is None:
            new_conf = FELSimulationConfig()
        else:
            new_conf = deepcopy(felconfig)
        
        new_conf.controller.components = {quad_name: {"k1": strength} for (quad_name, strength)
                                          in zip([x.id for x in elements], res)}

        return new_conf

        

class ElementNotFoundError(KeyError):
    pass


class SequenceController:
    def __init__(self, regex: str):
        self.regex = re.compile(regex)

    def match(self, element: OpticElement) -> bool:
        return re.match(self.regex, element.id)


class TDSControl(SequenceController):
    def __init__(
        self, regex: str, phi: Optional[float] = None, v: Optional[float] = None
    ):
        super().__init__(regex)
        self.phi = None
        self.v = None
        self.active = True

    @property
    def should_modify(self) -> bool:
        return not self.active or self.phi is not None or self.v is not None

    def apply(self, element: OpticElement) -> None:
        if not self.should_modify:
            return

        if not self.active:
            self.disable(element)
            return

        print(
            f"Applying Cavity control to element: {element}.  phi1 = {self.phi}, v1 = {self.v}"
        )
        if self.v is not None:
            element.v = self.v
        if self.phi is not None:
            element.phi = self.phi

    def disable(self, element: OpticElement) -> None:
        element.v = 0.0
        print(
            f"Applying Cavity control to element: {element}.  Disabling: {element.v=}"
        )

    def __repr__(self) -> str:
        tname = type(self).__name__
        return f"<{tname}: phi={self.phi}, v={self.v}>"


class ChicaneControl(SequenceController):
    def __init__(self, regex: str, rho: Optional[float] = None):
        super().__init__(regex)
        self.rho = None

    def apply(self, sequence: MachineSequence) -> None:
        if rho is None:
            return

    def _get_chicane(self, sequence: MachineSequence) -> None:
        dipoles = []
        for element in sequence:
            if self.match(element):
                dipoles.append(element)


class CavityControl(TDSControl):
    def __init__(
        self,
        regex: str,
        phi: Optional[float] = None,
        v: Optional[float] = None,
        coupler_kick: bool = False,
    ):
        super().__init__(regex, phi=phi, v=v)
        self.coupler_kick = coupler_kick

    def apply(self, element: OpticElement) -> None:
        super().apply(element)
        if not self.coupler_kick:
            print(f"Removing coupler kick: {repr(element)}")
            element.remove_coupler_kick()


class EuXFELController:
    def __init__(self, components: dict = None):
        self.components = components if components is not None else {}
        self.tds1 = TDSControl(regex=r"TDSA\.52\.I1")
        self.a1 = CavityControl(r"C\.A1\.1\.[0-8].I1")
        self.ah1 = CavityControl(r"C3\.AH1\.1\.[0-8].I1")
        self.a2 = CavityControl(regex=r"C\.A2\.[0-4]\.[0-8].L1")
        self.a3 = CavityControl(regex=r"C\.A[3-5]\.[0-4]\.[0-8].L2")
        self.bc0 = ChicaneControl(regex=r"BB\.(96|98|100|101)\.I1")

        self.global_coupler_kick = None

    def modify_sequence(self, sequence: Iterable[OpticElement]) -> MachineSequence:
        """Given an EuXFEL sequence, modify it based on the
        configuration of this EuXFELController instance."""
        new_sequence = deepcopy(sequence)
        for element in new_sequence:
            # Set TDS1
            if self.tds1.match(element):
                self.tds1.apply(element)
                continue

            # Set A1
            if self.a1.match(element):
                self.a1.apply(element)
                continue

            if self.ah1.match(element):
                self.ah1.apply(element)
                continue

            if self.a2.match(element):
                self.a2.apply(element)
                continue

            if self.a3.match(element):
                self.a3.apply(element)
                continue

        new_sequence = self.apply_magnet_conf_to_sequence(new_sequence)

        return MachineSequence(new_sequence)

    def apply_magnet_conf_to_sequence(
        self, sequence: Iterable[OpticElement]
    ) -> MachineSequence:
        sequence = deepcopy(sequence)
        magnet_conf = self.components
        for element in sequence:
            name = element.id
            if name in magnet_conf:
                for key, value in magnet_conf[name].items():
                    print(f"Setting: {name}->{key} to {value}")
                    setattr(element, key, value)

        return MachineSequence(sequence)


class FELSimulationConfig:
    def __init__(
        self,
        metadata: Optional[dict] = None,
        coupler_kick: bool = False,
        controller: Optional[EuXFELController] = None,
    ):
        self.metadata = metadata if metadata is not None else {}
        self.coupler_kick = False

        self.controller = EuXFELController() if controller is None else controller

        self.unit_step = 0.02

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""\
        FELSimulationConfig:
        - TDS1: phi = {self.controller.tds1.phi}, v = {self.controller.tds1.v}
        - A1 Cavities: phi = {self.controller.a1.phi}, v = {self.controller.a1.v}
        - AH1 Cavities: phi = {self.controller.ah1.phi}, v = {self.controller.ah1.v}"""
        )


    # def __ior__(self, other) -> FELSimulationConfig:
    #     # self.metadata |= other.metadata
    #     self.controller.components |= other.controller.components
    #     self.coupler_kick = other.coupler_kick
    #     self.unit_step = other.unit_step

    # def __or__(self, other):
    #     result = deepcopy(self)
    #     result |= other
    #     return result

    # def bc_analysis(self):
    #     # find positions
    #     L = 0.
    #     for elem in self.sequence:

    #         if elem in self.dipoles:
    #             elem.s_pos = L
    #         L += elem.l
    #     self.dipoles = np.array(self.dipoles)
    #     # sorting - put dipoles in right order
    #     sort_index = np.argsort(np.array([d.s_pos for d in self.dipoles]))
    #     #print(sort_index)
    #     self.dipoles = self.dipoles[sort_index]

    #     #self.bc_gap = self.dipoles[2].s_pos - self.dipoles[1].s_pos - self.dipoles[1].l
    #     #print("bc_gap", self.bc_gap)
    #     self.get_bc_shoulders()

    # def get_bc_shoulders(self):
    #     self.left_shoulder = []
    #     self.right_shoulder = []
    #     left_flag = False
    #     right_flag = False
    #     for i, elem in enumerate(self.sequence):
    #         if elem == self.dipoles[0]:
    #             left_flag = True
    #         elif elem == self.dipoles[1]:
    #             left_flag = False

    #         if elem == self.dipoles[2]:
    #             right_flag = True
    #         elif elem == self.dipoles[3]:
    #             right_flag = False

    #         if left_flag and elem != self.dipoles[0]:
    #             self.sequence[i] = copy.deepcopy(elem)
    #             self.left_shoulder.append(self.sequence[i])
    #         if right_flag and elem != self.dipoles[2]:
    #             self.sequence[i] = copy.deepcopy(elem)
    #             self.right_shoulder.append(self.sequence[i])

    #     self.bc_gap_left = np.sum([d.l for d in self.left_shoulder])
    #     right_len = np.sum([d.l for d in self.right_shoulder])

    #     for d in self.left_shoulder:
    #         d.len_coef = d.l / self.bc_gap_left
    #     for d in self.right_shoulder:
    #         d.len_coef = d.l / right_len

    # def change_bc_shoulders(self, drift):
    #     for d in self.left_shoulder:
    #         d.l = drift * d.len_coef

    #     for d in self.right_shoulder:
    #         d.l = drift * d.len_coef

    # def update_bunch_compressor(self, rho):
    #     if self.dipoles is None:
    #         print(self.__class__.__name__ + " No BC")
    #         return

    #     self.bc_analysis()

    #     if self.dipole_len is None:
    #         self.dipole_len = copy.copy(self.dipoles[0].l)

    #     if rho == 0:
    #         angle = 0
    #         ds = self.dipole_len
    #     else:
    #         angle = np.arcsin(self.dipole_len / rho)
    #         ds = angle * rho

    #     if self.bc_gap is None:
    #         self.bc_gap = self.bc_gap_left*np.cos(self.dipoles[0].angle)

    #     drift = self.bc_gap / np.cos(angle)
    #     self.change_bc_shoulders(drift)
    #     # d.l=drift
    #     for i, dip in enumerate(self.dipoles):
    #         dip.angle = angle * np.sign(dip.angle)
    #         dip.l = ds
    #         if i in [0, 2]:
    #             dip.e2 = angle * np.sign(dip.angle)
    #         else:
    #             dip.e1 = angle * np.sign(dip.angle)


class NoTwiss(PhysProc):
    MATCH = None

    def __init__(self, name):
        super().__init__()
        self.name = name
        self.twiss = None

    def __str__(self):
        return f"<{type(self).__name__}, name={self.name}>"


class TwissCalculator(PhysProc):
    MATCH = None
    """Calculate dispersion-free Twiss parameters"""

    def __init__(self, name):
        super().__init__()
        self.name = name
        self.twiss = None

    def apply(self, parray, dz):
        mean, cov = moments_from_parray(parray)
        self.twiss = optics_from_moments(mean, cov, parray.E)
        self.twiss.s = float(parray.s)

    def __str__(self):
        return f"<{type(self).__name__}, name={self.name}>"


class SliceTwissCalculator(PhysProc):
    MATCH = "Emax"
    """Calculate dispersion-free Twiss parameters"""

    def __init__(self, name):
        super().__init__()
        self.name = name
        self.twiss = None

    def apply(self, parray, dz):
        # mean, cov = moments_from_parray(parray)
        # self.twiss = optics_from_moments(mean, cov, parray.E)
        # self.twiss.s = parray.s
        # self.twiss.id = sel
        # print("!!!!!!!!!!!!!!!!!!!!!!!!", self.name)
        # if self.name == "MATCH.428.B2_after":
        #     import ipdb; ipdb.set_trace()

        self.twiss = twiss_parray_slice(
            parray,
            slice="Emax",
            nparts_in_slice=5000,
            smooth_param=0.05,
            filter_base=2,
            filter_iter=2,
        )
        self.twiss.s = float(parray.s)


def _add_markers_to_navi_for_optics(navi: Navigator, opticscls: Type = None):
    # Put markers everywhere.
    _, new_markers = insert_markers_by_predicate(navi.lat.sequence, lambda ele: True)

    if opticscls is None:
        opticscls = TwissCalculator

    twiss_processes = [opticscls(marker.id) for marker in new_markers]
    navi.add_physics_processes(twiss_processes, new_markers, new_markers)

    for proc in navi.process_table.proc_list:
        elem0 = proc.start_elem
        proc.indx0 = navi.lat.sequence.index(elem0)
        elem1 = proc.end_elem
        proc.indx1 = navi.lat.sequence.index(elem1)

    # Have to go by name for the reference table because it is by
    # definition a deepcopy of the original table (i.e. including the
    # elem0/elem1) so looking up by identity won't work.
    element_names = [ele.id for ele in navi.lat.sequence]
    for proc in navi.ref_process_table.proc_list:
        elem0 = proc.start_elem
        proc.indx0 = element_names.index(elem0.id)
        elem1 = proc.end_elem
        proc.indx1 = element_names.index(elem1.id)

    return twiss_processes


def _add_save_beams(navi, element_names, outdir):
    for element_name in element_names:
        fname = Path(outdir) / f"{element_name}.npz"
        proc = SaveBeam(fname)
        proc.element_name = element_name
        try:
            navi.add_physics_proc(proc, element_name, element_name)
        except KeyError:
            continue
        LOG.debug(f"Added SaveBeam process to {element_name}, writing to {fname}")


def _extract_twiss_paramters_from_twiss_processes(twiss_processes):
    # We filter Nones just because at the very end of the lattice
    # some markers are not tracked through due the tolerance check
    # in track.track (markers have 0. length).  Sort also by s.
    twiss_instances = []
    for proc in twiss_processes:
        twiss = proc.twiss
        if twiss is None:
            continue
        twiss_instances.append(twiss)
    twiss_instances.sort(key=lambda tws: tws.s)

    twissdf = twiss_iterable_to_df(twiss_instances)

    return twissdf

def _get_element_instances_from_mlat(mlat, names):
    elements = []
    def f(ele):
        return ele.id in names
    indices = mlat.find_indices_by_predicate(f)
    return [mlat.sequence[i] for i in indices]

# def _get_element_instance_from_mlat(mlat, name):
#     def f(ele):
#         return ele.id in names
#     indices = mlat.find_indices_by_predicate(f)
#     assert len(indices
#     return mlat.sequence[i]
    
