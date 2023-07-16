"""
Section class for s2e tracking.
S.Tomin. XFEL/DESY. 2017
"""
from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Optional, Union, Type, Iterable, TypeVar
from pathlib import Path
import logging
from copy import deepcopy
import textwrap
import re
import pandas as pd
from dataclasses import dataclass
import functools

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
from ocelot.cpbd.beam import twiss_parray_slice

from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.beam import twiss_iterable_to_df
from ocelot.cpbd.optics import Twiss, twiss as oce_calc_twiss
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.beam import optics_from_moments, moments_from_parray
from ocelot.cpbd.match import match, match_beam


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.INFO)


ElementT = TypeVar("ElementT", bound=OpticElement)
ElementAccessType = Optional[Union[int, str, ElementT]]
ElementSequenceT = list[ElementT]


class TwissMismatchError(RuntimeError):
    pass


class MachineSequence(Sequence):
    def __init__(self, sequence: list[ElementT], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._sequence = list(flatten(sequence))

    def __getitem__(self, key: ElementAccessType) -> Union[ElementT, MachineSequence]:
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
            return self.names().index(key)

        # if an element instance then get the key
        try:
            return self._sequence.index(key)
        except (TypeError, ValueError):
            pass

        raise ValueError(f"Unable to normalise key: {key}")

    def __len__(self) -> int:
        return len(self._sequence)

    def __add__(self, other: Iterable[ElementT]) -> MachineSequence:
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

    def element_s(self, key: Union[str, ElementT]) -> float:
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
    def __init__(self, sequence: MachineSequence):
        self.sequence = MachineSequence(sequence)
        self.physics_processes = []  # list of physics procs with start/stops
        self.setup_physics()

    def setup_physics(self):
        pass

    def add_physics_process(self, physics_process, start, stop):
        if not isinstance(start, str) and start is not None:
            raise TypeError(f"{start} should be a str , an element name, or None")
        if not isinstance(stop, str) and stop is not None:
            raise TypeError(f"{stop} should be a str, an element name, or None")

        self.physics_processes.append([physics_process, start, stop])

    def calculate_twiss(
        self,
        twiss0: Twiss,
        idstart: ElementAccessType = None,
        idstop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        sequence = self.get_sequence(
            idstart=idstart, idstop=idstop, felconfig=felconfig_net
        )
        mlat = MagneticLattice(sequence)
        all_twiss = oce_calc_twiss(mlat, twiss0, return_df=True)
        return all_twiss, mlat

    def get_sequence(self, idstart=None, idstop=None, felconfig=None):
        if felconfig is not None:
            return felconfig.new_sequence(self.sequence)
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
        self,
        sections: list[FELSection],
        twiss0: Twiss,
        outdir: os.PathLike = "./",
        felconfig=None,
    ):
        """

        :param sequence: list of SectionTrack()
        """
        self.outdir = Path(outdir)
        self.sections = sections
        self.twiss0 = twiss0
        self.felconfig = felconfig if felconfig else None

        self._check_sections_for_duplicate_start_stops()

    def _net_felconfig(self, felconfig2=None):
        felconfig = self.felconfig
        if felconfig is None:
            felconfig = FELSimulationConfig()

        if felconfig2 is None:
            felconfig2 = FELSimulationConfig()

        return felconfig | felconfig2

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

    def design_twiss(self, stop: ElementAccessType = None):
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=self.twiss0, stop=stop
        )
        return twiss, mlat

    def machine_twiss(
        self,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        felconfig = self._net_felconfig(felconfig)
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=self.twiss0, stop=stop, felconfig=felconfig
        )
        return twiss, mlat

    def calculate_twiss(
        self, twiss0, start, stop=None, felconfig=None
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        felconfig_net = self._net_felconfig(felconfig)
        LOG.debug(f'Calculating linear optics between "{start}" and "{stop}"')
        # Final stretch, between the matching point and the provided "stop"
        twiss0 = _coerce_to_ocelot_twiss(twiss0)
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
        physics=True,
    ) -> ParticleArray:
        parray1, _ = self._piecewise_matched_tracking(
            parray0,
            calculate_optics=False,
            start=start,
            stop=stop,
            felconfig=felconfig,
            dumps=dumps,
            physics=physics,
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

        navi = self.to_navigator(
            start=start, stop=stop, felconfig=felconfig, physics=physics
        )
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
        felconfig = self._net_felconfig(felconfig)

        # Make the navigator
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        mlat = MagneticLattice(sequence, method={"global": felconfig.tracking.method})
        navi = Navigator(mlat, unit_step=felconfig.tracking.unit_step)

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
            ) in section.physics_processes:
                assert isinstance(proc_start_name, str) or proc_start_name is None
                assert isinstance(proc_stop_name, str) or proc_stop_name is None

                # Special values meaning at the start/end of the section.
                if proc_start_name is None:
                    proc_start = None
                if proc_stop_name is None:
                    proc_stop = None

                if (
                    proc_start_name == proc_stop_name and process in sequence
                ):  # For occasion when you want to attach an element to the first
                    pass
                # If process stops before the sequene even begins, then skip.
                elif proc_stop_name in sequence_before:
                    LOG.debug(
                        f"From Section {type(self).__name__} - Skipping {process}: start={proc_start_name}, stop={proc_stop_name}. Reason: Process precedes sequence"
                    )
                    continue

                # If process doesn't even start until after the sequence ends, then skip
                elif proc_start_name in sequence_after:
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

                navi.add_physics_proc(process, proc_start, proc_stop)

                # LOG.info(
                #     f"From Section {type(section).__name__} - Adding {process}: start={proc_start_name}, stop={proc_stop_name}"
                # )

                LOG.info(
                    f"From Section {type(section).__name__} - Adding {process}: start={proc_start_name}, stop={proc_stop_name}:"
                )
                LOG.info(repr(process))

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

    def _get_csr_process_sigma_min(self):
        pass

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

        net_felconfig = self._net_felconfig(felconfig2=felconfig)
        return MachineSequence(net_felconfig.new_sequence(result))

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

    def get_element(self, name: str) -> ElementT:
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

    def match_beam(
        self,
        parray0: ParticleArray,
        *,
        elements: list[str],
        constraints: dict[str : dict[str, float]],
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
        physics=True,
        match="projected",
        **match_beam_kwargs,
    ) -> FELSimulationConfig:
        # Make navigator from given start, stop and felconfig
        navi = self.to_navigator(
            start=start, stop=stop, felconfig=felconfig, physics=physics
        )

        if match == "projected":
            twiss_function = None
        elif match.lower() == "emax":
            twiss_function = partial.functools(_twiss_central_slice, match_slice="EMax")
        elif match.lower() == "imax":
            twiss_function = partial.functools(_twiss_central_slice, match_slice="IMax")
        else:
            raise ValueError(f"Unknown match string: {match}")

        # Update match_beam_kwargs with function **kwargs.
        match_beam_kwargs = {
            "verbose": True,
            "max_iter": 1_000,
            "method": "simplex",
            "min_i5": False,
            "bounds": None,
            "vary_bend_angle": False,
            "tws_fn": twiss_function,
        } | match_beam_kwargs

        # Turn ElementAccessType into Element instances, if they are
        # not already, for elements and constraints arguments.

        # Turn element names into element instances in the the Navigator
        elements = _get_element_instances_from_mlat(navi.lat, elements)
        # Magic one liner that takes element names as keys and turns them into Element instances
        constraints = dict(
            zip(
                _get_element_instances_from_mlat(navi.lat, constraints.keys()),
                constraints.values(),
            )
        )
        # Do the matching

        res = match_beam(
            navi.lat, constraints, elements, parray0, navi, **match_beam_kwargs
        )

        # Apply the result to the simulation config.
        new_conf = deepcopy(self._net_felconfig(felconfig))

        new_conf.components.update(
            {
                quad_name: {"k1": strength}
                for (quad_name, strength) in zip([x.id for x in elements], res)
            }
        )

        return new_conf

    def match(
        self,
        twiss0: Twiss,
        *,
        elements: list[str],
        constraints: dict[str : dict[str, float]],
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: Optional[FELSimulationConfig] = None,
        verbose=True,
        **match_kwargs,
    ) -> FELSimulationConfig:
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        lat_matching = MagneticLattice(sequence)

        elements = _get_element_instances_from_mlat(lat_matching, elements)
        # Magic one liner that takes element names as keys and turns them into Element instances
        constraints = dict(
            zip(
                _get_element_instances_from_mlat(lat_matching, constraints.keys()),
                constraints.values(),
            )
        )

        res = match(
            lat_matching, constraints, elements, twiss0, verbose=verbose, **match_kwargs
        )

        # Apply the result to the simulation config.
        new_conf = deepcopy(self._net_felconfig(felconfig))

        matched_strengths = {
            quad_name: {"k1": strength}
            for (quad_name, strength) in zip([x.id for x in elements], res)
        }

        new_conf.components.update(matched_strengths)

        return new_conf


class ElementNotFoundError(KeyError):
    pass


class Controller:
    def __init__(self, regex: str):
        self.regex = re.compile(regex)

    def match(self, element: ElementT) -> bool:
        return re.match(self.regex, element.id)

    def should_modify(self) -> bool:
        pass

    def modify_element(self, element):
        pass

    def modify_sequence(self, sequence):
        pass


class CavityController(Controller):
    def __init__(
        self,
        regex: str,
        phi: Optional[float] = None,
        v: Optional[float] = None,
        active: bool = None,
        coupler_kick=None,
    ):
        super().__init__(regex)
        self.phi = phi
        self.v = v
        self.active = active
        self.coupler_kick = coupler_kick

    def should_modify(self) -> bool:
        return self.active is not None or self.phi is not None or self.v is not None

    def modify_element(self, element: ElementT) -> None:
        if not self.should_modify():
            return

        if self.active is False:
            self.disable(element)
            return

        if not self.coupler_kick:
            try:
                element.remove_coupler_kick()
                LOG.debug(f"Removing coupler kick: {repr(element)}")
            except AttributeError:
                pass

        LOG.debug(
            f"Applying Cavity control to element: {element}.  phi1 = {self.phi}, v1 = {self.v}"
        )
        if self.v is not None:
            element.v = self.v
        if self.phi is not None:
            element.phi = self.phi

    def disable(self, element):
        LOG.debug(f"Disabling: {element}")
        element.v = 0.0

    def __repr__(self) -> str:
        tname = type(self).__name__
        return f"<{tname}: phi={self.phi}, v={self.v}, active={self.active}>"

    def __or__(self, other):
        result = deepcopy(self)
        result |= other
        return result

    def __ior__(self, other):
        self.phi = _set_if_not_none(self.phi, other.phi)
        self.v = _set_if_not_none(self.v, other.v)
        self.active = _set_if_not_none(self.active, other.active)
        return self


class ChicaneController(Controller):
    MAX_DIPOLES = 4
    DEFAULT_BEND_POLARITY = -1
    # TODO: Does this work with sequences that are e.g. only one dipole?
    # Two/3 dipoles instead of full chicane?

    def __init__(self, regex, angle=None, r56=None):
        super().__init__(regex)
        self.angle = angle
        self.r56 = r56

        if self.angle is not None and self.r56 is not None:
            raise TypeError("angle and r56 are mutually exclusive")

    def _find_dipole_indices(self, sequence) -> list[int]:
        dipole_indices = []
        for i, element in enumerate(sequence):
            if self.regex.match(element.id):
                dipole_indices.append(i)
            if len(dipole_indices) == self.MAX_DIPOLES:
                break
        return sorted(dipole_indices)

    def should_modify(self) -> bool:
        return self.angle is not None or self.r56 is not None

    def _get_new_angle(self, sequence, dipole_indices):
        projected_outer_drift_length, *_ = self._get_projected_lengths_from_sequence(
            sequence, dipole_indices
        )

        if self.r56 is not None:
            dipole_length = sequence[dipole_indices[0]].l
            l0 = projected_outer_drift_length
            angle = np.sqrt((-self.r56 * 0.5) / (l0 + 2 * dipole_length / 3))
        else:
            angle = abs(self.angle)
        assert angle is not None

        return angle

    def _get_projected_lengths_from_sequence(
        self, sequence: ElementSequenceT, dipole_indices: list[int]
    ):
        first_drift = sequence[dipole_indices[0] + 1 : dipole_indices[1]]
        second_drift = sequence[dipole_indices[1] + 1 : dipole_indices[2]]

        first_drift_length = sum(x.l for x in first_drift)
        second_drift_length = sum(x.l for x in second_drift)

        first_dipole = sequence[dipole_indices[0]]
        first_projected_length = first_drift_length * np.cos(first_dipole.angle)

        return first_projected_length, second_drift_length, first_projected_length

    @staticmethod
    def _assert_consistency_in_dipole_angles(dipoles):
        angle = abs(dipoles[0].angle)
        for dipole in dipoles:
            if abs(dipole.angle) != angle:
                raise ValueError("Chicane has inconsistent dipole angles")

    def _update_angles(self, sequence, dipole_indices):
        # Get angle from attributes
        angle = abs(self._get_new_angle(sequence, dipole_indices))
        # Get dipoles whose angles we should update
        dipoles = [sequence[i] for i in sorted(dipole_indices)]
        # Get first dipole to get what polarity (-++-) or (+--+) we should use
        first_dipole = dipoles[0]
        bend_polarity = np.sign(first_dipole.angle)
        if not bend_polarity:
            bend_polarity = self.DEFAULT_BEND_POLARITY

        # Update the angles
        dipoles[0].angle = bend_polarity * angle
        dipoles[0].e1 = 0.0
        dipoles[0].e2 = bend_polarity * angle

        dipoles[1].angle = -bend_polarity * angle
        dipoles[1].e1 = -bend_polarity * angle
        dipoles[1].e2 = 0

        dipoles[2].angle = -bend_polarity * angle
        dipoles[2].e1 = 0.0
        dipoles[2].e2 = -bend_polarity * angle

        dipoles[3].angle = bend_polarity * angle
        dipoles[3].e1 = bend_polarity * angle
        dipoles[3].e2 = 0.0

    @staticmethod
    def _update_drifts(sequence, dipole_indices):
        # Two outer drifts of chicane that we need to modify
        # could be multiple elements because of markers etc...
        first_drift = sequence[dipole_indices[0] + 1 : dipole_indices[1]]
        third_drift = sequence[dipole_indices[2] + 1 : dipole_indices[3]]

        first_drift_length = sum(x.l for x in first_drift)
        third_drift_length = sum(x.l for x in third_drift)
        assert np.isclose(first_drift_length, third_drift_length, atol=1e-6)

        angle = abs(sequence[dipole_indices[0]].angle)
        projected_length = first_drift_length * np.cos(angle)
        new_drift_length = projected_length / np.cos(angle)

        ratio = new_drift_length / first_drift_length

        # finally update lengths.
        for element in first_drift + third_drift:
            element.l *= ratio

    def _assert_consistency_in_drift_lengths(self, sequence, dipole_indices):
        first_drift = sequence[dipole_indices[0] + 1 : dipole_indices[1]]
        third_drift = sequence[dipole_indices[2] + 1 : dipole_indices[3]]

        first_drift_length = sum(x.l for x in first_drift)
        third_drift_length = sum(x.l for x in third_drift)

        if not np.isclose(first_drift_length, third_drift_length, atol=1e-6):
            raise ValueError(
                "Outer chicane drifts have different lengths: {first_drift_length=}, {third_drift_length}="
            )

    def _assert_consistency_in_chicane(self, sequence, dipole_indices):
        dipoles = [sequence[i] for i in dipole_indices]
        # check all angles are the same in magnitude
        self._assert_consistency_in_dipole_angles(dipoles)
        # Check outer drifts have equal lengths
        self._assert_consistency_in_drift_lengths(sequence, dipole_indices)

    def modify_sequence(self, sequence):
        if not self.should_modify():
            return

        dipole_indices = self._find_dipole_indices(sequence)

        # Then this chicane isn't in the sequence anyway
        if not dipole_indices:
            return

        if len(dipole_indices) != self.MAX_DIPOLES:
            raise ValueError("Cannot modify partially defined chicane.")

        # Chicane should be consistent before we modify it
        self._assert_consistency_in_chicane(sequence, dipole_indices)

        # The order matters here.  Won't work if drifts are changed first.
        self._update_angles(sequence, dipole_indices)
        self._update_drifts(sequence, dipole_indices)

        # Chicane should be consistent still after we have modified it
        self._assert_consistency_in_chicane(sequence, dipole_indices)

    def __ior__(self, other):
        self.angle = _set_if_not_none(self.angle, other.angle)
        self.r56 = _set_if_not_none(self.r56, other.r56)
        return self


DEFAULT_CONTROLS = {
    "tds1": (CavityController, r"TDSA\.52\.I1"),
    "tds2": (CavityController, r"TDSB\.(428|430)\.B2"),
    "a1": (CavityController, r"C\.A1\.1\.[0-8].I1"),
    "ah1": (CavityController, r"C3\.AH1\.1\.[0-8].I1"),
    "a2": (CavityController, r"C\.A2\.[0-4]\.[0-8].L1"),
    "a3": (CavityController, r"C\.A[3-5]\.[0-4]\.[0-8].L2"),
    "lh": (ChicaneController, "BL\.(48I|48II|50I|50II)\.I1"),
    "bc0": (ChicaneController, r"BB\.(96|98|100|101)\.I1"),
    "bc1": (ChicaneController, r"BB\.(182|191|193|202)\.B1"),
    "bc2": (ChicaneController, r"BB\.(393|402|404|413)\.B2"),
}


def _make_default_controls():
    global DEFAULT_CONTROLS
    return {name: cls(regex) for (name, (cls, regex)) in DEFAULT_CONTROLS.items()}


class FELSimulationConfig:
    def __init__(
        self,
        metadata: Optional[dict] = None,
        controls: Optional[dict[str, Any]] = None,
        components: Optional[dict[str, dict[str, Any]]] = None,
        tracking: Optional[TrackingConfiguration] = None,
    ):
        self.metadata = metadata if metadata else {}
        self.components = components if components else {}
        self.controls = controls if controls else _make_default_controls()
        self.tracking = tracking if tracking else TrackingConfiguration()

    def new_sequence(self, sequence):
        new_sequence = self.apply_controls_to_sequence(sequence)
        new_sequence = self.apply_magnet_conf_to_sequence(new_sequence)
        return new_sequence

    def apply_controls_to_sequence(
        self, sequence: Iterable[ElementT]
    ) -> MachineSequence:
        new_sequence = deepcopy(sequence)

        for element in new_sequence:
            for control in self.controls.values():
                if control.should_modify() and control.match(element):
                    control.modify_element(element)

        for control in self.controls.values():
            control.modify_sequence(new_sequence)
        return new_sequence

    def apply_magnet_conf_to_sequence(
        self, sequence: Iterable[ElementT]
    ) -> MachineSequence:
        sequence = deepcopy(sequence)
        magnet_conf = self.components
        for element in sequence:
            name = element.id
            if name in magnet_conf:
                for key, value in magnet_conf[name].items():
                    previous = getattr(element, key)
                    if np.isclose(previous, value):
                        continue
                    setattr(element, key, value)
                    LOG.debug(f"Set: {name}->{key} to {value}")

        return MachineSequence(sequence)

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""\
        FELSimulationConfig:
        - TDS1: phi = {self.controller.tds1.phi}, v = {self.controller.tds1.v}
        - A1 Cavities: phi = {self.controller.a1.phi}, v = {self.controller.a1.v}
        - AH1 Cavities: phi = {self.controller.ah1.phi}, v = {self.controller.ah1.v}"""
        )

    def __ior__(self, other) -> FELSimulationConfig:
        self.metadata |= other.metadata
        self.components |= other.components
        for name, other_control in other.controls.items():
            try:
                self.controls[name] |= other_control
            except KeyError:
                self.controls[name] = deepcopy(other_control)

        return self

    def __or__(self, other):
        result = deepcopy(self)
        result |= other

        return result


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
        self.twiss = twiss_parray_slice(
            parray,
            slice="Emax",
            nparts_in_slice=5000,
            smooth_param=0.05,
            filter_base=2,
            filter_iter=2,
        )
        self.twiss.s = float(parray.s)


def _twiss_central_slice(parray, match_slice="Imax", **kwargs):
    return twiss_parray_slice(parray, match_slice=match_slice, **kwargs)


def _add_markers_to_navi_for_optics(navi: Navigator, opticscls: Type = None):
    # Put markers everywhere.
    _, new_markers = insert_markers_by_predicate(
        navi.lat.sequence, lambda ele: True, before=False, after_suffix=""
    )

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
        if isinstance(element_names, str):
            fname = Path(outdir) / f"{element_name}.npz"
        elif not isinstance(element_name, str):
            element_name, fname = element_name
            fname = outdir / f"{fname}"
            if fname.suffix != ".npz":
                fname = str(fname) + ".npz"

        proc = SaveBeam(fname)
        proc.element_name = element_name

        try:
            element = next(x for x in navi.lat.sequence if x.id == element_name)
        except StopIteration:
            raise ElementNotFoundError(f"{element_name} not found")

        navi.add_physics_proc(proc, element, element)
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
        twiss.id = proc.name
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


# class PhysicsConfiguration:
#     processes: list


@dataclass
class TrackingConfiguration:
    unit_step: float = 0.02
    method: Transformation = TransferMap


def _set_if_not_none(value, other):
    return value if other is None else other


def _coerce_to_ocelot_twiss(twisslike):
    if isinstance(twisslike, Twiss):
        return twisslike
    return Twiss.from_series(twisslike)
