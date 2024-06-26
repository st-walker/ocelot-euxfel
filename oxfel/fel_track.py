from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Type, Iterable, TypeVar, Callable, Any
from pathlib import Path
import logging
from copy import deepcopy
import textwrap
import re
import pandas as pd
from dataclasses import dataclass, field
import functools
import numbers
import pickle
from collections import defaultdict
from .processes import PhysicsList
from copy import deepcopy


import numpy as np
from ocelot.cpbd.csr import CSR
from ocelot.cpbd.beam import ParticleArray, global_slice_analysis
from ocelot.cpbd.physics_proc import PhysProc, CopyBeam
import ocelot.cpbd.track as track
from ocelot.cpbd.magnetic_lattice import (
    MagneticLattice,
    insert_markers_by_predicate,
    flatten,
)
from ocelot.cpbd.beam import twiss_parray_slice, get_envelope
from ocelot.cpbd.transformations import TransferMap
from ocelot.cpbd.transformations.transformation import Transformation
from ocelot.cpbd.beam import twiss_iterable_to_df
from ocelot.cpbd.optics import Twiss, twiss as oce_calc_twiss
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.beam import optics_from_moments, moments_from_parray
from ocelot.cpbd.match import match, match_beam


from .sequence import MachineSequence, ElementT, ElementAccessType, ElementAccessError


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.INFO)


TrackingStartPoint = ElementAccessType | float
ElementSequenceT = list[ElementT]


@dataclass
class TrackingResult:
    parray1: ParticleArray
    sequence: MachineSequence
    twiss: pd.DataFrame | None = None
    dumps: dict[str, ParticleArray] | None = None


class SimulationConfig:
    pass


class Linac:
    SIMCONF_CLASS = SimulationConfig
    def __init__(
        self,
        sequence: MachineSequence | list,
        twiss0: Twiss,
        outdir: os.PathLike = "./",
        felconfig=None,
        physics_list=None
    ):
        self.sequence = MachineSequence(sequence)
        self.twiss0 = twiss0
        self.felconfig = felconfig if felconfig else None
        self.physics_list = physics_list if physics_list else PhysicsList([])

    def _net_felconfig(self, felconfig2=None) -> SimulationConfig:
        try:
            conf_path = Path(felconfig2)
        except TypeError:
            pass
        else:
            with conf_path.open("rb") as f:
                felconfig2 = pickle.load(f)

        felconfig = self.felconfig
        if felconfig is None:
            felconfig = self.SIMCONF_CLASS()

        if felconfig2 is None:
            felconfig2 = self.SIMCONF_CLASS()

        return felconfig | felconfig2

    # def get_start_element_name(self, start: Optional[str] = None, twiss0: Optional[Twiss] = None) -> str:
    #     # If both have been provided then favour and explicit start parameter
    #     try:
    #         twiss0id = twiss0.id
    #     except AttributeError:
    #         twiss0id = ""
            
    #     if start:
    #         idstart = start
    #     elif twiss0id:
    #         idstart = twiss0id
    #     elif self.twiss0.id:
    #         idstart = self.twiss0.id
    #     else:
    #         idstart = self.sequence[0].id
    #     LOG.debug(f"Returning {idstart=} from {start=}, {twiss0=} and {self.twiss0.id=}")
    #     return idstart

    def _calculate_twiss_between_two_points(
            self, twiss0: Twiss, start: ElementAccessType | None = None,
            stop: ElementAccessType | None = None,
            felconfig: SimulationConfig | None = None
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        # self._check_twiss_position_matches(twiss0, start)
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        mlat = MagneticLattice(sequence)

        # mangle initial name a bit for clarity.
        twiss0 = deepcopy(twiss0)
        twiss0.id = "$fel_twiss_start"

        full_twiss = oce_calc_twiss(mlat, tws0=twiss0, return_df=True)
        return full_twiss, mlat

    def design_twiss(self, stop: ElementAccessType = None) -> tuple[Twiss, MagneticLattice]:
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=self.twiss0, stop=stop
        )
        return twiss, mlat

    def machine_twiss(
        self,
        stop: ElementAccessType | None = None,
        felconfig: SimulationConfig | None = None,
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        """Calculates the Twiss parameters and magnetic lattice
        configuration to the stop point, by default the end, using the
        member twiss0 instance as initial Twiss parameters.
        
        Parameters:
        - stop (ElementAccessType, optional): The stopping point in the machine for calculation.
        - felconfig (SimulationConfig, optional): Configuration for Free Electron Laser simulation.
        
        Returns:
        - tuple[pd.DataFrame, MagneticLattice]: A tuple containing a DataFrame of calculated Twiss parameters
        and the magnetic lattice configuration.

        """
        felconfig = self._net_felconfig(felconfig)
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=self.twiss0, stop=stop, felconfig=felconfig
        )
        return twiss, mlat

    def calculate_twiss(
        self, twiss0: Twiss, start: ElementAccessType, stop: ElementAccessType | None = None, felconfig: SimulationConfig | None = None
    ) -> tuple[pd.DataFrame, MagneticLattice]:
        """
        Calculates the Twiss parameters for the Linac instance over a specified section.

        This method computes the linear optics (Twiss parameters) between two points specified by `start` and `stop`
        within the magnetic lattice of a linear accelerator (Linac). It leverages an initial set of Twiss parameters
        (`twiss0`) and optionally considers an EuXFELSimXonfig during the computation.

        Parameters:
        - twiss0: A structure representing the initial Twiss parameters from which the calculation starts. This structure
                  must be compatible with the expected format for Twiss parameters in the Ocelot framework.
        - start (ElementAccessType): The name or identifier of the starting element in the magnetic lattice from which the calculation
                       begins.
        - stop (Optional[str]): The name or identifier of the stopping element in the magnetic lattice at which the
                                calculation ends. If `None`, the calculation proceeds to the end of the lattice or until
                                another stopping condition is met.
        - felconfig (Optional[SimulationConfig]): An optional configuration for the Linac section.

        Returns:
        A tuple containing two elements:
        - A pandas DataFrame containing the calculated Twiss parameters across the specified section of the Linac. Each row
          corresponds to a position within the magnetic lattice, and columns include relevant Twiss parameters and possibly
          additional information depending on the calculation specifics.
        - An instance of `MagneticLattice` representing the segment of the Linac over which the Twiss parameters were
          calculated. This object provides access to the magnetic lattice structure, including elements and their
          properties, over the specified range.  Useful for plotting.

        """
        # start = self.get_start_element_name(twiss0=twiss0, start=start)
        felconfig_net = self._net_felconfig(felconfig)
        LOG.debug(f'Calculating linear optics between "{start}" and "{stop}"')
        # Final stretch, between the matching point and the provided "stop"
        twiss0 = _coerce_to_ocelot_twiss(twiss0)
        twiss, mlat = self._calculate_twiss_between_two_points(
            twiss0=twiss0, start=start, stop=stop, felconfig=felconfig_net
        )

        return twiss, mlat

    def track(
        self,
            parray0: ParticleArray,
            start: TrackingStartPoint | None = None,
            stop: ElementAccessType = None,
            felconfig: SimulationConfig | None = None,
        dumps: Iterable[str] | None = None,
        physics: bool = True,
    ) -> TrackingResult:
        return self._do_tracking(
            parray0,
            calculate_optics=False,
            start=start,
            stop=stop,
            felconfig=felconfig,
            dumps=dumps,
            physics=physics,
        )

    def track_optics(
        self,
        parray0: ParticleArray,
        start: TrackingStartPoint | None = None,
        stop: ElementAccessType | None = None,
            felconfig: SimulationConfig | None = None,
        dumps: Iterable[str] | None = None,
        physics: bool = False,
        opticscls: type | str | None = None,
        outdir: os.PathLike = None,
    ) -> TrackingResult:
        tracking_result = self._do_tracking(
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
        return tracking_result

    def _do_tracking(
        self,
        parray0: ParticleArray,
        *,
        calculate_optics: bool,
        start: TrackingStartPoint = None,
        stop: ElementAccessType = None,
        dumps: Iterable[str] | None = None,
        physics: bool = False,
        opticscls: type | str | None = None,
        felconfig: SimulationConfig | None = None,
        **track_kwargs,
    ) -> TrackingResult:
        if isinstance(start, numbers.Number):
            zstart = start
            start = None
        else:
            zstart = 0.0

        navi = self.to_navigator(
            start=start, stop=stop, felconfig=felconfig, physics=physics
        )
        
        if calculate_optics:
            opticscls = _add_markers_to_navi_for_optics(navi, opticscls)

        if dumps is None:
            dumps = []

        copy_procs = _add_copy_beams(navi, dumps)
        # It is important to do this AFTER manipulating the navigator
        # (with copybeam and twiss calc instances being added) because
        # otherwise we end up with everything being in a strange and
        # bad state.  in principle could resolve it here but i don't know now.
        navi.jump_to(zstart)
        
        LOG.debug(f'Starting tracking between "{start}" and "{stop}"')
        _, parray1 = track.track(
            navi.lat,
            parray0.copy(),
            navi=navi,
            calc_tws=False,
            print_progress=True,
            overwrite_progress=False,
        )

        twiss = pd.DataFrame()
        if calculate_optics:
            twiss = _extract_twiss_parameters_from_twiss_processes(navi.inactive_processes)

        
        parrays_at_dumps = {}
        if dumps:
            for process in navi.inactive_processes:
                try:
                    parrays_at_dumps[process.name] = process.copied_parray
                except AttributeError:
                    pass

        return TrackingResult(parray1, navi.lat.sequence, twiss, parrays_at_dumps)

    def to_navigator(
        self,
        start: ElementAccessType | None = None,
        stop: ElementAccessType | None = None,
        felconfig: SimulationConfig | None = None,
        physics: bool = True,
    ) -> Navigator:
        felconfig = self._net_felconfig(felconfig)

        # Make the navigator
        sequence = self.get_sequence(start=start, stop=stop, felconfig=felconfig)
        mlat = MagneticLattice(sequence, method={"global": felconfig.physics.method})
        navi = Navigator(mlat, unit_step=felconfig.physics.unit_step)

        if not physics:
            return navi

        names_before = set(self._get_sequence_before_start(start).names())
        names_sequence = set(sequence.names())
        names_after = set(self._get_sequence_after_stop(stop).names())
        all_names = names_before | names_after | names_sequence

        for placed_process in self.physics_list:
            start = placed_process.start_name
            stop = placed_process.stop_name
            process = placed_process.process

            if start not in all_names:
                raise ElementAccessError(f"Element named {start} not found")
            if stop not in all_names:
                # from IPython import embed; embed()
                raise ElementAccessError(f"Element named {stop} not found")

            # If process stops before sequence starts, then skip.
            if stop in names_before:
                LOG.debug(f"Skipping process that stops before sequence starts: {placed_process}")
                continue

            # If process stops after sequence ends, then skip.
            if start in names_after:
                LOG.debug(f"Skipping process that starts after sequence stops: {placed_process}")
                continue

            # If process starts before sequence starts, then set start to sequence start.
            if start in names_before:
                start = sequence[0].id

            # If procss stops after sequence stops, then set stop to sequence stop.
            if stop in names_after:
                stop = sequence[-1].id

            # Set energy for CSR processes incase it needs it.
            if isinstance(process, CSR):
                self._set_csr_process_energy(process, start, stop, felconfig=felconfig)
            LOG.info(f"Adding {process}: start={start}, stop={stop}:")

            LOG.info(repr(process))

            istart = sequence.names().index(start)
            istop = sequence.names().index(stop)

            element_start = sequence[istart]
            element_stop = sequence[istop]

            navi.add_physics_proc(process, element_start, element_stop)

        return navi

    def _get_sequence_before_start(self, start: ElementAccessType) -> MachineSequence:
        if start is None:
            return MachineSequence([])
        else:
            return self.get_sequence(stop=start)[:-1]

    def _get_sequence_after_stop(self, stop: ElementAccessType) -> MachineSequence:
        if stop is None:
            return MachineSequence([])
        else:
            return self.get_sequence(start=stop)[1:]

    def _set_csr_process_energy(
        self,
        process: PhysProc,
        start_name: ElementAccessType,
        stop_name: ElementAccessType,
        felconfig: SimulationConfig,
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
        felconfig: SimulationConfig = None,
    ) -> MachineSequence:
        LOG.debug(f"Building sequence, {start=}, {stop=}, {felconfig=}")
        result = self.sequence.closed_interval(start, stop)
        net_felconfig = self._net_felconfig(felconfig2=felconfig)
        return MachineSequence(net_felconfig.new_sequence(result))

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

    def simple_match(
        self,
        parray0: ParticleArray,
        *,
        start: ElementAccessType,
        stop: ElementAccessType,
        quad_names: list[str],
        twiss_goal: Twiss,
        maxiter=1,
        felconfig: SimulationConfig | None = None,
        match: Callable[[ParticleArray], Twiss] | str = "projected",
    ):
        felconfig = self._net_felconfig(felconfig)

        navi = self.to_navigator(start=start, stop=stop, felconfig=felconfig)
        strengths, _ = match_with_backtracking(
            navi, parray0, twiss_goal, quad_names, maxiter=maxiter
        )

        cdict = {name: {"k1": k1} for name, k1 in zip(quad_names, strengths)}
        felconfig.components.update(cdict)
        return felconfig

    def match_beam(
        self,
        parray0: ParticleArray,
        *,
        elements: list[str],
        constraints: dict[str : dict[str, float]],
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: SimulationConfig | None = None,
        physics=True,
        match="projected",
        **match_beam_kwargs,
    ) -> SimulationConfig:
        # Make navigator from given start, stop and felconfig
        navi = self.to_navigator(
            start=start, stop=stop, felconfig=felconfig, physics=physics
        )

        if match == "projected":
            twiss_function = None
        elif match.lower() == "emax":
            twiss_function = functools.partial(_twiss_central_slice, match_slice="Emax")
        elif match.lower() == "imax":
            twiss_function = functools.partial(_twiss_central_slice, match_slice="Imax")
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

    def match_twiss(
        self,
        twiss0: Twiss,
        *,
        elements: list[str],
        constraints: dict[str : dict[str, float]],
        start: ElementAccessType = None,
        stop: ElementAccessType = None,
        felconfig: SimulationConfig | None = None,
        verbose: bool = True,
        **match_kwargs,
    ) -> SimulationConfig:
        """
        Match the optics for this Linac instance starting with the
            input optics twiss0 at position `start` and ending at `stop`.
            The elements kwarg should be a list of names of elements that
            will be varied (e.g. most likely quads).  

        Parameters:
        - twiss0 (Twiss): The initial, input twiss to be propagated through the lattice.
        - elements (list[str]): The names of the Quadrupoles whose strengths are to be varied.
        - constraints (dict[str: dict[str, float]]): The constraints, of the same form as required by `ocelot.match.match`.
        - start (ElementAccessType, optional): Name/index of first element to be tracked from.  If None, start from beginning of the sequence.
        - stop (ElementAccessType, optional): Name/index of last element to be tracked to.  If None, stop at end of the sequence.
        - felconfig (Optional[SimulationConfig], optional): Optional FELConfig instance to be applied to the sequence before matching.
        - verbose (bool, optional): Description of verbose. If True, provides detailed logging. Default is True.
        - match_kwargs (dict): Additional keyword arguments for matching.

        Returns:
        - SimulationConfig: The SimulationConfig with matched Quadrupoles.

        """

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

        strengths = match(
            lat_matching, constraints, elements, twiss0, verbose=verbose, **match_kwargs
        )


        # Apply the result to the simulation config.
        new_conf = deepcopy(self._net_felconfig(felconfig))

        matched_strengths = {
            quad_name: {"k1": strength}
            for (quad_name, strength) in zip([x.id for x in elements], strengths)
        }

        new_conf.components.update(matched_strengths)

        return new_conf

    # def _new_match_beam_linear_match(self,
    #                                  parray0: ParticleArray,
    #                                  *,
    #                                  quad_names: list[str],
    #                                  twiss_goal: Twiss,
    #                                  # constraints: dict[str : dict[str, float]],
    #                                  start: ElementAccessType = None,w
    #                                  stop: ElementAccessType = None,
    #                                  felconfig: Optional[FELSimulationConfig] = None,
    #                                  physics=True,
    #                                  match="projected"):
    #     sequence = self.get_sequence(start=start, stop=stop)

    #     twiss_function = _string_to_twiss_match_function(match)
    #     twiss0 = twiss_function(parray0)

    #     linear_matched_conf = fel.match(
    #         twiss0,
    #         # felconfig=felconfig,
    #         start=start,
    #         stop=stop,
    #         constraints=constraints,
    #         elements=,
    #         verbose=True,
    #         max_iter=100000000
    #     )


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
        phi: float | None = None,
        v: float | None = None,
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
                f"Outer chicane drifts have different lengths: {first_drift_length=}, {third_drift_length}="
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
        # self.twiss = get_envelope(parray)
        # self.twiss.s = float(parray.s)
        mean, cov = moments_from_parray(parray)
        self.twiss = optics_from_moments(mean, cov, parray.E)
        self.twiss.s = float(parray.s)
        self.twiss.tautau = parray.tau().std() ** 2

    def __str__(self):
        return f"<{type(self).__name__}, name={self.name}>"


class SliceTwissCalculator(PhysProc):
    """Calculate dispersion-free Twiss parameters"""

    def __init__(self, name: str):
        super().__init__()
        self.name = name
        self.twiss = None

    def apply(self, parray: ParticleArray, dz: float):
        slice_params = global_slice_analysis(
            parray,
            nparts_in_slice=5000,
            smooth_param=0.05,
            filter_base=2,
            filter_iter=2,
        )

        if self.MATCH.lower() == "imax":
            ind0 = np.argmax(slice_params.I)
        elif self.MATCH.lower() == "emax":
            ind0 = np.argmax(slice_params.me)

        self.twiss = slice_params.extract_slice(ind0)
        self.twiss.s = float(parray.s)
        self.twiss.E = float(parray.E)



class PeakEnergyTwissCalculator(SliceTwissCalculator):
    MATCH = "Emax"


class PeakCurrentTwissCalculator(SliceTwissCalculator):
    MATCH = "Imax"


class FullTwissCalculator(PhysProc):
    def __init__(self, name: str):
        super().__init__()
        self.name = name
        self.twiss = None

    def apply(self, parray: ParticleArray, dz: float) -> None:
        slice_params = global_slice_analysis(
            parray,
            nparts_in_slice=5000,
            smooth_param=0.05,
            filter_base=2,
            filter_iter=2,
        )

        if self.MATCH.lower() == "imax":
            ind0 = np.argmax(slice_params.I)
        elif self.MATCH.lower() == "emax":
            ind0 = np.argmax(slice_params.me)

        self.twiss = slice_params.extract_slice(ind0)
        self.twiss.s = float(parray.s)
        self.twiss.E = float(parray.E)


def _twiss_central_slice(parray: ParticleArray, match_slice: str = "Imax", **kwargs) -> Twiss:
    return twiss_parray_slice(parray, slice=match_slice, **kwargs)


def _add_markers_to_navi_for_optics(navi: Navigator, opticscls: type | str | None = None) -> type:
    """Add markers to the navigator with attached physics processes
    for the type of optics we are interested in calculating

    """

    if not opticscls:
        opticscls = TwissCalculator
        
    try:
        opticsname = opticscls.lower()
    except AttributeError:
        pass
    else:
        if opticsname == "projected":
            opticscls = TwissCalculator            
        elif opticsname == "emax":
            opticscls = PeakEnergyTwissCalculator
        elif opticsname == "imax":
            opticscls = PeakCurrentTwissCalculator
        else:
            raise ValueError(f"Unknown named class provided: {opticscls}")


    # Put markers everywhere.
    new_markers = insert_markers_by_predicate(
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

    return opticscls


def _add_copy_beams(navi: Navigator, element_names: list[str]) -> list[CopyBeam]:
    rprocs = []
    for element_name in element_names:
        # Make CopyBeam instance for this element name
        proc = CopyBeam(element_name)
        # find the element if it exists else bomb
        try:
            element = next(x for x in navi.lat.sequence if x.id == element_name)
        except StopIteration:
            raise ElementNotFoundError(f"{element_name} not found")

        # Add the CopyBeam instance to the navigator
        navi.add_physics_proc(proc, element, element)
        LOG.debug(f"Added CopyBeam process to {element_name}")
        rprocs.append(proc)
    return rprocs


def _extract_twiss_parameters_from_twiss_processes(processes: list[PhysProc]) -> pd.DataFrame:
    # We filter Nones just because at the very end of the lattice
    # some markers are not tracked through due the tolerance check
    # in track.track (markers have 0. length).  Sort also by s.
    twiss_instances = []
    for proc in processes:
        try:
            twiss = proc.twiss
        except AttributeError:
            continue

        if twiss is None:
            continue

        twiss.id = proc.name
        twiss_instances.append(twiss)
    twiss_instances.sort(key=lambda tws: tws.s)

    twissdf = twiss_iterable_to_df(twiss_instances)

    return twissdf


def _get_element_instances_from_mlat(mlat: MagneticLattice, names: list[str]) -> list[ElementT]:
    def f(ele):
        return ele.id in names

    indices = mlat.find_indices_by_predicate(f)
    return [mlat.sequence[i] for i in indices]


@dataclass
class PhysicsConfig:
    unit_step: float = 0.02
    method: Transformation = TransferMap
    beam_sizes: pd.DataFrame = None
    csr_sigma_factor: float = 0.1
    match_points: list[str] = field(default_factory=list)

    def get_csr_sigma_min(self, start, stop):
        subseq = self.beam_sizes.set_index("id").loc[start:stop]
        rms_bunch_length = np.sqrt(subseq.tautau)
        return self.csr_sigma_factor * min(rms_bunch_length)


def _set_if_not_none(value, other):
    return value if other is None else other


def _coerce_to_ocelot_twiss(twisslike):
    if isinstance(twisslike, Twiss):
        return twisslike
    return Twiss.from_series(twisslike)


def _string_to_twiss_match_function(match_string: str):
    if match_string == "projected":
        twiss_function = get_envelope
    elif match_string.lower() == "emax":
        twiss_function = functools.partial(_twiss_central_slice, match_slice="Emax")
    elif match_string.lower() == "imax":
        twiss_function = functools.partial(_twiss_central_slice, match_slice="Imax")
    else:
        raise ValueError(f"Unknown match string: {match_string}")
    return twiss_function
