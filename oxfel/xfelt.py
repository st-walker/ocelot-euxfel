from __future__ import annotations

from typing import Optional, Any, Iterable, TypeVar
from copy import deepcopy
import logging

import numpy as np
from ocelot.cpbd.beam import Twiss

from .sequence import ElementT, MachineSequence
from .fel_track import Linac, CavityController, ChicaneController, SimulationConfig, PhysicsConfig, Controller
from .optics import (get_match_point_constraints, get_match_point, INJECTOR_MATCHING_QUAD_NAMES, MATCH_37, MATCH_52, default_match_point_optics)
from .tracking import start_sim_to_match_37
from .matching import match_with_backtracking, MismatchSummary


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.INFO)



DEFAULT_CONTROLS = {
    "tdsi1": (CavityController, r"TDSA\.52\.I1"),
    "tds2": (CavityController, r"TDSB\.(428|430)\.B2"),
    "a1": (CavityController, r"C\.A1\.1\.[0-8].I1"),
    "ah1": (CavityController, r"C3\.AH1\.1\.[0-8].I1"),
    "a2": (CavityController, r"C\.A2\.[0-4]\.[0-8].L1"),
    "a3": (CavityController, r"C\.A[3-5]\.[0-4]\.[0-8].L2"),
    "lh": (ChicaneController, r"BL\.(48I|48II|50I|50II)\.I1"),
    "bc0": (ChicaneController, r"BB\.(96|98|100|101)\.I1"),
    "bc1": (ChicaneController, r"BB\.(182|191|193|202)\.B1"),
    "bc2": (ChicaneController, r"BB\.(393|402|404|413)\.B2"),
}


ControllerT = TypeVar("ControllerT", bound=Controller)

def _make_default_controls() -> dict[str, ControllerT]:
    global DEFAULT_CONTROLS
    return {name: cls(regex) for (name, (cls, regex)) in DEFAULT_CONTROLS.items()}


class EuXFELSimConfig(SimulationConfig):
    def __init__(
        self,
        metadata: dict | None = None,
        controls: dict[str, Any] | None = None,
        components: dict[str, dict[str, Any]] | None = None,
        physics: PhysicsConfig | None = None,
    ):
        self.metadata = metadata if metadata else {}
        self.components = components if components else {}
        self.controls = controls if controls else _make_default_controls()
        self.physics = physics if physics else PhysicsConfig()

    def new_sequence(self, sequence: Iterable[ElementT]) -> MachineSequence:
        new_sequence = self.apply_controls_to_sequence(sequence)
        new_sequence = self.apply_magnet_conf_to_sequence(new_sequence)
        return new_sequence

    def apply_controls_to_sequence(self, sequence: Iterable[ElementT]) -> MachineSequence:
        new_sequence = deepcopy(sequence)

        for element in new_sequence:
            for control in self.controls.values():
                if control.should_modify() and control.match(element):
                    control.modify_element(element)

        for control in self.controls.values():
            control.modify_sequence(new_sequence)
        return new_sequence

    def apply_magnet_conf_to_sequence(self, sequence: Iterable[ElementT]) -> MachineSequence:
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

    def __ior__(self, other: EuXFELSimConfig) -> EuXFELSimConfig:
        self.metadata |= other.metadata
        self.components |= other.components
        for name, other_control in other.controls.items():
            try:
                self.controls[name] |= other_control
            except KeyError:
                self.controls[name] = deepcopy(other_control)

        return self

    def __or__(self, other: EuXFELSimConfig):
        result = deepcopy(self)
        result |= other

        return result

    def to_pickle(self, fpath: os.PathLike) -> None:
        with open(fpath, "wb") as f:
            pickle.dump(self, f)
            LOG.debug(f"Written pickled config to {fpath}")

    def update_components(self, names: Iterable[str], values: Iterable[float | str], attribute: Iterable[float]) -> None:
        self.components.update({name: {attribute: val} for name, val in zip(names, values)})


EuXFELConfig = EuXFELSimConfig



class EuXFEL(Linac):
    SIMCONF_CLASS = EuXFELSimConfig
    def match_injector_linear(self, twiss0: Twiss | None = None, felconfig: EuXFELSimConfig | None=None, start: ElementAccessType | None = None, **match_kwargs) -> EuXFELSimConfig:
        constraints = get_match_point_constraints()
        twiss0 = twiss0 or self.twiss0
        constraints = {"MATCH.52.I1": constraints["MATCH.52.I1"]}
        felconf = self.match_twiss(twiss0,
                                   elements=INJECTOR_MATCHING_QUAD_NAMES,
                                   constraints=constraints,
                                   start=start,
                                   stop="MATCH.52.I1",
                                   felconfig=felconfig,
                                   **match_kwargs)
        return felconf

    def match_injector(self,parray032: ParticleArray | None = None,
                       parray37: ParticleArray | None = None,
                       felconfig: EuXFELSimConfig | None = None,
                       match: str ="projected",
                       maxiter: int = 1)-> tuple[EuXFELSimConfig, MismatchSummary]:
        """Match the EuXFEL Injector for a range of possible starting
        points. `parray37` starts at MATCH.37.I1, `parray032` starts
        at 3.2m.

        """

        if parray032 and parray37:
            raise TypeError("Provide parray032 or parray37, not both.")

        if parray032:
            start = "start_ocelot"
            parray37 = start_sim_to_match_37(self, parray032, felconfig=felconfig)
        elif parray37:
            start = "MATCH.37.I1"

        navi = self.to_navigator(start=MATCH_37, stop=MATCH_52, felconfig=felconfig)
        match52 = default_match_point_optics().set_index("id").loc["MATCH.52.I1"]
        goal_twiss = Twiss(**dict(alpha_x=match52.alpha_x,
                                  alpha_y=match52.alpha_y,
                                  beta_x=match52.beta_x,
                                  beta_y=match52.beta_y,
                                  id="MATCH.52.I1"))

        strengths, mismatch = match_with_backtracking(navi, parray37,
                                                      goal_twiss,
                                                      INJECTOR_MATCHING_QUAD_NAMES,
                                                      maxiter=maxiter,
                                                      match=match)

        if felconfig is None:
            felconfig = EuXFELSimConfig()
        
        felconfig.update_components(INJECTOR_MATCHING_QUAD_NAMES, strengths, "k1")

        return felconfig, mismatch        
