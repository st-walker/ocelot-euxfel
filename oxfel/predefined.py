import toml
from pathlib import Path
import pickle
from enum import Enum, auto
from typing import Union
from importlib.resources import files

import pandas as pd
from ocelot import Twiss

from oxfel.accelerator.lattice import (i1, i1d, l1, b1d, l2, b2d, l3,
                                       cl, tld, tl1tl2, tl3tl4, tl5,
                                       t1, sa2, t3, t5, t5d)
from oxfel.fel_track import MachineSequence
from oxfel.conversion import TRACKING_CONF_NAME, REAL_CONF_NAME, get_bunchsizerc_path
from oxfel.processes import load_default_physics_lists, PhysicsList
from .xfelt import EuXFELSimConfig, EuXFEL

# I1D_SECTIONS = [G1, A1, AH1, LH, I1D]
# B1D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, B1D]
# B2D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D]

# T4D_SECTIONS = []
# T5D_SECTIONS = []

_TWISS0 = i1.tws0

I1D_SUBSEQUENCES = ["i1", "i1d"]
B1D_SUBSEQUENCES = ["i1", "l1", "b1d"]
B2D_SUBSEQUENCES = ["i1", "l1",  "l2", "b2d"]
TLD_SUBSEQUENCES = ["i1", "l1", "l2", "l3", "cl", "tl1tl2", "tld"]
# T4D_SUBSEQUENCES = ["i1", "l1", "l2", "l3", "cl"]
T5D_SUBSEQUENCES = ["i1", "l1", "l2", "l3", "cl", "tl1tl2", "tl3tl4", "t1", "sa2", "t3", "t5", "t5d"]

def _init_module_level_cells(module_names):
    result = []
    for mname in module_names:
        result.extend(globals()[mname].cell)
    return MachineSequence(result)


class ModelType(str, Enum):
    DESIGN = auto()
    REAL = auto()
    TRACKING = auto()


def _make_real_felconfig() -> EuXFELSimConfig:
    conf = files("oxfel.accelerator.lattice") / REAL_CONF_NAME
    with conf.open("r") as f:
        conf = toml.load(conf)
    return EuXFELSimConfig(components=conf["components"])


def _make_tracking_felconfig() -> EuXFELSimConfig:
    conf = files("oxfel.accelerator.lattice") / TRACKING_CONF_NAME
    with conf.open("r") as f:
        conf = toml.load(conf)
    return EuXFELSimConfig(components=conf["components"])


class ModelBuilder:
    def __init__(self, twiss0: Twiss, subsequence_names: list[str], destination: str):
        self.twiss0 = twiss0
        self.subsequence_names = subsequence_names
        self.destination = destination

    def __call__(
        self, *, model_type: Union[str, ModelType] = ModelType.DESIGN
    ) -> EuXFEL:
        felconfig = self._get_fel_config(model_type)
        physics_list = self._get_physics_list()
        fel = EuXFEL(
            _init_module_level_cells(self.subsequence_names), twiss0=self.twiss0, felconfig=felconfig,
            physics_list=physics_list
        )
        return fel

    def _get_fel_config(self, model_type: Union[str, ModelType]):
        if not isinstance(model_type, ModelType):
            try:
                model_type = ModelType[model_type.upper()]
            except AttributeError:
                raise TypeError(f"Unparsable model_type: {model_type}")

        if model_type is ModelType.DESIGN:
            return None
        elif model_type is ModelType.REAL:
            felconfig = _make_real_felconfig()
        elif model_type is ModelType.TRACKING:
            felconfig = _make_tracking_felconfig()
        else:
            raise TypeError(f"Unknown model_type: {model_type}")

        path = get_bunchsizerc_path(self.destination)
        tracking_optics = pd.read_pickle(path)
        felconfig.physics.beam_sizes = tracking_optics

        return felconfig

    def _get_physics_list(self):
        lists = load_default_physics_lists()
        sequence_lists = [lists[name] for name in self.subsequence_names]
        result = PhysicsList()
        for sublist in sequence_lists:
            result.extend(sublist)
        return result


cat_to_i1d = ModelBuilder(_TWISS0, I1D_SUBSEQUENCES, "i1d")
cat_to_b1d = ModelBuilder(_TWISS0, B1D_SUBSEQUENCES, "b1d")
cat_to_b2d = ModelBuilder(_TWISS0, B2D_SUBSEQUENCES, "b2d")
cat_to_tld = ModelBuilder(_TWISS0, TLD_SUBSEQUENCES, "tld")
cat_to_t5d = ModelBuilder(_TWISS0, T5D_SUBSEQUENCES, "t5d")



def all_models(model_type) -> dict:
    models = {}
    for name, obj in globals().items():
        if name.startswith("cat_to_"):
            models[name] = obj(model_type=model_type)
    return models

