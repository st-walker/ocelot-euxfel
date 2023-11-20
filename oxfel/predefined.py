import toml
from pathlib import Path
import pickle
from enum import Enum, auto
from typing import Union
from importlib_resources import files
import pandas as pd

from oxfel.accelerator.lattice import i1, i1d
from oxfel.fel_track import Linac, EuXFELSimConfig, MachineSequence
from oxfel.conversion import TRACKING_CONF_NAME, REAL_CONF_NAME, get_bunchsizerc_path
from oxfel.processes import load_default_physics_lists, PhysicsList

# I1D_SECTIONS = [G1, A1, AH1, LH, I1D]
# B1D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, B1D]
# B2D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D]
# TLD_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, L3, CL1, CL2, CL3, TLD]
# T4D_SECTIONS = []
# T5D_SECTIONS = []

TWISS0 = i1.tws0

I1D_SUBSEQUENCES = ["i1", "i1d"]
B2D_SUBSEQUENCES = []


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
    def __init__(self, twiss0, subsequence_names, destination):
        self.twiss0 = twiss0
        self.subsequence_names = subsequence_names
        self.destination = destination

    def __call__(
        self, *, model_type: Union[str, ModelType] = ModelType.DESIGN
    ) -> Linac:
        felconfig = self._get_fel_config(model_type)
        physics_list = self._get_physics_list()
        fel = Linac(
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


cat_to_i1d = ModelBuilder(TWISS0, I1D_SUBSEQUENCES, "i1d")
# cat_to_b1d = ModelBuilder(TWISS0, B1D_SECTIONS, "b1d")
cat_to_b2d = ModelBuilder(TWISS0, B2D_SUBSEQUENCES, "b2d")
# cat_to_tld = ModelBuilder(TWISS0, TLD_SECTIONS, "tld")
# cat_to_t4d = ModelBuilder(TWISS0, T4D_SECTIONS, "t4d")



def all_models(model_type) -> dict:
    models = {}
    for name, obj in globals().items():
        if name.startswith("cat_to_"):
            models[name] = obj(model_type=model_type)
    return models

