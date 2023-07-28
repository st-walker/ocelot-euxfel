import toml
from pathlib import Path
import pickle
from enum import Enum, auto
from typing import Union
from importlib_resources import files
import pandas as pd

from oxfel.accelerator.lattice import i1
from oxfel.accelerator.sections import (
    G1,
    A1,
    AH1,
    LH,
    I1D,
    DL,
    BC0,
    L1,
    BC1,
    B1D,
    L2,
    BC2,
    B2D,
    L3,
    CL1,
    CL2,
    CL3,
    TLD,
)
from oxfel.fel_track import SectionedFEL, FELSimulationConfig, FELSection
from oxfel.conversion import TRACKING_CONF_NAME, REAL_CONF_NAME, get_bunchsizerc_path

I1D_SECTIONS = [G1, A1, AH1, LH, I1D]
B1D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, B1D]
B2D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D]
TLD_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, L3, CL1, CL2, CL3, TLD]
T4D_SECTIONS = []
T5D_SECTIONS = []

TWISS0 = i1.tws0


class ModelType(str, Enum):
    DESIGN = auto()
    REAL = auto()
    TRACKING = auto()


def _make_real_felconfig() -> FELSimulationConfig:
    conf = files("oxfel.accelerator.lattice") / REAL_CONF_NAME
    with conf.open("r") as f:
        conf = toml.load(conf)
    return FELSimulationConfig(components=conf["components"])


def _make_tracking_felconfig() -> FELSimulationConfig:
    conf = files("oxfel.accelerator.lattice") / TRACKING_CONF_NAME
    with conf.open("r") as f:
        conf = toml.load(conf)
    return FELSimulationConfig(components=conf["components"])


def _init_sections(sections) -> list[FELSection]:
    return [x() for x in sections]


class ModelBuilder:
    def __init__(self, twiss0, sections, destination):
        self.twiss0 = twiss0
        self.sections = sections
        self.destination = destination

    def _init_sections(self):
        return [x() for x in self.sections]

    def __call__(
        self, *, model_type: Union[str, ModelType] = ModelType.DESIGN
    ) -> SectionedFEL:
        felconfig = self._get_fel_config(model_type)
        fel = SectionedFEL(
            self._init_sections(), twiss0=self.twiss0, felconfig=felconfig
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


cat_to_i1d = ModelBuilder(TWISS0, I1D_SECTIONS, "i1d")
cat_to_b1d = ModelBuilder(TWISS0, B1D_SECTIONS, "b1d")
cat_to_b2d = ModelBuilder(TWISS0, B2D_SECTIONS, "b2d")
cat_to_tld = ModelBuilder(TWISS0, TLD_SECTIONS, "tld")
cat_to_t4d = ModelBuilder(TWISS0, T4D_SECTIONS, "t4d")


def all_models(model_type) -> dict:
    models = {}
    for name, obj in globals().items():
        if name.startswith("cat_to_"):
            models[name] = obj(model_type=model_type)
    return models
