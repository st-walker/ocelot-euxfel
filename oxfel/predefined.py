import toml
from pathlib import Path
import pickle
from enum import Enum, auto
from typing import Union
from importlib_resources import files

from oxfel.accelerator.lattice import i1
from oxfel.accelerator.sections import (
    G1,
    A1,
    AH1,
    LH,
    I1D_Screen,
    DL,
    BC0,
    L1,
    BC1,
    L2,
    BC2,
    B2D,
)
from oxfel.fel_track import SectionedFEL, FELSimulationConfig, FELSection
from oxfel.conversion import TRACKING_CONF_NAME, REAL_CONF_NAME

I1D_SECTIONS = [G1, A1, AH1, LH, I1D_Screen]
B2D_SECTIONS = [G1, A1, AH1, LH, DL, BC0, L1, BC1, L2, BC2, B2D]
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


def _init_i1d_sections() -> list[FELSection]:
    return [x() for x in I1D_SECTIONS]


def _init_b2d_sections() -> list[FELSection]:
    return [x() for x in B2D_SECTIONS]


def cat_to_i1d(*, model_type: Union[str, ModelType] = ModelType.DESIGN) -> SectionedFEL:
    fel = SectionedFEL(_init_i1d_sections(), twiss0=TWISS0)
    _add_config_to_fel(fel, model_type)
    return fel


def cat_to_b2d(*, model_type: Union[str, ModelType] = ModelType.DESIGN) -> SectionedFEL:
    fel = SectionedFEL(_init_b2d_sections(), twiss0=TWISS0)
    _add_config_to_fel(fel, model_type)
    return fel


def _add_config_to_fel(fel, model_type: Union[str, ModelType]):
    if not isinstance(model_type, ModelType):
        try:
            model_type = ModelType[model_type.upper()]
        except AttributeError:
            raise TypeError(f"Unparsable model_type: {model_type}")

    if model_type is ModelType.DESIGN:
        return fel
    if model_type is ModelType.REAL:
        fel.felconfig = _make_real_felconfig()
    elif model_type is ModelType.TRACKING:
        fel.felconfig = _make_tracking_felconfig()
    else:
        raise TypeError(f"Unknown model_type: {model_type}")
    return fel
