from typing import Optional
import os

import numpy as np
import numpy.typing as npt
from ocelot.cpbd.beam import Twiss

import pandas as pd

from .longlist import get_default_component_list, XFELComponentList

INJECTOR_MATCHING_QUAD_NAMES: list[str] = [
    "Q.37.I1",
    "Q.38.I1",
    "QI.46.I1",
    "QI.47.I1",
    "QI.50.I1",
]
B2_MATCHING_QUAD_NAMES: list[str] = [
    "Q.333.L2",
    "Q.345.L2",
    "Q.357.L2",
    "Q.369.L2",
    "Q.381.L2",
]

#
START_SIM: str = "start_ocelot"
MATCH_37: str = "MATCH.37.I1"
MATCH_52: str = "MATCH.52.I1"


# Not a tru "fixed" point, just in front of the TDS2
MATCH_428: str = "MATCH.428.B2"


FIXED_MATCH_POINTS: list[str] = [
    # In front of TDS must be fixed
    "MATCH.52.I1",
    # In front of dogleg must be fixed
    "MATCH.73.I1",
    # ??? Why must this be fixed ???
    "MATCH.104.I1",
    # In front of BC1 TDS must be fixed
    "MATCH.218.B1",
    "MATCH.446.B2",
    # Entrance to L3 must be fixed
    "MATCH.525.L3",
    # Entrance to collimation dogleg must be fixed
    "MATCH.1673.CL"
]

ALL_INTERESTING_MATCH_POINTS = FIXED_MATCH_POINTS + [MATCH_428]

FIXED_MATCH_POINTS = ALL_INTERESTING_MATCH_POINTS

_OCELOT_OPTICS_NAMES = ["beta_x", "alpha_x", "beta_y", "alpha_y"]
_OCELOT_OTHER_NAMES = ["id", "s"]
_OCELOT_TWISS_NAMES = _OCELOT_OTHER_NAMES + _OCELOT_OPTICS_NAMES



def get_name2_fixed_match_points(additional_names: Optional[list[str]] = None):
    if additional_names is None:
        additional_names = []
    ll = get_default_component_list()
    return [
        ll.name1_to_name2(name1).item()
        for name1 in FIXED_MATCH_POINTS + additional_names
    ]


def get_match_point_constraints(longlist: Optional[XFELComponentList] = None) -> dict[str, dict[str, float]]:
    if longlist is None:
        longlist = get_default_component_list()
    points = ALL_INTERESTING_MATCH_POINTS

    tups = [x for x in longlist.df.itertuples() if x.NAME1 in points]

    constraints = {}
    for tup in tups:
        odict = {"alpha_x": tup.ALFX,
                 "beta_x": tup.BETX,
                 "alpha_y": tup.ALFY,
                 "beta_y": tup.BETY}
        
        constraints[tup.NAME1] = odict

    return constraints


def _normalise_twiss_df(df: pd.DataFrame) -> pd.DataFrame:
    # Change columns to ocelot column names
    # Change NAME2 to NAME1s.
    df = df.rename(
        mapper={
            "BETX": "beta_x",
            "SUML": "s",
            "BETY": "beta_y",
            "ALFX": "alpha_x",
            "ALFY": "alpha_y",
            "ENERGY": "E",
            "S": "s",
            "NAME1": "id",
            "NAME": "id",
        },
        axis=1,
    )
    return df


def get_match_point_optics(twiss_df: pd.DataFrame, additional_names: Optional[list[str]] = None) -> pd.DataFrame:
    twiss_df = _normalise_twiss_df(twiss_df)

    if additional_names is None:
        additional_names = []

    point_names = FIXED_MATCH_POINTS + additional_names
    names = [name for name in point_names if name in np.array(twiss_df["id"])]

    try:
        return twiss_df.set_index("id").loc[names].reset_index()[_OCELOT_TWISS_NAMES]
    except TypeError:
        raise ValueError("Unable to extract, possibly duplicate names in twiss_df id column")
    
    

def default_match_point_optics():
    df = get_default_component_list().longlist
    df = _normalise_twiss_df(df)

    match_points = df.set_index("id").loc[FIXED_MATCH_POINTS].reset_index()
    twiss = match_points[_OCELOT_TWISS_NAMES]
    return twiss


def get_default_match_point(match_name: str) -> Twiss:
    df = default_match_point_optics()
    twiss_series = df[df.id == match_name]
    
    return Twiss.from_series(twiss_series)

def print_match_point_analysis(twiss_or_twiss_df: pd.DataFrame,
                               additional_names: Optional[list[str]] = None) -> None:
    try:
        print(get_match_point(twiss_or_twiss_df, additional_names=additional_names))
        return
    except AttributeError:
        pass
    except:
        raise TypeError(f"Unknown Twiss type: {twiss_or_twiss_df}")

    print(get_match_point(twiss_or_twiss_df.to_series().to_frame().T, additional_names=additional_names))

def get_match_point(twiss_df: pd.DataFrame, additional_names: Optional[list[str]] = None) -> pd.DataFrame:
    twiss_match = get_match_point_optics(twiss_df, additional_names=additional_names)
    twiss_match_reference = default_match_point_optics()

    # Calculate relative difference
    comparison = twiss_match.copy()[_OCELOT_TWISS_NAMES]

    # use id as index here so that the correct rows are used with each other.
    twiss_match = twiss_match.set_index("id")
    twiss_match_reference = twiss_match_reference.set_index("id")

    gamma_x = (1 + twiss_match.alpha_x * twiss_match.alpha_x) / twiss_match.beta_x
    gamma_y = (1 + twiss_match.alpha_y * twiss_match.alpha_y) / twiss_match.beta_y

    twiss_match["gamma_x"] = gamma_x
    twiss_match["gamma_y"] = gamma_y

    bmag_x = bmag(
        twiss_match.beta_x,
        twiss_match.alpha_x,
        twiss_match_reference.beta_x,
        twiss_match_reference.alpha_x,
    )
    bmag_y = bmag(
        twiss_match.beta_y,
        twiss_match.alpha_y,
        twiss_match_reference.beta_y,
        twiss_match_reference.alpha_y,
    )

    twiss_match["bmag_x"] = bmag_x
    twiss_match["bmag_y"] = bmag_y

    return twiss_match.sort_values("s")


def read_mad8(fname: os.PathLike) -> pd.DataFrame:
    import pand8

    df8 = pand8.read(fname)
    df8 = df8.rename(mapper={"NAME": "NAME2"}, axis=1)

    ll = get_default_component_list()

    di = dict(zip(ll.df.NAME2, ll.df.NAME1))

    name1s = [di.get(name2, name2) for name2 in df8.NAME2]

    df8["NAME1"] = name1s

    return _normalise_twiss_df(df8)


def bmag(beta: npt.ArrayLike, alpha: npt.ArrayLike,
         beta_design: npt.ArrayLike, alpha_design: npt.ArrayLike) -> np.ndarray:
    bmag = 0.5 * (
        (beta / beta_design + beta_design / beta)
        + (beta * beta_design * ((alpha_design / beta_design) - (alpha / beta)) ** 2)
    )
    return np.squeeze(bmag)
