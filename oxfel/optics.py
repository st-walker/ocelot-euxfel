from .longlist import make_default_longlist
import numpy as np


FIXED_MATCH_POINTS = [
    "MATCH.52.I1",
    "MATCH.73.I1",
    "MATCH.104.I1",
    "MATCH.218.B1",
    "MATCH.446.B2",
]


_OCELOT_OPTICS_NAMES = ["beta_x", "alpha_x", "beta_y", "alpha_y"]
_OCELOT_OTHER_NAMES = ["id", "s"]
_OCELOT_TWISS_NAMES = _OCELOT_OTHER_NAMES + _OCELOT_OPTICS_NAMES


def get_name2_fixed_match_points(additional_names=None):
    if additional_names is None:
        additional_names = []
    ll = make_default_longlist()
    return [
        ll.name1_to_name2(name1).item()
        for name1 in FIXED_MATCH_POINTS + additional_names
    ]


def _normalise_twiss_df(df):
    # Change columns to ocelot column names
    # Change NAME2 to NAME1s.
    df = df.rename(mapper={"BETX": "beta_x",
                           "SUML": "s",
                           "BETY": "beta_y",
                           "ALFX": "alpha_x",
                           "ALFY": "alpha_y",
                           "S": "s",
                           "NAME1": "id",
                           "NAME": "id"}, axis=1)
    return df

def get_match_point_optics(twiss_df, additional_names=None):
    twiss_df = _normalise_twiss_df(twiss_df)

    if additional_names is None:
        additional_names = []

    point_names = FIXED_MATCH_POINTS + additional_names
    names = [name for name in point_names if name in np.array(twiss_df["id"])]
    return twiss_df.set_index("id").loc[names].reset_index()[_OCELOT_TWISS_NAMES]


def default_match_point_optics():
    df = make_default_longlist().df
    df = _normalise_twiss_df(df)

    match_points = df.set_index("id").loc[FIXED_MATCH_POINTS].reset_index()
    twiss =  match_points[_OCELOT_TWISS_NAMES]

    return twiss

def print_match_point_analysis(twiss_df, additional_names=None):
    twiss_match = get_match_point_optics(twiss_df, additional_names=additional_names)
    twiss_match_reference = default_match_point_optics()

    # Calculate relative difference
    comparison = twiss_match.copy()[_OCELOT_TWISS_NAMES]
    # comparison = comparison[_OCELOT_OPTICS_NAMES] - twiss_match_reference[_OCELOT_OPTICS_NAMES]
    comparison[_OCELOT_OPTICS_NAMES] = 100 * np.abs((twiss_match[_OCELOT_OPTICS_NAMES] - twiss_match_reference[_OCELOT_OPTICS_NAMES])
                                        / twiss_match_reference[_OCELOT_OPTICS_NAMES])
    comparison = comparison.dropna()

    print("\nRelative errors (%) with respect to the design MATCH settings:")
    print(comparison)
    print("\n")


def read_mad8(fname):
    import pand8
    df8 = pand8.read(fname)
    df8 = df8.rename(mapper={"NAME": "NAME2"}, axis=1)

    ll = make_default_longlist()

    di = dict(zip(ll.df.NAME2, ll.df.NAME1))

    name1s = [di.get(name2, name2) for name2 in df8.NAME2]

    df8["NAME1"] = name1s

    return _normalise_twiss_df(df8)
