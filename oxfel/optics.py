from .longlist import make_default_longlist


FIXED_MATCH_POINTS = [
    "MATCH.52.I1",
    "MATCH.73.I1",
    "MATCH.104.I1",
    "MATCH.218.B1",
    "MATCH.446.B2",
]


_TWISS_NAMES = ["beta_x", "alpha_x", "beta_y", "alpha_y"]
_MAD8_TWISS_NAMES = ["BETX", "ALFX", "BETY", "ALFY"]


def get_name2_fixed_match_points(additional_names=None):
    if additional_names is None:
        additional_names = []
    ll = make_default_longlist()
    return [
        ll.name1_to_name2(name1).item()
        for name1 in FIXED_MATCH_POINTS + additional_names
    ]


def get_match_point_optics(twiss_df, additional_names=None):
    if additional_names is None:
        additional_names = []
    try:
        return twiss_df.set_index("id").loc[FIXED_MATCH_POINTS + additional_names][
            _TWISS_NAMES
        ]
    except KeyError:
        name2s = get_name2_fixed_match_points(additional_names=additional_names)
        return twiss_df.set_index("NAME").loc[name2s][_MAD8_TWISS_NAMES]
