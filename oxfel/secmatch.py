from __future__ import annotations

from .tracking import start_sim_to_match_37
from .optics import MATCH_37, MATCH_52, INJECTOR_MATCHING_QUAD_NAMES, default_match_point_optics
from .fel_track import EuXFELSimConfig
from .matching import get_unary_twiss_function, match_with_backtracking


def match_injector(fel, parray032, felconfig=None, match="projected"):
    parray37 = start_sim_to_match_37(fel, parray032, felconfig=felconfig)
    navi = fel.to_navigator(start=MATCH_37, stop=MATCH_52)

    match52 = default_match_point_optics().set_index("id").loc["MATCH.52.I1"]

    goal_twiss = Twiss(**dict(alpha_x=match52.alpha_x,
                              alpha_y=match52.alpha_y,
                              beta_x=match52.beta_x,
                              beta_y=match52.beta_y,
                              id="MATCH.52.I1"))

    strengths, mismatch = match_with_backtracking(navi, parray37,
                                                  goal_twiss,
                                                  INJECTOR_MATCHING_QUAD_NAMES,
                                                  maxiter=1,
                                                  match=match)

    if felconfig is None:
        felconfig = EuXFELSimConfig()

    felconfig.update_components(INJECTOR_MATCHING_QUAD_NAMES, strengths, "k1")

    return felconfig, mismatch

def get_injector_matching_convergence(fel, parray032, felconfig=None, match="projected"):
    parray37 = start_sim_to_match_37(fel, parray032, felconfig=felconfig)
    navi = fel.to_navigator(start=MATCH_37, stop=MATCH_52)
    match52 = default_match_point_optics().set_index("id").loc["MATCH.52.I1"]
    goal_twiss = Twiss(**dict(alpha_x=match52.alpha_x,
                              alpha_y=match52.alpha_y,
                              beta_x=match52.beta_x,
                              beta_y=match52.beta_y,
                              id="MATCH.52.I1"))

    twiss_function = get_unary_twiss_function(match)
    matcher = BacktrackingLinearMatcher(navi,
                                        parray37,
                                        goal_twiss,
                                        INJECTOR_MATCHING_QUAD_NAMES,
                                        twiss_function=twiss_function)

    matcher.initial_match()
    matcher.track_forwards()

    bmags = matcher.bmags()

    bmagxs = [bmags[0]]
    bmagys = [bmags[1]]
    l2losses = [matcher.l2loss()]


    for i in range(5):
        from IPython import embed; embed()
        matcher.rematch()

        bmags = matcher.bmags()
        l2loss = matcher.l2loss()

        bmagxs.append(bmags[0])
        bmagys.append(bmags[1])
        l2losses.append(l2loss)

    return bmagxs, bmagys, l2losses




# def match_injector_slice(fel, parray032, felconfig=None, match="emax"):
#     parray37 = start_sim_to_match_37(fel, parray032, felconfig=felconfig)
#     navi = fel.to_navigator(start=MATCH_37, stop=MATCH_52)

#     match52 = default_match_point_optics().set_index("id").loc["MATCH.52.I1"]

#     goal_twiss = Twiss(**dict(alpha_x=match52.alpha_x,
#                               alpha_y=match52.alpha_y,
#                               beta_x=match52.beta_x,
#                               beta_y=match52.beta_y,
#                               id="MATCH.52.I1"))

#     # from IPython import embed; embed()

#     strengths, mismatch = match_with_backtracking(navi, parray37,
#                                                   goal_twiss,
#                                                   INJECTOR_MATCHING_QUAD_NAMES,
#                                                   maxiter=4
#                                                   match=match)

#     if felconfig is None:
#         felconfig = FELSimulationConfig()

#     felconfig.update_components(INJECTOR_MATCHING_QUAD_NAMES, strengths, "k1")

#     return felconfig, mismatch
