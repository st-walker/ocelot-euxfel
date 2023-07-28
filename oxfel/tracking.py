from .optics import START_SIM, MATCH_37, MATCH_52


def start_sim_to_match_37(fel, parray032, felconfig=None):
    result37 = fel.track(parray032,
                         start=START_SIM,
                         stop=MATCH_37,
                         felconfig=felconfig)
    return result37.parray1
