from ocelot.cpbd.io import load_particle_array


def load_reference_0320_distribution():
    return load_particle_array("/Users/stuartwalker/repos/oxfel/oxfel/astra/euxfel-injector-100k.0320.001")

def load_reference_0320_10k_distribution():
    return load_particle_array("/Users/stuartwalker/repos/oxfel/oxfel/astra/euxfel-injector-10k.0320.001")

