from ocelot.cpbd.io import load_particle_array
from importlib_resources import files


ASTRA_FILES_DIR = files("oxfel.astra")

def load_reference_0320_100k_distribution():
    return load_particle_array(ASTRA_FILES_DIR / "euxfel-injector-100k.0320.001")


def load_reference_0320_10k_distribution():
    return load_particle_array(ASTRA_FILES_DIR / "euxfel-injector-10k.0320.001")
