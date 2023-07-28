from ocelot.cpbd.io import load_particle_array
from importlib_resources import files

from .converter import injector_to_astra

ASTRA_FILES_DIR = files("oxfel.astra")


FILES = {10_000: ASTRA_FILES_DIR / "euxfel-injector-10k.0320.001",
         100_000: ASTRA_FILES_DIR / "euxfel-injector-100k.0320.001"}

def load_reference_0320_100k_distribution():
    return load_reference_0320_distribution(100_000)


def load_reference_0320_10k_distribution():
    return load_reference_0320_distribution(10_000)

def load_reference_0320_distribution(nint):
    return load_particle_array(FILES[nint])
