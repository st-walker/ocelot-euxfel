import glob
from pathlib import Path

def find_wake_file(filename):
    here = Path(__file__).parent
    for path in here.rglob("*.dat"):
        if path.name == filename:
            return here / path

    raise FileNotFoundError(filename)
