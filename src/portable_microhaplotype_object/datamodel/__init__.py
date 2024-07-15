from pathlib import Path
from .portable_microhaplotype_object import *

THIS_PATH = Path(__file__).parent

SCHEMA_DIRECTORY = THIS_PATH.parent / "schema"
MAIN_SCHEMA_PATH = SCHEMA_DIRECTORY / "portable_microhaplotype_object.yaml"
