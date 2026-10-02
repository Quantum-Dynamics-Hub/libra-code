"""Manuscript DAC parameters. Edit here; run.py needs no arguments."""
from pathlib import Path
from math import sqrt

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
CASES = {"": {
    "title": "Tully 2: dual avoided crossing", "model": "tully2",
    "ldr_bounds": [[-20.0, 22.0]], "ldr_spacing": [0.05],
    "dvr_bounds": [[-20.0, 22.0]], "dvr_spacing": [0.05],
    "mass": [2000.0], "q0": [-6.0], "p0": [30.0],
    "sigma": [1 / sqrt(2)], "istate": 0,
    "dt": 1.0, "duration": 1500.0, "save_every": 10,
    "snapshots": [0.0, 500.0, 1000.0, 1500.0],
}}
