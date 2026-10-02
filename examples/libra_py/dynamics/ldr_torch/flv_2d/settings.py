"""Manuscript regular-grid FLV parameters for both coupling strengths."""
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
COMMON = {
    "model": "flv",
    "ldr_bounds": [[1.0, 7.0], [-1.5, 1.5]], "ldr_spacing": [0.05, 0.05],
    # Larger DVR box used in the archived reference calculations.
    "dvr_bounds": [[0.0, 10.0], [-2.0, 2.0]], "dvr_spacing": [0.05, 0.05],
    "mass": [20000.0, 6667.0], "q0": [2.0, 0.0], "p0": [0.0, 0.0],
    "sigma": [0.150, 0.197], "istate": 1,
    "dt": 1.0, "duration": 3000.0, "save_every": 10,
    "snapshots": [0.0, 1000.0, 2000.0],
}
CASES = {
    "weak": dict(COMMON, title="FLV: weak coupling", gamma=0.01),
    "strong": dict(COMMON, title="FLV: strong coupling", gamma=0.08),
}
