"""
Calculate and plot the three-state phenol model along an O-H scan.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/models/example_phenol_models.py
"""

from pathlib import Path
import os

OUT_DIR = Path("phenol_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")

import numpy as np

from libra_py.dyn.models import PhenolModel
from libra_py.units import Angst

from _example_utils import compute_curves, plot_model_curves


def main():
    r = np.linspace(0.80, 1.60, 601) * Angst
    theta = np.full_like(r, 0.60)
    q = np.vstack((r, theta))
    curves = compute_curves(PhenolModel(), q)
    data_path, plot_path = plot_model_curves(
        "phenol_theta_0p60",
        r / Angst,
        curves,
        OUT_DIR,
    )
    print(f"phenol_theta_0p60: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
