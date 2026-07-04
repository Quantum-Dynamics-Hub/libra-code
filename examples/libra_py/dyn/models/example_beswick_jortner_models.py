"""
Calculate and plot the Beswick-Jortner ICN model.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/models/example_beswick_jortner_models.py
"""

from pathlib import Path
import os

OUT_DIR = Path("beswick_jortner_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")

import numpy as np

from libra_py.dyn.models import BeswickJortnerModel
from libra_py.units import Angst

from _example_utils import compute_curves, plot_model_curves


def main():
    r = np.linspace(0.95, 1.60, 501) * Angst
    R = np.full_like(r, 2.50 * Angst)
    q = np.vstack((r, R))
    curves = compute_curves(BeswickJortnerModel(), q)
    data_path, plot_path = plot_model_curves(
        "beswick_jortner",
        r / Angst,
        curves,
        OUT_DIR,
    )
    print(f"beswick_jortner: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
