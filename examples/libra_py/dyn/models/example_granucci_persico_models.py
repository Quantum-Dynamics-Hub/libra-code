"""
Calculate and plot Granucci-Persico model Hamiltonian quantities.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/models/example_granucci_persico_models.py
"""

from pathlib import Path
import os

OUT_DIR = Path("granucci_persico_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")

import numpy as np

from libra_py.dyn.models import GranucciPersicoModel1, GranucciPersicoModel2

from _example_utils import compute_curves, plot_model_curves


def main():
    x = np.linspace(-5.0, 12.0, 801)

    curves1 = compute_curves(GranucciPersicoModel1(), np.asarray([x]))
    data_path, plot_path = plot_model_curves(
        "granucci_persico_model1",
        x,
        curves1,
        OUT_DIR,
    )
    print(f"granucci_persico_model1: saved {data_path} and {plot_path}")

    y = np.full_like(x, 0.5)
    curves2 = compute_curves(GranucciPersicoModel2(), np.vstack((x, y)))
    data_path, plot_path = plot_model_curves(
        "granucci_persico_model2_y_0p5",
        x,
        curves2,
        OUT_DIR,
    )
    print(f"granucci_persico_model2_y_0p5: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
