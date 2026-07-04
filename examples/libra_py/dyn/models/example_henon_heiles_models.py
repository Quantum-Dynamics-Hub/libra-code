from pathlib import Path
import os

OUT_DIR = Path("henon_heiles_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import HenonHeilesModel


def main():
    x = np.linspace(-3.0, 3.0, 801)
    y = np.zeros_like(x)
    curves = compute_curves(HenonHeilesModel(), np.asarray([x, y]))
    data_path, plot_path = plot_model_curves("henon_heiles", x, curves, OUT_DIR)
    print(f"saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
