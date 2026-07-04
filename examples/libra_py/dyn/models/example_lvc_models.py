from pathlib import Path
import os

OUT_DIR = Path("lvc_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import LVCModel


def main():
    x = np.linspace(-4.0, 4.0, 801)
    y = np.zeros_like(x)
    model = LVCModel(params={"Delta1": 0.0, "Delta2": 0.1, "omega": [0.01, 0.02], "d1": [0.0, 0.0], "d2": [0.01, -0.01], "coup": [0.001, 0.002], "mass": [1.0, 2.0]})
    curves = compute_curves(model, np.asarray([x, y]))
    data_path, plot_path = plot_model_curves("lvc_model", x, curves, OUT_DIR)
    print(f"saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
