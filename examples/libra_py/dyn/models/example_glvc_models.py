from pathlib import Path
import os

OUT_DIR = Path("glvc_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import GLVCModel


def main():
    x = np.linspace(-4.0, 4.0, 801)
    y = np.zeros_like(x)
    model = GLVCModel(params={"nstates": 2, "num_osc": 2, "Ham": [[0.0, 0.002], [0.002, 0.01]], "omega": [[0.01, 0.02], [0.015, 0.025]], "coupl": [[0.001, -0.002], [-0.001, 0.003]], "mass": [1.0, 2.0]})
    curves = compute_curves(model, np.asarray([x, y]))
    data_path, plot_path = plot_model_curves("glvc_model", x, curves, OUT_DIR)
    print(f"saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
