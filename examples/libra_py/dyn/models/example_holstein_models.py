from pathlib import Path
import os

OUT_DIR = Path("holstein_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import Holstein2Model, Holstein3Model, Holstein4Model, Holstein5Model


def main():
    q = np.linspace(-2.0, 4.0, 801)
    params = {"E_n": [0.0, 0.001, 0.002], "x_n": [0.0, 1.0, 2.0], "k_n": [0.01, 0.012, 0.014]}
    models = {
        "holstein2_model": Holstein2Model(params={**params, "V": 0.001}),
        "holstein3_model": Holstein3Model(params={**params, "V_n": [0.001, 0.0002]}),
        "holstein4_model": Holstein4Model(params={**params, "V": [[0.0, 0.001, 0.0002], [0.001, 0.0, 0.001], [0.0002, 0.001, 0.0]]}),
        "holstein5_model": Holstein5Model(params={**params, "V": [[0.0, 0.001, 0.0002], [0.001, 0.0, 0.001], [0.0002, 0.001, 0.0]], "alpha": [[0.0, 0.2, 0.1], [0.2, 0.0, 0.2], [0.1, 0.2, 0.0]], "x_nm": [[0.0, 0.5, 1.0], [0.5, 0.0, 1.5], [1.0, 1.5, 0.0]]}),
    }
    for name, model in models.items():
        curves = compute_curves(model, q)
        data_path, plot_path = plot_model_curves(name, q, curves, OUT_DIR)
        print(f"{name}: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
