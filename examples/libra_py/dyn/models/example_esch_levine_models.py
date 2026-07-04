from pathlib import Path
import os

OUT_DIR = Path("esch_levine_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import EschLevineJCP2020Model, EschLevineLinearModel


def main():
    q = np.linspace(-2.0, 2.0, 801)
    models = {
        "esch_levine_linear": EschLevineLinearModel(params={"nstates": 2, "V": [[0.0, 0.005], [0.005, 0.0]], "w": [[-0.1, 0.0], [0.0, 0.1]]}),
        "esch_levine_jcp2020": EschLevineJCP2020Model(params={"nstates": 4}),
    }
    for name, model in models.items():
        curves = compute_curves(model, q)
        data_path, plot_path = plot_model_curves(name, q, curves, OUT_DIR)
        print(f"{name}: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
