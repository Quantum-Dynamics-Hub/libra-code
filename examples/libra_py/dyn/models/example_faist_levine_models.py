from pathlib import Path
import os

OUT_DIR = Path("faist_levine_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import FaistLevineModel, faist_levine_lii_params, faist_levine_nai_params


def main():
    q = np.linspace(3.0, 12.0, 801)
    for name, model in {
        "faist_levine_nai": FaistLevineModel(params=faist_levine_nai_params()),
        "faist_levine_lii": FaistLevineModel(params=faist_levine_lii_params()),
    }.items():
        curves = compute_curves(model, q)
        data_path, plot_path = plot_model_curves(name, q, curves, OUT_DIR)
        print(f"{name}: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
