from pathlib import Path
import os

OUT_DIR = Path("zhu_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import ZhuDualLZSModel, ZhuDualRZDModel


def main():
    q = np.linspace(-8.0, 8.0, 801)
    for name, model in {"zhu_dual_rzd": ZhuDualRZDModel(), "zhu_dual_lzs": ZhuDualLZSModel()}.items():
        curves = compute_curves(model, q)
        data_path, plot_path = plot_model_curves(name, q, curves, OUT_DIR)
        print(f"{name}: saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
