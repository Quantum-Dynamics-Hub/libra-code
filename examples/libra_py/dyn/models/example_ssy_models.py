from pathlib import Path
import os

OUT_DIR = Path("ssy_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import SSYModel


def main():
    x = np.linspace(-4.0, 4.0, 801)
    y = np.zeros_like(x)
    curves = compute_curves(SSYModel(), np.asarray([x, y]))
    data_path, plot_path = plot_model_curves("ssy_model", x, curves, OUT_DIR)
    print(f"saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
