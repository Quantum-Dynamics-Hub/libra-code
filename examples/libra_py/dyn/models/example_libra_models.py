from pathlib import Path
import os

OUT_DIR = Path("libra_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import numpy as np

from _example_utils import compute_curves, plot_model_curves
from libra_py.dyn.models import LibraModel1


def main():
    q = np.linspace(-4.0, 4.0, 801)
    curves = compute_curves(LibraModel1(), q)
    data_path, plot_path = plot_model_curves("libra_model1", q, curves, OUT_DIR)
    print(f"saved {data_path} and {plot_path}")


if __name__ == "__main__":
    main()
