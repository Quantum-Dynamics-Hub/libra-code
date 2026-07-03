"""
Educational example: fault-tolerant step saving into an explicit folder.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/savers/example_saver.py
"""

from pathlib import Path

import numpy as np

from libra_py.dyn.savers import (
    FaultTolerantSaver,
    HDF5Saver,
    load_hdf5_steps,
    load_saved_steps,
)


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def main():
    _clean_example_outputs()
    npz_saver = FaultTolerantSaver(output_dir=OUTPUT_DIR)

    npz_saver.save_step(
        0,
        {
            "time": 0.0,
            "populations": np.array([1.0, 0.0]),
            "total_energy": np.array([0.25]),
        },
    )
    npz_saver.save_step(
        1,
        {
            "time": 0.5,
            "populations": np.array([0.9, 0.1]),
            "total_energy": np.array([0.251]),
        },
    )

    with HDF5Saver(output_dir=OUTPUT_DIR, filename="data.hdf") as hdf_saver:
        hdf_saver.save_step(
            0,
            {
                "time": 0.0,
                "populations": np.array([1.0, 0.0]),
                "total_energy": np.array([0.25]),
            },
        )
        hdf_saver.save_step(
            1,
            {
                "time": 0.5,
                "populations": np.array([0.9, 0.1]),
                "total_energy": np.array([0.251]),
            },
        )

    print("npz records:", len(load_saved_steps(OUTPUT_DIR)))
    print("hdf records:", len(load_hdf5_steps(output_dir=OUTPUT_DIR, filename="data.hdf")))
    print("directory:", OUTPUT_DIR)


def _clean_example_outputs():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for pattern in ("manifest.jsonl", "step_*.npz", "data.hdf", ".step_*.tmp"):
        for path in OUTPUT_DIR.glob(pattern):
            path.unlink()


if __name__ == "__main__":
    main()
