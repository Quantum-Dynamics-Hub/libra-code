# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Explicit QTAG persistence; no files are created by the dynamics runner.

NPZ is the dependency-free default. HDF5 is optional (requires h5py). Arrays
are detached/copied to CPU only here; the live result and its graphs survive.
The archive contains time, Etot, dEtot, norm, pops, coeffs, q, p, a, s,
states and target_states. No pickle or C++ objects are serialized.
"""
from pathlib import Path
import numpy as np
from ._backend import to_numpy


def save_qtag_data(result, path):
    """Write result.history (or a history dict) to .npz, .h5, or .hdf5.

    Creates parent directories and returns the resolved output Path. Existing
    output at the explicit path is replaced. Full precision is retained,
    including complex coefficients and per-DOF phases.
    """
    path = Path(path)
    history = result.history if hasattr(result, "history") else result
    arrays = {k: to_numpy(v) for k, v in history.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".npz":
        np.savez_compressed(path, **arrays)
    elif path.suffix in (".h5", ".hdf5"):
        import h5py
        with h5py.File(path, "w") as f:
            f.attrs["format"] = "libra-qtag-1"
            for k, v in arrays.items():
                f.create_dataset(k, data=v)
    else:
        raise ValueError("output extension must be .npz, .h5 or .hdf5")
    return path.resolve()


def load_qtag_data(path, backend="numpy", device="cpu"):
    """Load a history dictionary; optionally convert numeric arrays to tensors."""
    path = Path(path)
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as data:
            arrays = {k: data[k] for k in data.files}
    elif path.suffix in (".h5", ".hdf5"):
        import h5py
        with h5py.File(path, "r") as data:
            arrays = {k: data[k][()] for k in data}
    else:
        raise ValueError("input extension must be .npz, .h5 or .hdf5")
    if backend == "torch":
        import torch
        return {k: torch.as_tensor(v, device=device) for k, v in arrays.items()}
    if backend != "numpy":
        raise ValueError("backend must be numpy or torch")
    return arrays


def save_text(result, directory):
    """Export readable time, energy, norm and populations (one row per time).

    Use NPZ/HDF5 for restarting or reconstructing wavefunctions; flattened
    legacy text layouts are deliberately not used for multidimensional data.
    """
    h = result.history if hasattr(result, "history") else result
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    data = np.column_stack([to_numpy(h[k]) for k in ("time", "Etot", "norm", "pops")])
    header = "time Etot norm "+" ".join(f"pop_{int(n)}" for n in to_numpy(h["target_states"]))
    np.savetxt(directory/"observables.txt", data, header=header)
    return directory/"observables.txt"
