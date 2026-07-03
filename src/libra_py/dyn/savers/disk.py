"""
Fault-tolerant on-disk saving for dynamics snapshots.

Each step is written to a temporary NPZ file, fsynced, and atomically renamed.
A JSON-lines manifest is appended after the step file is durable. If a run
crashes, all manifest records point to complete step files.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover - exercised only without optional h5py
    h5py = None


class FaultTolerantSaver:
    """
    Save selected arrays step-by-step with atomic file replacement.

    Parameters
    ----------
    path
        Output directory. Kept for backward compatibility.
    output_dir
        Explicit output directory. If supplied, it takes precedence over path.
    compressed
        Use ``np.savez_compressed``. Leave False for faster writes.
    manifest_name
        JSON-lines manifest filename inside ``path``.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        output_dir: str | Path | None = None,
        compressed: bool = False,
        manifest_name: str = "manifest.jsonl",
    ):
        if path is None and output_dir is None:
            raise ValueError("Either path or output_dir must be supplied")
        self.path = Path(output_dir if output_dir is not None else path)
        self.output_dir = self.path
        self.compressed = compressed
        self.manifest_path = self.path / manifest_name
        self.path.mkdir(parents=True, exist_ok=True)

    def save_step(
        self,
        step: int,
        data: dict[str, Any],
        metadata: dict[str, Any] | None = None,
    ) -> Path:
        """Atomically save one step and append a durable manifest record."""

        filename = f"step_{int(step):08d}.npz"
        final_path = self.path / filename
        tmp_path = self.path / f".{filename}.tmp"

        arrays = {
            key: _as_saveable_array(value)
            for key, value in data.items()
            if _is_array_payload(value)
        }
        scalars = {
            key: _as_json_value(value)
            for key, value in data.items()
            if not _is_array_payload(value)
        }
        arrays["__scalars__"] = np.array(json.dumps(scalars, sort_keys=True))

        opener = np.savez_compressed if self.compressed else np.savez
        with open(tmp_path, "wb") as handle:
            opener(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, final_path)
        _fsync_dir(self.path)

        record = {
            "step": int(step),
            "file": filename,
            "compressed": bool(self.compressed),
            "metadata": metadata or {},
        }
        self._append_manifest(record)
        return final_path

    def save_observables(
        self,
        storage: Any,
        traj: Any,
        config: Any = None,
        step: int | None = None,
        time: float | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> Path:
        """Compute observables and save them as one durable step record."""

        from ..observables import ObservableConfig, compute_observables

        config = config or ObservableConfig()
        step_value = storage.timestep if step is None else step
        data = compute_observables(storage, traj, config, step=step_value, time=time)
        return self.save_step(step_value, data, metadata=metadata)

    def _append_manifest(self, record: dict[str, Any]) -> None:
        line = json.dumps(record, sort_keys=True) + "\n"
        with open(self.manifest_path, "a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())


class HDF5Saver:
    """
    Save dynamics snapshots into one HDF5 file.

    Each step is stored under ``/steps/step_XXXXXXXX`` and the file is flushed
    after every step. A step group is marked complete only after all datasets
    and metadata have been written, so readers can ignore incomplete groups if
    a process stops mid-write.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        output_dir: str | Path | None = None,
        filename: str = "data.hdf",
        mode: str = "a",
        compression: str | None = None,
        compression_opts: int | None = None,
    ):
        if h5py is None:
            raise ImportError("HDF5Saver requires the optional 'h5py' package")
        if path is None and output_dir is None:
            raise ValueError("Either path or output_dir must be supplied")
        self.path = _resolve_hdf5_path(path, output_dir, filename)
        self.output_dir = self.path.parent
        self.filename = self.path.name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.compression = compression
        self.compression_opts = compression_opts
        self.file = h5py.File(self.path, mode)
        self.steps = self.file.require_group("steps")
        self.file.attrs.setdefault("format", "libra_py.dyn.savers.hdf5")

    def save_step(
        self,
        step: int,
        data: dict[str, Any],
        metadata: dict[str, Any] | None = None,
    ) -> Path:
        """Save one complete step group and flush it to disk."""

        name = f"step_{int(step):08d}"
        if name in self.steps:
            del self.steps[name]

        group = self.steps.create_group(name)
        group.attrs["complete"] = False
        group.attrs["step"] = int(step)
        group.attrs["metadata_json"] = json.dumps(metadata or {}, sort_keys=True)

        scalars = {}
        for key, value in data.items():
            if _is_array_payload(value):
                self._write_dataset(group, key, _as_saveable_array(value))
            else:
                scalars[key] = _as_json_value(value)
        group.attrs["scalars_json"] = json.dumps(scalars, sort_keys=True)
        group.attrs["complete"] = True

        self.file.flush()
        _fsync_file(self.file)
        return self.path

    def save_observables(
        self,
        storage: Any,
        traj: Any,
        config: Any = None,
        step: int | None = None,
        time: float | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> Path:
        """Compute observables and save them as one complete HDF5 step."""

        from ..observables import ObservableConfig, compute_observables

        config = config or ObservableConfig()
        step_value = storage.timestep if step is None else step
        data = compute_observables(storage, traj, config, step=step_value, time=time)
        return self.save_step(step_value, data, metadata=metadata)

    def close(self) -> None:
        """Flush and close the HDF5 file."""

        if self.file:
            self.file.flush()
            _fsync_file(self.file)
            self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def _write_dataset(self, group, key: str, value: Any) -> None:
        kwargs = {}
        if self.compression is not None and np.asarray(value).shape != ():
            kwargs["compression"] = self.compression
            if self.compression_opts is not None:
                kwargs["compression_opts"] = self.compression_opts
        group.create_dataset(key, data=value, **kwargs)


def load_saved_steps(path: str | Path, manifest_name: str = "manifest.jsonl") -> list[dict[str, Any]]:
    """Load all complete step records listed in the manifest."""

    path = Path(path)
    manifest = path / manifest_name
    records = []
    if not manifest.exists():
        return records

    for record in _manifest_records(manifest):
        step_path = path / record["file"]
        if not step_path.exists():
            continue
        with np.load(step_path, allow_pickle=False) as payload:
            data = {key: payload[key] for key in payload.files if key != "__scalars__"}
            scalars = json.loads(str(payload["__scalars__"]))
        data.update(scalars)
        records.append({"record": record, "data": data})
    return records


def load_hdf5_steps(
    path: str | Path | None = None,
    output_dir: str | Path | None = None,
    filename: str = "data.hdf",
) -> list[dict[str, Any]]:
    """Load complete step groups from an HDF5 dynamics file."""

    if h5py is None:
        raise ImportError("load_hdf5_steps requires the optional 'h5py' package")

    path = _resolve_hdf5_path(path, output_dir, filename)
    records = []
    with h5py.File(path, "r") as handle:
        steps = handle.get("steps")
        if steps is None:
            return records
        for name in sorted(steps):
            group = steps[name]
            if not bool(group.attrs.get("complete", False)):
                continue
            data = {key: group[key][()] for key in group.keys()}
            scalars = json.loads(group.attrs.get("scalars_json", "{}"))
            metadata = json.loads(group.attrs.get("metadata_json", "{}"))
            data.update(scalars)
            records.append({
                "record": {
                    "step": int(group.attrs["step"]),
                    "group": name,
                    "metadata": metadata,
                },
                "data": data,
            })
    return records


def _manifest_records(path: Path) -> Iterable[dict[str, Any]]:
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield json.loads(line)


def _resolve_hdf5_path(
    path: str | Path | None,
    output_dir: str | Path | None,
    filename: str,
) -> Path:
    if output_dir is not None:
        return Path(output_dir) / filename
    if path is None:
        raise ValueError("Either path or output_dir must be supplied")
    resolved = Path(path)
    if resolved.suffix:
        return resolved
    return resolved / filename


def _is_array_payload(value: Any) -> bool:
    return hasattr(value, "shape") or isinstance(value, (list, tuple))


def _as_saveable_array(value: Any) -> Any:
    return np.asarray(value)


def _as_json_value(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_as_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _as_json_value(val) for key, val in value.items()}
    return str(value)


def _fsync_dir(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _fsync_file(h5_file) -> None:
    try:
        fd = h5_file.id.get_vfd_handle()
    except Exception:
        return
    if isinstance(fd, int) and fd >= 0:
        os.fsync(fd)
