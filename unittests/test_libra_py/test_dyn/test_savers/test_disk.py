from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.observables import ObservableConfig
from libra_py.dyn.savers import (
    FaultTolerantSaver,
    HDF5Saver,
    HDF5TimeSeriesSaver,
    JSONLinesSaver,
    load_hdf5_steps,
    load_hdf5_timeseries,
    load_json_steps,
    load_saved_steps,
    stack_saved_observables,
)


OUTPUT_ROOT = ROOT / "unittests" / "test_libra_py" / "test_dyn" / "test_savers" / "output"


def test_fault_tolerant_saver_writes_manifest_after_npz():
    output_dir = _prepare_output_dir("npz_steps")
    saver = FaultTolerantSaver(output_dir=output_dir)

    path0 = saver.save_step(0, {"time": 0.0, "populations": np.array([1.0, 0.0])})
    path1 = saver.save_step(1, {"time": 0.5, "populations": np.array([0.8, 0.2])})

    assert path0.exists()
    assert path1.exists()
    assert (output_dir / "manifest.jsonl").exists()

    records = load_saved_steps(output_dir)
    assert len(records) == 2
    assert records[0]["data"]["time"] == 0.0
    np.testing.assert_allclose(records[1]["data"]["populations"], [0.8, 0.2])


def test_fault_tolerant_saver_write_mode_starts_fresh_series():
    output_dir = _prepare_output_dir("npz_write_mode")
    FaultTolerantSaver(output_dir=output_dir).save_step(7, {"time": 7.0})

    saver = FaultTolerantSaver(output_dir=output_dir, mode="w")
    saver.save_step(0, {"time": 0.0})

    assert not (output_dir / "step_00000007.npz").exists()
    records = load_saved_steps(output_dir)
    assert [record["record"]["step"] for record in records] == [0]


def test_save_observables_computes_and_saves_selected_fields():
    output_dir = _prepare_output_dir("npz_observables")
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.p[0, 0] = [2.0]
    storage.iM[0, 0] = [0.5]
    storage.ham_adi[0, 0] = np.diag([0.1, 0.2])

    saver = FaultTolerantSaver(output_dir=output_dir)
    saver.save_observables(
        storage,
        traj,
        ObservableConfig(energies=True, populations=True),
        step=3,
        time=1.5,
    )

    records = load_saved_steps(output_dir)
    assert records[0]["record"]["step"] == 3
    assert records[0]["data"]["step"] == 3
    np.testing.assert_allclose(records[0]["data"]["mean_populations"], [1.0, 0.0])
    np.testing.assert_allclose(records[0]["data"]["kinetic_energy"], [1.0])


def test_hdf5_saver_writes_complete_step_groups():
    output_dir = _prepare_output_dir("hdf_steps")

    with HDF5Saver(output_dir=output_dir, filename="data.hdf") as saver:
        saver.save_step(
            0,
            {"time": 0.0, "populations": np.array([1.0, 0.0])},
            metadata={"kind": "initial"},
        )
        saver.save_step(
            1,
            {"time": 0.5, "populations": np.array([0.8, 0.2])},
        )

    records = load_hdf5_steps(output_dir=output_dir, filename="data.hdf")

    assert len(records) == 2
    assert records[0]["record"]["step"] == 0
    assert records[0]["record"]["metadata"] == {"kind": "initial"}
    assert records[1]["data"]["time"] == 0.5
    np.testing.assert_allclose(records[1]["data"]["populations"], [0.8, 0.2])


def test_hdf5_saver_saves_observables():
    output_dir = _prepare_output_dir("hdf_observables")
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.p[0, 0] = [2.0]
    storage.iM[0, 0] = [0.5]
    storage.ham_adi[0, 0] = np.diag([0.1, 0.2])

    with HDF5Saver(output_dir=output_dir, filename="observables.hdf") as saver:
        saver.save_observables(
            storage,
            traj,
            ObservableConfig(energies=True, populations=True),
            step=4,
            time=2.0,
        )

    records = load_hdf5_steps(output_dir=output_dir, filename="observables.hdf")

    assert records[0]["record"]["step"] == 4
    assert records[0]["data"]["step"] == 4
    np.testing.assert_allclose(records[0]["data"]["mean_populations"], [1.0, 0.0])
    np.testing.assert_allclose(records[0]["data"]["kinetic_energy"], [1.0])


def test_hdf5_saver_accepts_directory_path_without_filename():
    output_dir = _prepare_output_dir("hdf_directory_path")

    with HDF5Saver(output_dir) as saver:
        saver.save_step(0, {"time": 0.0, "populations": np.array([1.0, 0.0])})

    assert (output_dir / "data.hdf").exists()
    assert len(load_hdf5_steps(output_dir)) == 1


def test_hdf5_saver_accepts_less_frequent_flush_options():
    output_dir = _prepare_output_dir("hdf_less_frequent_flush")

    with HDF5Saver(
        output_dir=output_dir,
        filename="less_flush.hdf",
        durable=False,
        flush_stride=10,
    ) as saver:
        saver.save_step(0, {"time": 0.0})
        saver.save_step(1, {"time": 0.5})

    records = load_hdf5_steps(output_dir=output_dir, filename="less_flush.hdf")
    assert len(records) == 2
    np.testing.assert_allclose([record["data"]["time"] for record in records], [0.0, 0.5])


def test_hdf5_timeseries_saver_appends_chunked_datasets():
    output_dir = _prepare_output_dir("hdf_timeseries")

    with HDF5TimeSeriesSaver(output_dir=output_dir, filename="timeseries.hdf") as saver:
        saver.save_step(0, {"time": 0.0, "populations": np.array([1.0, 0.0])})
        saver.save_step(1, {"time": 0.5, "populations": np.array([0.8, 0.2])})

    data = load_hdf5_timeseries(output_dir=output_dir, filename="timeseries.hdf")
    np.testing.assert_allclose(data["time"], [0.0, 0.5])
    np.testing.assert_allclose(data["populations"], [[1.0, 0.0], [0.8, 0.2]])


def test_json_lines_saver_writes_human_readable_records():
    output_dir = _prepare_output_dir("json_lines")

    with JSONLinesSaver(output_dir=output_dir, filename="observables.jsonl") as saver:
        saver.save_step(0, {"time": 0.0, "populations": np.array([1.0, 0.0])})
        saver.save_step(1, {"time": 0.5, "populations": np.array([0.8, 0.2])})

    records = load_json_steps(output_dir=output_dir, filename="observables.jsonl")
    assert len(records) == 2
    assert records[1]["data"]["time"] == 0.5
    assert records[1]["data"]["populations"] == [0.8, 0.2]


def test_stack_saved_observables_handles_loaded_records():
    records = [
        {"data": {"time": 0.0, "populations": np.array([1.0, 0.0]), "q": np.array([[0.0]])}},
        {"data": {"time": 0.5, "populations": np.array([0.8, 0.2]), "q": np.array([[1.0]])}},
    ]

    stacked = stack_saved_observables(records, keys=("time", "populations", "q"))

    np.testing.assert_allclose(stacked["time"], [0.0, 0.5])
    np.testing.assert_allclose(stacked["populations"], [[1.0, 0.0], [0.8, 0.2]])
    assert stacked["q"].shape == (2, 1, 1)


def _prepare_output_dir(name):
    output_dir = OUTPUT_ROOT / name
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("manifest.jsonl", "step_*.npz", "*.hdf", "*.jsonl", ".step_*.tmp"):
        for path in output_dir.glob(pattern):
            path.unlink()
    return output_dir
