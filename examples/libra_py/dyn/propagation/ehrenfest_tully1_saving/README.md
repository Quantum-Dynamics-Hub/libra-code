# Ehrenfest Tully1 Saving Examples

Run any script from the repository root, for example:

```bash
PYTHONPATH=src python examples/libra_py/dyn/propagation/ehrenfest_tully1_saving/01_hdf5_stride.py
```

Use a shorter run for quick checks:

```bash
EHRENFEST_TULLY1_NSTEPS=30 PYTHONPATH=src python examples/libra_py/dyn/propagation/ehrenfest_tully1_saving/03_hdf5_chunked_timeseries.py
```

Each script writes to its own top-level output directory, such as
`01_hdf5_stride_outputs` or `06_json_saver_outputs`.

The examples demonstrate:

- `01_hdf5_stride.py`: step-wise HDF5 with save stride.
- `02_hdf5_less_frequent_flush.py`: step-wise HDF5 with less frequent flush/fsync.
- `03_hdf5_chunked_timeseries.py`: chunked extendable HDF5 datasets.
- `04_hdf5_fewer_fields.py`: save fewer observable fields.
- `05_npz_saver.py`: per-step NPZ files plus manifest.
- `06_json_saver.py`: JSON-lines summary output.
