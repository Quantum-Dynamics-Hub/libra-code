"""
Ehrenfest Tully1 ensemble: chunked extendable HDF5 datasets.
"""

from common import (
    add_mean_energies,
    full_observable_config,
    make_output_dir,
    plot_energy_and_populations,
    plot_phase_portrait,
    print_summary,
    run_with_saver,
    save_npz_copy,
)
from libra_py.dyn.savers import HDF5TimeSeriesSaver, load_hdf5_timeseries


out_dir = make_output_dir("03_hdf5_chunked_timeseries")
filename = "observables_timeseries.hdf"

with HDF5TimeSeriesSaver(
    output_dir=out_dir,
    filename=filename,
    mode="w",
    flush_stride=100,
) as saver:
    run_with_saver(saver, full_observable_config())

observables = load_hdf5_timeseries(output_dir=out_dir, filename=filename)
add_mean_energies(observables)
saved = [out_dir / filename, save_npz_copy(out_dir, "observables_timeseries.npz", observables)]
saved += plot_energy_and_populations(out_dir, observables, prefix="timeseries")
saved.append(plot_phase_portrait(out_dir, observables, prefix="timeseries"))
print_summary("Chunked HDF5 time-series saving", out_dir, observables, saved)
