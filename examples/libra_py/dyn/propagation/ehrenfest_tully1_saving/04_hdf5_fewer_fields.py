"""
Ehrenfest Tully1 ensemble: save fewer fields for smaller output.
"""

from common import (
    add_mean_energies,
    compact_observable_config,
    make_output_dir,
    plot_energy_and_populations,
    print_summary,
    run_with_saver,
    save_npz_copy,
)
from libra_py.dyn.savers import HDF5TimeSeriesSaver, load_hdf5_timeseries


out_dir = make_output_dir("04_hdf5_fewer_fields")
filename = "observables_compact.hdf"

with HDF5TimeSeriesSaver(output_dir=out_dir, filename=filename, mode="w") as saver:
    run_with_saver(saver, compact_observable_config())

observables = load_hdf5_timeseries(output_dir=out_dir, filename=filename)
add_mean_energies(observables)
saved = [out_dir / filename, save_npz_copy(out_dir, "observables_compact.npz", observables)]
saved += plot_energy_and_populations(out_dir, observables, prefix="compact")
print_summary("Compact HDF5 saving", out_dir, observables, saved)
