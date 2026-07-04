"""
Ehrenfest Tully1 ensemble: HDF5 saving with a stride.
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
from libra_py.dyn.savers import HDF5Saver, load_hdf5_steps, stack_saved_observables


out_dir = make_output_dir("01_hdf5_stride")
filename = "observables_stride.hdf"
save_stride = 10

with HDF5Saver(output_dir=out_dir, filename=filename, mode="w") as saver:
    run_with_saver(saver, full_observable_config(), save_stride=save_stride)

records = load_hdf5_steps(output_dir=out_dir, filename=filename)
observables = stack_saved_observables(records)
add_mean_energies(observables)
saved = [out_dir / filename, save_npz_copy(out_dir, "observables_stride.npz", observables)]
saved += plot_energy_and_populations(out_dir, observables, prefix="stride")
saved.append(plot_phase_portrait(out_dir, observables, prefix="stride"))
print_summary("HDF5 stride saving", out_dir, observables, saved)
