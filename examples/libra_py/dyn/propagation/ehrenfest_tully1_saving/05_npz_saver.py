"""
Ehrenfest Tully1 ensemble: per-step NPZ files with manifest.
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
from libra_py.dyn.savers import FaultTolerantSaver, load_saved_steps, stack_saved_observables


out_dir = make_output_dir("05_npz_saver")

saver = FaultTolerantSaver(output_dir=out_dir, compressed=False)
run_with_saver(saver, full_observable_config(), save_stride=10)

records = load_saved_steps(out_dir)
observables = stack_saved_observables(records)
add_mean_energies(observables)
saved = [out_dir / "manifest.jsonl", save_npz_copy(out_dir, "observables_from_npz_steps.npz", observables)]
saved += plot_energy_and_populations(out_dir, observables, prefix="npz")
saved.append(plot_phase_portrait(out_dir, observables, prefix="npz"))
print_summary("NPZ step saving", out_dir, observables, saved)
