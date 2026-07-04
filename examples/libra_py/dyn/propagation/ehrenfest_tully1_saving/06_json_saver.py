"""
Ehrenfest Tully1 ensemble: JSON-lines summary saving.

JSON is best for compact summaries, so this example omits trajectory-resolved
coordinates and momenta.
"""

import numpy as np

from common import (
    add_mean_energies,
    compact_observable_config,
    make_output_dir,
    plot_energy_and_populations,
    print_summary,
    run_with_saver,
    save_npz_copy,
)
from libra_py.dyn.savers import JSONLinesSaver, load_json_steps


out_dir = make_output_dir("06_json_saver")
filename = "observables.jsonl"

with JSONLinesSaver(output_dir=out_dir, filename=filename, flush_stride=25) as saver:
    run_with_saver(saver, compact_observable_config(), save_stride=10)

records = load_json_steps(output_dir=out_dir, filename=filename)
observables = {
    key: np.stack([np.asarray(record["data"][key]) for record in records])
    for key in (
        "time",
        "populations",
        "mean_populations",
        "kinetic_energy",
        "potential_energy",
        "total_energy",
    )
}
add_mean_energies(observables)
saved = [out_dir / filename, save_npz_copy(out_dir, "observables_from_json.npz", observables)]
saved += plot_energy_and_populations(out_dir, observables, prefix="json")
print_summary("JSON-lines summary saving", out_dir, observables, saved)
