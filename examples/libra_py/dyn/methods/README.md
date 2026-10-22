# DynamicsEngine method examples

| Folder | Model | Engine selector |
| --- | --- | --- |
| `1_qtag/1_basics` | Holstein / Ferretti | `qtag` |
| `2_tsh/1_basics` | Tully-1 | `tsh` (`tsh_method=0`) |
| `3_ehrenfest/1_basics` | Tully-1 | `ehrenfest` |

Each basics directory provides `run.py`: configure parameters, initialize TBFs,
run DynamicsEngine, save snapshots, and write individual PNG figures. Run with
`PYTHONPATH=/absolute/path/to/libra/src` from any directory. Shared CLI options
are `--steps`, `--output`, and `--save-format npz|hdf5`. QTAG additionally
supports NumPy/PyTorch and model selection; Tully examples use NumPy and expose
`--ntraj` and `--save-stride`.

Default outputs live beside each script in `output/<model>_<backend>`.
All basics examples print method/backend/model, completed steps, final
populations, norm, energy change, and output location, and save energy,
energy-change, norm, population, position, and momentum plots. QTAG adds
widths, phases, and density snapshots; FSSH adds active-state diagnostics.
QTAG populations are coherent diabatic populations; Tully populations are
adiabatic ensemble averages, with hopping fractions reported separately.

`2_tsh/2_decoherence` retains the extended FSSH/GFSH/FSSH2 comparison with
coherent, ID-A, and SDM dynamics. Its self-contained `run.py` runs all cases and plots them;
`--compute-only` and `--plot-only` select either stage separately.
