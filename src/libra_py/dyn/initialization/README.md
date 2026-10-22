# Ensemble initialization

`ensembles.py` owns shared nuclear sampling and both initialization paths:

- `sample_gaussian_initial_conditions`: returns `(npackets,ndof)` NumPy samples.
  Default `draw_order="trajectory"` preserves the original API and seeded
  results. Optional `draw_order="dof"` preserves QTAG's historical RNG ordering.
- `make_independent_trajectory_ensemble`: retains its existing arguments and
  `(storage, trajectory)` return value, normalizing electronic amplitudes per path.
- `make_gaussian_ensemble`: registers coherent Gaussian TBFs and projects the
  initial wavefunction using its overlap metric. Returns a `TrajectoryManager`.
  Uses the shared sampler for random placement and the shared nuclear assignment
  helper for coordinates, momenta and inverse masses.

QTAG-specific grid construction, state replication, widths/phases, and projection
helpers (`grid`, `gaussian`, `gbfs`, `wfc_centers`, `coeffs`) also live in
`ensembles.py`. `gaussian` here is a sampling function, not a separate module.
Import `make_gaussian_ensemble` from `libra_py.dyn.initialization` or directly
from `libra_py.dyn.initialization.ensembles`; the former `initialization.gaussian`
module has been removed.

Initial `wfc_a0` uses the GWP exponent `-a*(r-q)^2`; basis `alpha` uses
`-alpha*(r-q)^2/2`. Coherent coefficients must not be normalized independently
per TBF. These physical differences are kept explicit rather than routed through
independent-path amplitude normalization. The former standalone implementation lives in the test-only `qtag_reference` package.
