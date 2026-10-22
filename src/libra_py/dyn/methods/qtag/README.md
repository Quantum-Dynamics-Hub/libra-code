# NumPy/PyTorch QTAG

This package implements quantum trajectory adaptive Gaussian (QTAG) dynamics
using NumPy arrays or PyTorch tensors. It ports the algebra in `src/dyn/qtag`
and the workflow in `libra_py.dynamics.qtag`, using `dyn.gwp.integrals` for
Gaussian overlaps and kinetic integrals. It does not import `liblibra_core`,
`util.libutil`, or SymPy. The old package remains available for C++ callers.

## Modules

| Module | Responsibilities |
| --- | --- |
| `core` | Wavefunctions, overlaps, kinetic integrals, BAT/LHA potential approximations, exact Gaussian coupling, quantum momentum and generalized electronic propagation |
| `initialize` | Cartesian or Gaussian-sampled bases, initial wavepacket parameters and projection coefficients |
| `compute` | Momentum fitting, basis updates/synchronization, re-expansion, observables, time propagation and grid reconstruction |
| `save` | Explicit NPZ/HDF5 persistence and readable observable tables |
| `plot` | Headless individual PNG plots of densities, energy, populations and basis trajectories |

All functions have docstrings describing their array contracts. Imports:

```python
from libra_py.dyn.methods.qtag import core, initialize, compute, save, plot
```

## Shapes, units and the width convention

Atomic units are used throughout. A QTAG basis function is

```text
G_j(r) = product_d (a_dj/pi)^(1/4)
         * exp(sum_d [-a_dj*(r_d-q_dj)^2/2
                      + i*p_dj*(r_d-q_dj) + i*s_dj]).
```

| Quantity | Shape | Meaning |
| --- | --- | --- |
| `q,p,alp,s` | `(ndof,ntraj)` | Centers (bohr), momenta, positive widths (bohr^-2), phases per DOF |
| `states` | list of `ntraj` integers | Electronic state label of each packet; ordering can be interleaved |
| `coeff` | `(ntraj,1)` or `(ntraj,)` | One complex expansion coefficient per labeled packet |
| `invM` | `(ndof,)` or `(ndof,1)` | Positive inverse mass per coordinate |
| `S,H` | `(ntraj,ntraj)` | Electronic-state-aware overlap and Hamiltonian |
| wavefunction evaluation points | `(ndof,npoints)` | One point per column |
| `core.qtag_psi` coefficients | `(nstates,ntraj)` | One coefficient row per electronic state |
| `core.qtag_psi` result | `(nstates,npoints)` | Complex wavefunction value, including phases |

**QTAG basis widths are twice GWP widths.** Core integral wrappers pass
`alp/2` to `dyn.gwp.integrals`. The original tutorial's `wfc_a0` describes the
initial wavepacket in the **GWP** convention `exp(-a0*(r-q0)^2)`; its basis
widths remain `wfc_a0*alp_scl` in the **QTAG** convention. `initialize.coeffs`
keeps this distinction explicitly. Phases are per coordinate, not one phase
for the entire vector.

## Model interface

Models operate on all centers simultaneously. A callable `model(q, params)`
returns `core.HamiltonianData(H,dH,d2H)` or a tuple with those entries:

- `H.shape == (ntraj,nstates,nstates)`.
- `dH.shape == (ntraj,ndof,nstates,nstates)`.
- Optional `d2H.shape == (ntraj,ndof,ndof,nstates,nstates)`.

The last two axes are electronic indices. Gradients are potential derivatives,
not forces. Matrices must be finite and Hermitian. LHA requires the full
Hessian, including mixed derivatives; the code raises if it is missing.
There is no hidden numerical Hessian approximation.

Existing `libra_py.dyn.models` objects exposing `evaluate(q, params)` can be
passed directly for BAT/BATe. Their result keys are `H_dia`, `dH_dia`, and
optionally `d2H_dia`. Choose a model backend matching the coordinates. A model
that returns C++ `MATRIX`/`CMATRIX` objects is not supported by this interface.
The tutorial example adds analytic Hessians to the existing Holstein5 and
Ferretti array models and demonstrates the callback contract.

## Potential approximations

`qtag_pot_approx_method` selects:

| Value | Diagonal electronic elements | Off-diagonal elements |
| --- | --- | --- |
| 0 | BAT | BAT |
| 1 | LHA | LHA |
| 2 | LHA | Exact Gaussian coupling |
| 3 | BAT | Exact Gaussian coupling |

`core.BAT`, `core.LHA`, `core.BATe`, `core.LHAe` and `core.qtag_potential`
return the **potential divided by the spatial overlap**. Matrix assembly
multiplies by the overlap; kinetic terms are added only for equal electronic
states. Direct `BAT`/`LHA` functions follow the C++ positional argument order,
using arrays instead of pointer lists. Direct LHA Hessians can also be flattened
as `(ndof**2,nstates,nstates)`, like the C++ lists.

For QTAG widths, the Gaussian product distribution has complex mean
`mu=(a1*q1+a2*q2+i*(p2-p1))/(a1+a2)` and diagonal covariance `1/(a1+a2)`.
BAT averages the first-order Taylor expansions at the two centers. LHA adds
one quarter of each endpoint Hessian contracted with
`E[(r-q_endpoint)(r-q_endpoint)^T]`, including mixed-coordinate terms. Thus
LHA is exact for quadratic potentials.

Methods 2/3 require `ex_cpl_A`, `ex_cpl_B`, `ex_cpl_C` in `model_params`.
Following C++, these represent
`sum_d A*exp(-B*(r_d-C)^2)`, with scalar `B>=0`, for every off-diagonal pair.
This is a **sum**, not a multidimensional product or a generic exact coupling.
The Holstein tutorial coupling fits this form; the Ferretti coupling does not.
For each coordinate the exact ratio is

```text
A / sqrt(1 + 2*B/(a1+a2))
  * exp[-B*(mu-C)^2/(1 + 2*B/(a1+a2))].
```

`core.qtag_overlap(active_states, spatial_overlap, nstates)` preserves the
C++ sparse `(ntraj*nstates,ntraj*nstates)` embedding for direct comparisons.
The runner uses the compact `ntraj`-square overlap instead.

## Initialization and propagation

`initialize.gbfs(params)` returns `(ntraj,q,p,alp,s,states)`. Important keys:

- `backend`: `numpy` (default) or `torch`; optional `device` and `dtype`.
- `init_placement`: Cartesian grid (0) or Gaussian sampling (1).
- `states`: electronic labels to populate with basis functions.
- `wfc_q0`, `wfc_p0`, `wfc_a0`, `wfc_s0`: one value per coordinate.
- `grid_dims`: per-coordinate grid counts, or counts whose product gives the
  sampled packet count per state.
- `alp_scl`: positive scaling of initial wavepacket widths for the basis.
- For grids, `0<rho_cut<1` gives density-relative bounds; `rho_cut<=0` uses
  `grid_min`/`grid_max`. For sampling, `rho_cut>0` is the standard deviation
  (the legacy meaning); `seed` controls a local reproducible NumPy RNG.

Sampling creates the same initial basis on each electronic state, with a clear
one-to-one synchronization map. Nonzero initial phases are preserved.
`initialize.wfc_centers` returns the initial GWP parameters.
`initialize.coeffs` projects it onto the basis and normalizes the projection
by default; set `normalize=False` to retain the legacy unnormalized projection.

```python
result = compute.run_qtag(q, p, alp, s, states, coeff, invM,
                          dyn_params, model, model_params)
q, p, alp, s, states, coeff = result  # optional legacy-style final-state unpacking
history = result.history
```

The step is electronic half-step, moving-basis update/re-expansion, electronic
half-step. Electronic evolution solves the Hermitian generalized eigenproblem
via the positive spectral subspace of S. `overlap_rcond` (default `1e-12`)
discards eigenvalues relative to the largest overlap eigenvalue; no explicit
matrix inverse is formed. Re-expansion solves `S_new*C_new = S_new,old*C_old`
in that retained subspace. Results are not renormalized after each step.

| `dyn_params` key | Default | Meaning |
| --- | --- | --- |
| `dt`, `nsteps` | 1, 1 | Positive time step and nonnegative number of complete steps |
| `qtag_pot_approx_method` | 0 | Potential approximation from the table above |
| `states` | Sorted basis labels | Population order used for trajectory synchronization |
| `target_states` | `states` | Population columns to save; does not restrict propagation |
| `q_update_method` | 1 | 0 freezes q; 1 applies `q += dt*invM*momentum` |
| `p_update_method` | 1 | 0 freezes p; 1 replaces p with the quantum momentum |
| `a_update_method` | 0 | 1 applies `a -= 2*dt*invM*a*gradient(momentum)` |
| `s_update_method` | 0 | Only 0 is supported, as no phase evolution is implemented in legacy QTAG |
| `q_sync_method`, `p_sync_method` | 1, 1 | Copy parameters from the most populated surface when below threshold |
| `a_sync_method`, `s_sync_method` | 0, 0 | Same synchronization control for widths/phases |
| `decpl_den` | 0.25 | Population threshold for independent surface motion |
| `mom_calc_type` | 1 | 0 raw momentum/zero slopes; 1 fitted; 2 raw momentum/analytic slopes |
| `d_weight` | 1 | Density-weighted fit (0 disables weighting) |
| `fit_rcond` | `1e-12` | Spectral cutoff in the direct linear fit |
| `momentum_fit` | see below | `direct` or `legacy` iterative fit |
| `linfit_beta` | `0.1` in legacy mode | Iterative coefficient-change tolerance |
| `node_tol` | `1e-14` | Relative wavefunction-amplitude threshold below which raw ratios are zeroed |
| `overlap_rcond` | `1e-12` | Relative metric cutoff for electronic evolution and projection |

Linear fitting defaults to a direct weighted least-squares solve. Supplying
`linfit_beta` selects the legacy iterative fit unless `momentum_fit` is explicitly
set. The tutorial uses `momentum_fit="legacy", linfit_beta=1e-3`. This reproduces
the C++ Gauss-Seidel iterations and stopping criterion using arrays only. The
tolerance affects basis motion and subsequent densities; a direct solve or a
tighter tolerance does not reproduce the tutorial trajectory. Stable population ordering handles
ties. The most populated nonzero surface is always available as a synchronization
reference, even when every population is below `decpl_den`. Synchronization
requires equal packet counts on the affected surfaces. Width updates use the
legacy Euler expression; nonpositive results raise and suggest reducing dt.
Unsupported phase updates raise rather than being silently ignored.

PyTorch graphs are preserved through calculations at fixed spectral rank and
fixed synchronization choices. Threshold decisions, state ordering and rank
selection are discrete and are not differentiable. Degenerate eigenspaces can
make eigendecomposition derivatives undefined. Use `torch.no_grad()` for long
runs without parameter-gradient requirements; histories otherwise retain graphs.
Initialization's grid/random sampling is construction, not a differentiable
operation with respect to initialization-dictionary values.

## Results, saving and plots

Histories contain `nsteps+1` snapshots: the initial state at t=0 and each fully
completed step. Coordinates and observables at a snapshot refer to the same
time (correcting the old half-step observable timing).

- `time`, `Etot`, `dEtot`, `norm`: `(ntimes,)`.
- `pops`: `(ntimes,len(target_states))`.
- `q,p,a,s`: `(ntimes,ndof,ntraj)`.
- `coeffs`: `(ntimes,ntraj)` complex.
- `states`, `target_states`: label arrays.

`Etot` is `C^H H C` and `norm` is `C^H S C`, matching the legacy expectation
convention; energy is not divided by the norm. `dEtot` is its difference from
the initial value. Finite moving-basis projection may lose norm and drift in
energy; both are reported so users can assess convergence.

```python
save.save_qtag_data(result, "out/trajectory.npz")
save.save_text(result, "out")
loaded = save.load_qtag_data("out/trajectory.npz", backend="numpy")
wave = compute.wf_calc_nD(loaded, {
    "xmin": [-4.], "xmax": [4.], "npoints": [100],
    "snaps": [0, 10], "which_states": [0, 1],
})
plot.wf_plot(wave, "out")
plot.energy_and_pops(loaded, "out")
plot.trajectories(loaded, "out")
```

NPZ requires only NumPy; `.h5`/`.hdf5` additionally require h5py. Matplotlib is
imported only when plotting; its Agg canvas needs no GUI. Each plot is an
individual PNG. Wavefunction reconstruction uses saved phases and actual state
labels, not assumed packet block sizes. `save`/`plot` detach tensors to CPU
only at their output boundary. Data stays in its original backend in memory.

The old saver setup helpers (`init_qtag_savers`, `save_qtag_hdf5_*`, output-level
flags) are replaced by explicit archive/table writes. Old plot dictionaries
that inferred data from flattened text files are replaced by the result/history
and reconstructed-grid objects. There is no hidden file I/O in `run_qtag`.

## C++ audit and verification

The port corrected these defects in the corresponding C++ code:

1. LHA used the bra Hessian twice instead of using both endpoints.
2. Mixed LHA moments used the other coordinate's width denominator and paired
   the mean shifts with the wrong coordinate centers.
3. Quantum momentum used the evaluation packet's normalization width instead
   of the expansion packet's width and wrote `(itraj,idof)` into an
   `(ndof,ntraj)` matrix. Real output is the logarithmic amplitude gradient;
   the imaginary output is the quantum momentum.
4. Hamiltonian assembly grouped contiguous runs rather than actual state
   labels, and left stale output overlap entries when labels were reordered.
5. The C++ model callback copied only `ndof` entries of an `ndof**2` Hessian.

See the tests under `unittests/test_libra_py/test_dyn/test_methods/test_qtag`.
They compare against the rebuilt extension and independently check Gaussian
quadrature, mixed-coordinate LHA moments, log-gradient finite differences,
metric-unitary propagation, state permutation invariance, NumPy/PyTorch
trajectories, scalar autograd, file round trips and headless plots.

The pure-Python tutorial lives in
`examples/libra_py/dyn/methods/qtag/1_basics`. Its README includes commands,
output descriptions, and measured limitations of the notebook parameters.
