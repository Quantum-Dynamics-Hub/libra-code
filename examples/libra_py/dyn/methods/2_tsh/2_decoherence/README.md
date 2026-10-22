# Tully-1 hopping and decoherence

Run `python run.py --nsteps 1200 --ntraj 12 --save-stride 10` to compute
and plot all nine combinations of FSSH/GFSH/FSSH2 with coherent, ID-A,
and SDM dynamics. All simulation and plotting functions are in `run.py`.
Use `--compute-only` to save snapshots without plotting, or `--plot-only`
to regenerate figures from existing snapshots. `--steps` is an alias for
`--nsteps`; a quick check is `python run.py --steps 2 --ntraj 2 --save-stride 1`.
Set PYTHONPATH to the repository's absolute src path before running.

All hopping cases select `method="tsh"`; `tsh_method` selects FSSH (0),
GFSH (1), or FSSH2 (7). Each case
uses the same seeded nuclear ensemble and hopping RNG.

NPZ snapshots go to `output/tsh_decoherence_comparison/<case>`; PNGs go
to `plots/tsh_decoherence_comparison`. Plots compare electronic and hopping
populations, energy drift, and active-state histories. All units are atomic.
