# Holstein reference densities

`holstein_reference.npz` contains `steps` (5), `x` (100), and `density`
(5, 2, 100). It was generated from the original
`libra_py.dynamics.qtag.compute.run_qtag` using the Python-exposed C++
implementation and the Holstein5 parameters in tutorial
`6_dynamics/5_qtag/1_basics/tutorial.ipynb` (case 1), on 2026-10-21.
The parameters are reproduced in the example `run.py`: 25 packets per state,
800 steps of 0.01, BATe, density-weighted iterative fit with beta=0.001.

Generation captured q, p, a, s, and coefficients passed to the legacy
`save_qtag_data` after each completed step. At steps 1, 200, 400, 600, and
800, densities were independently evaluated on `np.linspace(-4,4,100)` as

```
abs(sum_j c_j * (a_j/pi)**0.25
    * exp(-a_j*(x-q_j)**2/2 + 1j*(p_j*(x-q_j)+s_j)))**2
```

The sum uses the 25 packets belonging to each electronic state. This matches
the source tutorial's reconstruction (its phases are zero). No new array
propagation or integral code was used to generate the reference trajectory.
The source notebook and compiled extension are not needed to run the regression.
The test tolerance (2e-6 absolute, 2e-5 relative) allows accumulated platform
roundoff while easily rejecting the earlier direct-fit trajectory.
