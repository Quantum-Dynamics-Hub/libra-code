# State-tracking and projector examples

These examples use three electronic states and run directly from this
directory. They assume the active environment already provides `libra_py` and
do not modify `sys.path`.

```bash
cd examples/libra_py/dyn/transformations
python 01_state_tracking.py
python 02_projectors.py
```

`01_state_tracking.py` demonstrates the C++ permutation convention,
energy-aware cost matrices, deterministic and stochastic assignment, and
active-state remapping.

`02_projectors.py` demonstrates permutation-plus-phase projectors and
local-diabatization projectors imported through `projectors.py`.

The relevant state-tracking selectors are:

- `1`: greedy maximum-overlap tracking;
- `2`: Munkres–Kuhn/Hungarian tracking;
- `21`: alternative Hungarian implementation;
- `3`, `32`, `33`: stochastic tracking variants;
- `4`: force-based assignment.

Projection-update selectors additionally include `-1` for local
diabatization, `0` for no update, and `5`/`6` for the SVD-based LD variants.
