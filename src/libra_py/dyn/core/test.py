import numpy as np
import storage as st

storage = st.TensorStorage(
    backend=np,
    ntraj=100,
    ndof=300,
    nstates=4,
    ntbf_initial=2,
)

storage.q[5, 0] += 0.1

new_tbf = storage.spawn(traj=5)

storage.q[5, new_tbf] = storage.q[5, 0]
storage.ampl_adi[5, new_tbf, 2] = 1.0

storage.prune(10, 1)
storage.compact()
