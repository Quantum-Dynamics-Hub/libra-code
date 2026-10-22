"""Shared Gaussian initialization preserves the reference sampling/projection."""
import numpy as np
import pytest
from libra_py.dyn.initialization import ensembles
from qtag_reference import initialize as reference
from libra_py.dyn.backends.arrays import to_numpy
from libra_py.dyn.gwp.storage import read_packets


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
@pytest.mark.parametrize('placement', [0, 1])
def test_initialization_reference(backend, placement):
    if backend == 'torch':
        pytest.importorskip('torch')
    pars = dict(backend=backend, wfc_q0=[.2,-.1], wfc_p0=[.4,.6], wfc_a0=[1.2,.8],
                wfc_s0=[.3,-.2], grid_dims=[3,4], alp_scl=[4.,5.], states=[2,0],
                init_placement=placement, seed=13, rho_cut=.05 if placement==0 else .8)
    expected = reference.gbfs(pars)
    actual = ensembles.gbfs(pars)
    for a, b in zip(actual[1:5], expected[1:5]):
        np.testing.assert_array_equal(to_numpy(a),to_numpy(b))
    manager = ensembles.make_gaussian_ensemble(pars,[.5,.7])
    slots, packet, labels, c, mass = read_packets(manager.storage,manager.trajectories[0])
    q0,p0,a0,s0 = reference.wfc_centers(pars)
    ref = reference.coeffs(q0,p0,a0,s0,[2],*expected[1:5],expected[-1])
    np.testing.assert_allclose(to_numpy(c),to_numpy(ref[:,0]),atol=1e-12)
    assert labels == expected[-1]
    assert len(manager.tbfs) == 24
