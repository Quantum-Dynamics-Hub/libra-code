"""Shared Gaussian kernels match the untouched QTAG kernels for unequal widths."""
import numpy as np
import pytest
from libra_py.dyn.gwp import packets, potentials
from qtag_reference import core as reference
from libra_py.dyn.backends.arrays import to_numpy


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_gaussian_kernels(backend):
    rng=np.random.default_rng(41)
    q,p,s=[rng.normal(size=(2,4)) for _ in range(3)]
    a=rng.uniform(.7,2.1,size=q.shape)
    points=rng.normal(size=(2,9))
    c=rng.normal(size=4)+1j*rng.normal(size=4)
    if backend=='torch':
        torch=pytest.importorskip('torch')
        q,p,a,s,points=[torch.as_tensor(x,dtype=torch.float64) for x in (q,p,a,s,points)]
        c=torch.as_tensor(c)
    functions=[(packets.spatial_overlap,reference.qtag_overlap_elementary,(q,p,a,s)),
               (packets.spatial_kinetic,reference.qtag_kinetic_elementary,(q,p,a,s,[.2,.3])),
               (packets.wavefunction,reference.qtag_psi,(points,q,p,a,s,c[None])),
               (packets.momentum_fields,reference.momentum_fields,(q,p,a,s,c))]
    for new,old,args in functions:
        actual,expected=new(*args),old(*args)
        if not isinstance(actual,tuple): actual,expected=(actual,),(expected,)
        for x,y in zip(actual,expected):
            np.testing.assert_allclose(to_numpy(x),to_numpy(y),atol=1e-12)
    actual=potentials.gaussian_coupling_ratio(q,p,a,q,p,a,(.4,.3,.8))
    expected=reference._exact_coupling(q,p,a,q,p,a,(.4,.3,.8))
    np.testing.assert_allclose(to_numpy(actual),to_numpy(expected),atol=1e-12)
