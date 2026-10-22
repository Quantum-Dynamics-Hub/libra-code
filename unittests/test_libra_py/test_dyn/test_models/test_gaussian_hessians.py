"""Analytic model Hessians support LHA without example-local model copies."""
import numpy as np
import pytest
from libra_py.dyn.models import Holstein5Model, FerrettiModel
from libra_py.dyn.backends.arrays import to_numpy


@pytest.mark.parametrize('backend',['numpy','torch'])
@pytest.mark.parametrize('case',['holstein','ferretti'])
def test_hessian_finite_difference(backend,case):
    if case=='holstein':
        q=np.array([[-.3,0.,1.2]])
        model=Holstein5Model(backend=backend,params=dict(E_n=[0.,.1],x_n=[0.,1.],k_n=[1.,2.],
            V=[[0.,.2],[.2,0.]],alpha=[[0.,.7],[.7,0.]],x_nm=[[0.,.4],[.4,0.]]))
    else:
        q=np.array([[2.5,3.,3.2],[-.3,0.,.6]])
        model=FerrettiModel(backend=backend)
    if backend=='torch':
        q=pytest.importorskip('torch').as_tensor(q)
    hess=to_numpy(model.evaluate(q)['d2H_dia'])
    eps=1e-5
    for d in range(q.shape[0]):
        dq=q*0; dq[d]=eps
        fd=(to_numpy(model.evaluate(q+dq)['dH_dia'])-to_numpy(model.evaluate(q-dq)['dH_dia']))/(2*eps)
        np.testing.assert_allclose(hess[:,:,d],fd,atol=2e-9,rtol=2e-7)
    np.testing.assert_allclose(hess,hess.swapaxes(1,2),atol=1e-14)
