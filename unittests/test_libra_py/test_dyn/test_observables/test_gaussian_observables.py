"""Coherent observables survive generic savers and retain cross-packet terms."""
import numpy as np
import pytest
from libra_py.dyn.initialization.ensembles import make_gaussian_ensemble
from libra_py.dyn.methods.qtag import QTAGEngine
from libra_py.dyn.models import FerrettiModel
from libra_py.dyn.observables.gaussian import reconstruct_density
from libra_py.dyn.observables import plotting
from libra_py.dyn.savers import FaultTolerantSaver, load_saved_steps, stack_saved_observables
from libra_py.dyn.backends.arrays import to_numpy


@pytest.mark.parametrize('backend',['numpy','torch'])
def test_save_and_reconstruct(tmp_path,backend):
    if backend=='torch':
        pytest.importorskip('torch')
    manager=make_gaussian_ensemble(dict(backend=backend,wfc_q0=[0.,1.],grid_dims=[3,3],
                                        wfc_s0=[.3,.7],alp_scl=[4.,4.],rho_cut=.1),[.5,.7])
    engine=QTAGEngine(manager,FerrettiModel(backend=backend),dict(dt=.001))
    result=engine.run(2,savers=FaultTolerantSaver(tmp_path/'snapshots'))
    h=result.histories[0]
    loaded=stack_saved_observables(load_saved_steps(tmp_path/'snapshots'))
    for key in h:
        np.testing.assert_array_equal(loaded[key],to_numpy(h[key]))
    # C†SC includes interference; sum |C_i|² is generally not the norm.
    assert not np.isclose(np.sum(abs(to_numpy(h['coeffs'][0]))**2),float(h['norm'][0]),rtol=1e-3)
    grid=dict(xmin=[-2.,-1.],xmax=[2.,3.],npoints=[11,11],snaps=[0,2])
    original=reconstruct_density(h,**grid)
    restored=reconstruct_density(loaded,**grid)
    np.testing.assert_allclose(to_numpy(original['density']),restored['density'],atol=1e-12)
    paths=plotting.energy_and_pops(h,tmp_path)+plotting.trajectories(h,tmp_path)+plotting.wf_plot(original,tmp_path)
    assert len(paths)==16
    assert all(path.exists() and path.stat().st_size>1000 for path in paths)
    continued=engine.run(1,keep_history=False)
    assert continued.histories=={}
    assert continued.storage.timestep==3
