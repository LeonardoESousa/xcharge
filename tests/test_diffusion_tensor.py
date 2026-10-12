"""Transport tensors checked against analytic periodic walks and log output."""
import numpy as np
import pytest
from scipy.sparse import coo_matrix
from kmc.diffusion_average import DiffusionAverager, effective_diffusion, effective_transport, stationary_distribution
from kmc.particles import Interstitial, Vacancy, diffusion_log_values, DIFFUSION_LOG_COLUMNS
from kmc.rates import Migration, Formation
from test_defect_events import make_system
from test_diffusion_average import mixed


def ring(bond_rates):
    n = len(bond_rates)
    rows, cols, rates, vectors = [], [], [], []
    for i in range(n):
        for j, k, dx in [((i+1)%n, bond_rates[i], 1),
                         ((i-1)%n, bond_rates[(i-1)%n], -1)]:
            rows.append(i); cols.append(j); rates.append(k); vectors.append([dx,0,0])
    # Align vectors with the sorted sparse edge ordering.
    ordered = sorted(zip(rows,cols,rates,vectors))
    graph = coo_matrix(([x[2] for x in ordered],
                       ([x[0] for x in ordered],[x[1] for x in ordered])),shape=(n,n)).tocsr()
    return graph, np.asarray([x[3] for x in ordered],dtype=float)


def test_homogeneous_ring_and_heterogeneous_harmonic_mean():
    for rates in [[3]*8, [1,9]*4, [1]*4+[9]*4]:
        graph, vectors = ring(rates)
        pi = np.ones(len(rates))/len(rates)
        tensor = effective_diffusion(graph, pi, vectors)
        expected = np.zeros((3,3))
        expected[0,0] = len(rates)/np.sum(1/np.asarray(rates))
        np.testing.assert_allclose(tensor,expected,atol=1e-12)
    # A local-rate arithmetic average would incorrectly give 5 instead of 1.8.
    assert tensor[0,0] == pytest.approx(1.8)


def test_directional_tensor_and_offdiagonal_covariance():
    graph, vectors = ring([2]*6)
    vectors[:,1] = vectors[:,0]
    tensor = effective_diffusion(graph,np.ones(6)/6,vectors)
    np.testing.assert_allclose(tensor, [[2,2,0],[2,2,0],[0,0,0]])


def test_destination_barrier_profile_on_periodic_ring():
    s,a = mixed(n=8,migration=3)
    s.periodic = True
    s.set_cell(np.diag([8.,20.,20.]))
    s.migration_barrier_change = np.r_[np.zeros(4),np.full(4,.04)]
    s.migration_kbt = .025
    s._cached_translations = None
    s.build_neighbor_topology()
    model = a.build_model()
    g = np.exp(-s.migration_barrier_change/s.migration_kbt)
    pi = g/g.sum()
    conductances = pi*3*np.roll(g,-1)
    expected = 8**2/np.sum(1/conductances)
    for group in model.groups:
        np.testing.assert_allclose(group.probability,pi)
        np.testing.assert_allclose(group.diffusion,np.diag([expected,0,0]),atol=1e-12)


def cube():
    coordinates = [[x,y,z] for x in range(4) for y in range(4) for z in range(4)]
    s = make_system(coordinates,[0]*len(coordinates),periodic=True,cell=np.diag([4.,4.,4.]))
    s.animation_mode = s.bimolec = False
    s.time_limit = 10.
    s.creation_threshold = 0.
    s.particle_condition = True
    s.generation = []
    s.monomolecular = {}
    s.processes = {'interstitial':[Migration(k={(0,0):3})],
                   'vacancy':[Migration(k={(0,0):1})]}
    s.set_particles([Interstitial(0),Vacancy(len(coordinates)-1)])
    a = DiffusionAverager(s,lambda p,s: None)
    s.diffusion_averager = a
    return s,a


def test_three_periodic_directions_and_accumulated_tensor():
    s,a = cube()
    # At z=0, z=3 is a wrapped neighbor with displacement -1.
    cut,*_ = s.neighbor_data(0,'first')
    assert 3 in cut
    assert 1 in cut
    # x=3 wraps to x=0 with displacement -1 along x.
    assert s.edge_displacement(0,48) == (-1.,0.,0.)
    assert s.edge_displacement(0,3) == (0.,0.,-1.)
    model = a.build_model()
    for group in model.groups:
        k = 3 if group.species=='interstitial' else 1
        np.testing.assert_allclose(group.diffusion,np.diag([k,k,k]),atol=1e-12)
    assert not a.advance(lambda _:pytest.fail('unexpected direct step'))
    for p in s.particles:
        values = dict(zip(DIFFUSION_LOG_COLUMNS,diffusion_log_values(p,s)))
        k = 3 if p.species=='interstitial' else 1
        assert values['AveragedDiffusionTime']==10
        assert values['D_est_xx']==pytest.approx(k)
        assert values['D_est_zz']==pytest.approx(k)
        assert values['MSD_est_x']==pytest.approx(20*k)
        assert values['MSD_est_z']==pytest.approx(20*k)


def test_accumulation_over_changed_environment_and_particle_reset():
    s,a = cube()
    first = a.build_model()
    a._accumulate_diffusion(first,2)
    # A retained obstacle changes the graph after a slow event.
    s.set_material(1,999)
    second = a.build_model()
    a._accumulate_diffusion(second,3)
    assert second is not first
    for old,new in zip(first.groups,second.groups):
        for p in new.particles:
            np.testing.assert_allclose(p.averaged_diffusion_integral,2*old.diffusion+3*new.diffusion)
            assert p.averaged_diffusion_time==5
    newborn = Interstitial(0)
    assert newborn.averaged_diffusion_time==0
    assert not np.any(newborn.averaged_diffusion_integral)


def test_biased_ring_has_centered_diffusion_and_drift():
    graph,vectors = ring([1]*6)
    edges = graph.tocoo()
    rates = np.where(vectors[:,0]>0,2.,1.)
    graph = coo_matrix((rates,(edges.row,edges.col)),shape=graph.shape).tocsr()
    tensor, drift = effective_transport(graph,np.ones(6)/6,vectors)
    np.testing.assert_allclose(tensor,np.diag([1.5,0,0]))
    np.testing.assert_allclose(drift,[1,0,0])


def test_confined_walk_has_zero_long_time_tensor():
    s,a = mixed()
    for group in a.build_model().groups:
        np.testing.assert_allclose(group.diffusion,np.zeros((3,3)),atol=1e-12)


def test_large_reversible_iterative_corrector():
    # Exercises the solver used to avoid LU fill-in in large morphologies.
    graph,vectors = ring([1,9]*3000)
    tensor = effective_diffusion(graph,np.ones(6000)/6000,vectors)
    np.testing.assert_allclose(tensor,np.diag([1.8,0,0]),atol=1e-7)


def test_nonuniform_biased_cycle_against_tilted_generator():
    # Independent reference: first two derivatives of the displacement CGF.
    from scipy.sparse import diags
    graph,vectors = ring([1]*8)
    edges = graph.tocoo()
    rates = np.where(vectors[:,0]>0,1+edges.row,0.7)
    graph = coo_matrix((rates,(edges.row,edges.col)),shape=graph.shape).tocsr()
    exits = np.asarray(graph.sum(axis=1)).ravel()
    pi,_ = stationary_distribution(graph-diags(exits))
    tensor,drift = effective_transport(graph,pi,vectors)
    assert abs(drift[0])>0.1
    def eigenvalue(tilt):
        matrix = coo_matrix((rates*np.exp(tilt*vectors[:,0]),
                             (edges.row,edges.col)),shape=graph.shape).toarray()
        matrix -= np.diag(exits)
        return np.linalg.eigvals(matrix).real.max()
    h = 1e-4
    plus,zero,minus = eigenvalue(h),eigenvalue(0),eigenvalue(-h)
    assert drift[0] == pytest.approx((plus-minus)/(2*h),rel=1e-6)
    assert tensor[0,0] == pytest.approx((plus-2*zero+minus)/(2*h*h),rel=1e-5)


def test_drift_accumulation_log_separates_msd_and_variance():
    s,a = cube()
    model = a.build_model()
    # Two slow intervals with reversed macroscopic drift.
    for group in model.groups:
        group.drift = np.array([2.,-1.,0.])
    a._accumulate_diffusion(model,2)
    for group in model.groups:
        group.drift = np.array([-1.,2.,0.])
    a._accumulate_diffusion(model,3)
    for p in s.particles:
        values = dict(zip(DIFFUSION_LOG_COLUMNS,diffusion_log_values(p,s)))
        k = 3 if p.species=='interstitial' else 1
        assert values['MeanDelta_est_x']==pytest.approx(1)
        assert values['MeanDelta_est_y']==pytest.approx(4)
        assert values['Drift_est_x']==pytest.approx(.2)
        assert values['Drift_est_y']==pytest.approx(.8)
        assert values['Var_est_x']==pytest.approx(10*k)
        assert values['Var_est_y']==pytest.approx(10*k)
        assert values['MSD_est_x']==pytest.approx(10*k+1)
        assert values['MSD_est_y']==pytest.approx(10*k+16)
        assert values['Cov_est_xy']==pytest.approx(0)
    newborn=Vacancy(0)
    assert not np.any(newborn.averaged_drift_integral)


def test_biased_lattice_integration_continues_without_zero_drift_error():
    # Three contiguous material regions, with directional inter-region rates.
    s,a = mixed(n=12,migration=2)
    s.periodic=True
    s.set_cell(np.diag([12.,20.,20.]))
    s._cached_translations=None
    s.mats[:]=np.repeat([0,1,2],4)
    s.uniq=np.unique(s.mats)
    s.build_neighbor_topology()
    rates={(i,j):(4. if (j-i)%3==1 else 1.) for i in range(3) for j in range(3)}
    s.processes={'interstitial':[Migration(k=rates),Formation(k=1)],
                 'vacancy':[Migration(k=rates)]}
    s.diffusion_averager=a
    s.set_rng(324)
    model=a.build_model()
    assert any(np.linalg.norm(group.drift)>1e-3 for group in model.groups)
    assert a.advance(lambda _:pytest.fail('unexpected direct event'))
    assert s.IT==1
    for particle in s.dead:
        values=dict(zip(DIFFUSION_LOG_COLUMNS,diffusion_log_values(particle,s)))
        assert np.isfinite(values['D_est_xx'])
        assert abs(values['Drift_est_x'])>1e-3
        assert values['MSD_est_x']==pytest.approx(
            values['Var_est_x']+values['MeanDelta_est_x']**2)
