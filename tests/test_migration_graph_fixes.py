"""Regression cases for reciprocal migration, image channels, and balance."""
import pickle
import numpy as np
import pytest
from scipy.sparse import coo_matrix, diags
from kmc.diffusion_average import DiffusionAverager, effective_transport, stationary_distribution
from kmc.particles import Interstitial
from kmc.potentials import MigrationBarrier
from kmc.rates import Migration
from kmc.system import System
from test_defect_events import make_system


def cube(n=2, disorder=0, seed=0):
    R=np.array([[x,y,z] for x in range(n) for y in range(n) for z in range(n)],float)
    R+=np.random.default_rng(seed).normal(0,disorder,R.shape)
    s=make_system(R,[0]*len(R),periodic=True,cell=np.diag([float(n)]*3),rtol=.15)
    s.animation_mode=s.bimolec=False
    s.creation_threshold=0.;s.time_limit=1.;s.particle_condition=True
    s.generation=[];s.monomolecular={};s.processes={'interstitial':[Migration(k={(0,0):3.})]}
    s.set_particles([Interstitial(0)])
    return s,DiffusionAverager(s,lambda p,s:None)


@pytest.mark.parametrize('seed',[0,1,2,4,5])
def test_irregular_migration_graph_is_reciprocal_and_has_no_drift(seed):
    s,a=cube(n=4,disorder=.04,seed=seed)
    edges={}
    for i in range(len(s.R)):
        cut,dx,dy,dz,_=s.neighbor_data(i,'migration')
        for j,x,y,z in zip(cut,dx,dy,dz):
            edges[(i,int(j),*np.round([x,y,z],12))]=True
    for i,j,x,y,z in edges:
        assert (j,i,*np.round([-x,-y,-z],12)) in edges
    group,=a.build_model().groups
    np.testing.assert_allclose(group.probability,np.ones(len(s.R))/len(s.R))
    np.testing.assert_allclose(group.drift,[0,0,0],atol=1e-12)


@pytest.mark.parametrize('n',[1,2,4])
def test_small_periodic_cells_keep_all_six_migration_images(n):
    s,a=cube(n)
    for i in range(len(s.R)):
        cut,dx,dy,dz,r=s.neighbor_data(i,'migration')
        assert len(cut)==6
        assert sorted(r)==pytest.approx([1]*6)
        assert set(zip(dx,dy,dz))=={(-1.,0.,0.),(1.,0.,0.),(0.,-1.,0.),(0.,1.,0.),(0.,0.,-1.),(0.,0.,1.)}
    group,=a.build_model().groups
    np.testing.assert_allclose(group.drift,[0,0,0],atol=1e-12)
    np.testing.assert_allclose(group.diffusion,np.diag([3.,3.,3.]),atol=1e-12)
    if n==2:
        cut,*_=s.neighbor_data(0,'first')
        assert len(cut)==3  # Each reaction partner is still counted once.
        np.testing.assert_allclose(np.asarray(s.edge_displacement(0,1)),
                                   -np.asarray(s.edge_displacement(1,0)))


def test_explicit_migration_samples_both_images_of_selected_destination():
    s,a=cube(2)
    p=s.particles[0];process=s.processes['interstitial'][0]
    s.set_rng(1234)
    samples={0:[],1:[]}
    for _ in range(2000):
        origin=p.position;old=p.Dz
        process.action(p,s,1-origin)
        samples[origin].append(p.Dz-old)
    for increments in samples.values():
        assert set(increments)=={-1.,1.}
        assert .45<np.mean(np.asarray(increments)>0)<.55


def test_default_site_profile_has_no_implicit_directional_weight():
    profile=MigrationBarrier(lambda x,L:np.zeros_like(x),temperature=300)
    assert profile.directional is False


@pytest.mark.parametrize('coupling',['destination','symmetric'])
@pytest.mark.parametrize('directional',[False,True])
def test_angular_profile_obeys_detailed_balance_on_reciprocal_cycle(coupling,directional):
    s=System();s.set_morph([0,1,2],[0,0,1],[0,0,0],[0]*3)
    U=np.array([0.,.03,.07]);s.migration_barrier_change=U;s.migration_kbt=.025852
    s.migration_barrier_axis=0;s.migration_barrier_angular_exponent=2
    s.migration_barrier_coupling=coupling;s.migration_barrier_directional=directional
    pairs=[(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    vectors=np.array([[1,0,0],[-1,1,0],[-1,0,0],[1,1,0],[1,-1,0],[-1,-1,0]],float)
    rates=[]
    for (i,j),d in zip(pairs,vectors):
        rates.append(s._migration_factor_for(i,np.array([j]),*[np.array([v]) for v in d],np.array([np.linalg.norm(d)]))[0])
    graph=coo_matrix((rates,([i for i,j in pairs],[j for i,j in pairs])),shape=(3,3)).tocsr()
    pi,_=stationary_distribution(graph-diags(np.asarray(graph.sum(axis=1)).ravel()))
    expected=np.exp(-U/s.migration_kbt);expected/=expected.sum()
    np.testing.assert_allclose(pi,expected)
    lookup=dict(zip(pairs,rates))
    for i,j in pairs:
        assert pi[i]*lookup[i,j]==pytest.approx(pi[j]*lookup[j,i])
    _,drift=effective_transport(graph,pi,vectors)
    np.testing.assert_allclose(drift,[0,0,0],atol=1e-12)


def test_image_specific_barrier_cache_survives_worker_pickling():
    s,a=cube(2)
    MigrationBarrier(lambda x,L:.04*x,temperature=300,directional=True).assign_to_system(s)
    s.build_neighbor_topology()
    worker=pickle.loads(pickle.dumps(s))
    for i in range(len(s.R)):
        cut,dx,dy,dz,r=worker.neighbor_data(i,'migration')
        expected=worker._migration_factor_for(i,cut,dx,dy,dz,r)
        np.testing.assert_allclose(worker.migration_rate_factors(i,cut),expected)


@pytest.mark.parametrize('directional',[False,True])
def test_grain_profile_on_irregular_reciprocal_graph_has_no_drift(directional):
    s,a=cube(n=4,disorder=.04,seed=0)
    MigrationBarrier(lambda x,L:-.03*np.exp(-(np.minimum(x,L-x)/.6)**2),
                     temperature=300,directional=directional).assign_to_system(s)
    s.build_neighbor_topology()
    group,=a.build_model().groups
    expected=np.exp(-s.migration_barrier_change/s.migration_kbt);expected/=expected.sum()
    np.testing.assert_allclose(group.probability,expected)
    np.testing.assert_allclose(group.drift,[0,0,0],atol=1e-12)


def test_half_cell_minimum_image_ties_are_reciprocal_for_scalar_and_array_queries():
    s,a=cube(2)
    for i in range(len(s.R)):
        forward=np.array(s.minimum_image_displacement(i))
        for j in range(len(s.R)):
            np.testing.assert_allclose(forward[:,j],-np.array(s.minimum_image_displacement(j,i)))
