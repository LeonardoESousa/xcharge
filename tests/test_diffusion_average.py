"""Analytic and executed-event checks for spatial diffusion averaging."""
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from kmc.diffusion_average import DiffusionAverager, stationary_distribution
from kmc.particles import Interstitial, Vacancy
from kmc.rates import Formation, Migration, DissociationFP
from test_defect_events import make_system


def mixed(n=20, migration=1e8, formation=2.0, both=False):
    s = make_system([[i, 0, 0] for i in range(n)], [0]*n)
    s.animation_mode = s.bimolec = False
    s.time_limit = np.inf
    s.creation_threshold = 0.0
    s.particle_condition = True
    s.generation = []
    s.monomolecular = {}
    s.processes = {'interstitial': [Migration(k={(0,0): migration}), Formation(k=formation)],
                   'vacancy': [Migration(k={(0,0): migration})]}
    if both:
        s.processes['vacancy'].append(Formation(k=formation))
    s.set_particles([Interstitial(0), Vacancy(n-1)])
    return s, DiffusionAverager(s, lambda p,s: None)


def test_stationary_symmetric_and_asymmetric_generators():
    for q, expected in [([[-2,2],[2,-2]], [.5,.5]),
                        ([[-2,2],[1,-1]], [1/3,2/3]),
                        ([[-1,1,0],[0,-1,1],[1,0,-1]], [1/3]*3)]:
        pi, residual = stationary_distribution(csr_matrix(q, dtype=float))
        np.testing.assert_allclose(pi, expected)
        assert residual < 1e-10


def test_single_pair_probability_and_two_initiator_channels():
    for both in [False, True]:
        s,a = mixed(both=both)
        model = a.build_model()
        # Open chain: 2(N-1) ordered neighboring configurations / N^2.
        expected = 2*19/20**2
        assert model.total_rate == pytest.approx(2*expected*(1+both))
        np.testing.assert_allclose(model.groups[0].probability, np.ones(20)/20)
        assert a.build_model() is model


def test_localized_barrier_and_local_formation_weights():
    s,a = mixed(n=12)
    s.migration_barrier_change = np.r_[np.zeros(6), np.full(6, .04)]
    s.migration_kbt = .025
    s.mats[6:] = 1
    migration = Migration(k={(i,j): 1e8 for i in [0,1] for j in [0,1]})
    formation = Formation(k={(i,j): (5 if i == 1 else 1) for i in [0,1] for j in [0,1]})
    s.processes = {'interstitial':[migration,formation], 'vacancy':[migration]}
    model = a.build_model()
    expected_pi = np.exp(-s.migration_barrier_change/s.migration_kbt)
    expected_pi /= expected_pi.sum()
    for group in model.groups:
        np.testing.assert_allclose(group.probability, expected_pi)
    channel = model.channels[0]
    expected = expected_pi[channel.origins]*expected_pi[channel.destinations]*np.where(channel.origins>=6,5,1)
    np.testing.assert_allclose(np.diff(np.r_[0,channel.cumulative_rates]), expected)
    rng = np.random.default_rng(321)
    selected = np.array([channel.select(rng)[0] for _ in range(10000)])
    probability = expected[channel.origins>=6].sum()/expected.sum()
    assert abs(np.mean(selected>=6)-probability) < .02


def test_executed_formation_clock_product_and_no_hops():
    times = []
    for seed in range(400):
        s,a = mixed()
        s.set_rng(seed)
        rate = a.build_model().total_rate
        assert a.advance(lambda _: pytest.fail('unexpected direct step'))
        assert s.IT == 1
        assert a.stats['averaged_events'] == 1
        pair, = s.particles
        assert pair.species == 'frenkelpair0'
        assert abs(pair.position-pair.ghost_site) == 1
        assert s.mats[pair.ghost_site] == 999
        assert len(s.dead) == 2
        times.append(s.time*rate)
        DissociationFP(k={0:1}).action(pair,s,pair.position)
        assert a.build_model().total_rate == pytest.approx(rate)
        assert a.stats['model_builds'] == 2
    assert abs(np.mean(times)-1) < .15
    assert abs(np.mean(np.asarray(times)<=1)-(1-np.exp(-1))) < .08


def test_time_boundary_and_stationary_vacancy():
    s,a = mixed()
    s.processes['vacancy'] = []
    model = a.build_model()
    vacancy = next(g for g in model.groups if g.species=='vacancy')
    assert vacancy.probability[-1] == 1
    assert model.total_rate == pytest.approx(2/20)
    s.time_limit = 1e-12
    assert not a.advance(lambda _: pytest.fail('direct step'))
    assert s.time == s.time_limit
    assert s.IT == 0
    assert len(s.particles) == 2


def test_unsupported_spatial_process_rejected():
    s,a = mixed()
    s.processes['interstitial'].append(object())
    with pytest.raises(ValueError, match='only Migration and Formation'):
        a.build_model()


def test_against_exact_two_particle_exclusion_chain():
    """Exact absorption mean retains every hop and forbids shared sites.

    The independent approximation permits virtual shared-site occupancy:
    its dilute error persists even as migration becomes infinitely fast.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve
    n = 30
    states = [(i,j) for i in range(n) for j in range(n) if i != j]
    lookup = {state:index for index,state in enumerate(states)}
    for hopping in [1e3, 1e6]:
        rows,cols,values = [],[],[]
        for index,(i,j) in enumerate(states):
            total = 2.0 if abs(i-j)==1 else 0.0
            for species in [0,1]:
                origin = (i,j)[species]
                for destination in [origin-1,origin+1]:
                    if not 0<=destination<n or destination==(j,i)[species]:
                        continue
                    target = (destination,j) if species==0 else (i,destination)
                    rows.append(index); cols.append(lookup[target]); values.append(-hopping)
                    total += hopping
            rows.append(index); cols.append(index); values.append(total)
        minus_q = coo_matrix((values,(rows,cols)), shape=(len(states),len(states))).tocsr()
        exact_mean = float(spsolve(minus_q, np.ones(len(states))).mean())
        s,a = mixed(n=n,migration=hopping)
        reduced_mean = 1/a.build_model().total_rate
        relative_error = reduced_mean/exact_mean-1
        assert abs(relative_error) < .04
        if hopping == 1e6:
            assert relative_error == pytest.approx(1/(n-1), abs=1e-4)


def test_cli_averaging_enabled(tmp_path):
    import os
    import subprocess
    import sys
    from pathlib import Path
    from test_cli import CONFIG
    configuration = CONFIG.replace('num_sites=4', 'num_sites=64').replace(
        'vector=[1.0, 0.0, 0.0]', 'vector=[1.0, 1.0, 1.0]') + '''
from kmc.rates import Migration, Formation
from kmc.rates import Annihilation

del exciton
diffusion_average = True
periodic = True
rounds = 2
pair_i = morphology.Create_Particles("interstitial", 1, morphology.randomized, mat=[0])
pair_v = morphology.Create_Particles("vacancy", 1, morphology.randomized, mat=[0])
migration = Migration(k={(0,0):1e12})
formation = Formation(k=2)
processes = {"interstitial":[migration,formation], "vacancy":[migration], "frenkelpair0":[]}
monomolecular = {"interstitial":[], "vacancy":[], "frenkelpair0":[Annihilation(k={0:1})]}
'''
    path = tmp_path/'input.py'
    path.write_text(configuration)
    env = dict(os.environ)
    env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1])+os.pathsep+env.get('PYTHONPATH','')
    result = subprocess.run([sys.executable,'-m','kmc',str(path)],cwd=tmp_path,
                            env=env,capture_output=True,text=True,timeout=20)
    assert result.returncode == 0, result.stdout+result.stderr
    output = (tmp_path/'Simulation_smoke.txt').read_text()
    assert 'formation' in output
    assert 'annihilation' in output
    assert 'nan' in output
    import csv
    import io
    from kmc.particles import DIFFUSION_LOG_COLUMNS
    rows = [row for row in csv.DictReader(io.StringIO('\n'.join(
        line for line in output.splitlines() if not line.startswith('#'))))
        if row['Time'] != 'END']
    assert all(column in rows[0] for column in DIFFUSION_LOG_COLUMNS)
    assert all(None not in row for row in rows)
    for row in rows:
        if row['Type'] in ('interstitial','vacancy'):
            duration = float(row['AveragedDiffusionTime'])
            assert duration == pytest.approx(float(row['Time']))
            assert float(row['D_est_xx']) == pytest.approx(1e12)
            assert float(row['MSD_est_x']) == pytest.approx(2e12*duration)
            assert float(row['D_est_yy']) == pytest.approx(1e12)
            assert float(row['D_est_zz']) == pytest.approx(1e12)
            assert float(row['MSD_est_z']) == pytest.approx(2e12*duration)
        else:
            assert np.isnan(float(row['D_est_xx']))
    # Appending a file produced with the previous schema must fail clearly.
    output_path = tmp_path/'Simulation_smoke.txt'
    output_path.write_text('Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,FinalX,FinalY,FinalZ,CausaMortis,Status\n')
    path.write_text(configuration+'\nappend_output = True\n')
    result = subprocess.run([sys.executable,'-m','kmc',str(path)],cwd=tmp_path,
                            env=env,capture_output=True,text=True,timeout=20)
    assert result.returncode != 0
    assert 'different column schema' in result.stderr


def test_biased_grain_regions_cli_with_multiprocessing(tmp_path):
    import csv
    import os
    import subprocess
    import sys
    from pathlib import Path
    from test_cli import CONFIG
    config = CONFIG + '''
from kmc.rates import Migration, Formation, Annihilation

del exciton
diffusion_average = True
periodic = True
n_proc = 2
rounds = 4
class GrainLattice:
    def assign_to_system(self, system):
        system.set_morph(np.arange(12),np.zeros(12),np.zeros(12),
                         np.repeat([0,1,2],4),cell=np.diag([12.,20.,20.]))
lattice_func = GrainLattice()
s1 = morphology.Gaussian_energy({0:(2.5,0.),1:(2.5,0.),2:(2.5,0.),"level":"s1"})
pair_i = morphology.Create_Particles("interstitial",1,morphology.randomized,mat=[0,1,2])
pair_v = morphology.Create_Particles("vacancy",1,morphology.randomized,mat=[0,1,2])
rates = {(i,j):(4. if (j-i)%3==1 else 1.) for i in range(3) for j in range(3)}
migration = Migration(k=rates)
formation = Formation(k=2.)
processes = {"interstitial":[migration,formation],"vacancy":[migration],"frenkelpair0":[]}
monomolecular = {"interstitial":[],"vacancy":[],"frenkelpair0":[Annihilation(k={0:1.,1:1.,2:1.})]}
'''
    path=tmp_path/'input.py'
    path.write_text(config)
    env=dict(os.environ)
    env['PYTHONPATH']=str(Path(__file__).resolve().parents[1])+os.pathsep+env.get('PYTHONPATH','')
    result=subprocess.run([sys.executable,'-m','kmc',str(path)],cwd=tmp_path,
                          env=env,capture_output=True,text=True,timeout=30)
    assert result.returncode==0, result.stdout+result.stderr
    with (tmp_path/'Simulation_smoke.txt').open() as handle:
        rows=[r for r in csv.DictReader(line for line in handle if not line.startswith('#'))
              if r['Type'] in ('interstitial','vacancy')]
    assert len(rows)==8
    for row in rows:
        assert abs(float(row['Drift_est_x']))>1e-3
        assert float(row['D_est_xx'])>0
        assert float(row['MSD_est_x'])==pytest.approx(
            float(row['Var_est_x'])+float(row['MeanDelta_est_x'])**2)
