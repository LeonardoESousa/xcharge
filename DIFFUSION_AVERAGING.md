# Spatial diffusion averaging

Enable in the existing input file:

```python
diffusion_average = True
animation_mode = False
bimolec = False
```

The default is `False`. Retain the original Migration, Formation,
DissociationFP and Annihilation definitions and their time units. There is
no new fitted rate or acceleration parameter.

For each free species and occupied closed communicating class, build its
migration generator on the actual lattice and solve `pi Q = 0`. Retained
particles and material 999 sites remain obstacles. Formation channels have
weights `n_I n_V pi_I(i) pi_V(j) k_ij` on actual eligible neighbor edges.
If both species define formation, both initiating channels are retained,
exactly as in the explicit scheduler. Do not add a second channel merely to
activate averaging. Material labels and localized migration barrier profiles
are preserved; material fractions alone are never used to replace geometry.

The next slow event competes with retained monomolecular events and supported
FP generation. Its time is exponential with the total averaged propensity.
Its reactive edge is sampled proportional to that edge's weight; the original
Formation.action places the FP0 and reserves the vacancy ghost site.
A reaction, dissociation, generation or changed retained obstacle invalidates
the model. If no free defects exist, the original scheduler resumes. An
immobile vacancy is a singleton class, retaining its actual location.

Fast endpoints are sampled at slow events, with exclusion imposed when
realizing positions. Fast trajectories are not reconstructed: free-defect
Dx/Dy/Dz are NaN, and animation is unsupported. Iterations count slow events,
not skipped migration hops. Runtime is controlled by model construction and
slow-event count; rapid repeated formation/dissociation can still be costly.

## Assumptions

1. Migration reaches its stationary distribution much sooner than any slow
   event, including mixing across the localized regions. A hop-rate gap alone
   does not prove global mixing. The code does not estimate the mixing time.
2. This separation holds again after creation, dissociation and obstacle
   changes. Initial, spatially correlated transients are discarded.
3. I and V positions factorize. Birth correlations, reaction-induced
   correlations and mutual migration blocking are neglected in propensities.
4. Occupancy is dilute on each accessible class, including high-probability
   sites. Endpoint exclusion does not make the averaged propensities exact.
5. Geometry, temperature, migration profiles and process rate parameters are
   fixed between slow events. Runtime mutation of these inputs is unsupported.
6. The original local formation channels describe successful formation;
   contact alone does not force a reaction. All rates use consistent units.
7. Averaged migration classes are closed. Transient directed classes and
   unsupported free-defect spatial processes are rejected.

## Validation and its limits

`tests/test_diffusion_average.py` checks stationary balance for reversible
and nonreversible generators, finite-lattice single-pair dilution, both
initiating channels, contiguous spatial regions with migration and formation
contrast, sampled locations, executed formation and its exponential clock,
ghost restoration and model rebuilding after dissociation, fixed vacancies,
time boundaries, and rejection of unsupported processes.

An independent exact continuous-time two-particle chain test includes every
migration hop and site exclusion. On a 30-site uniform open chain, the
independent averaged mean waiting time differs from the exclusion model by
approximately 3.45% in the fast-migration limit. The discrepancy is the known
`N/(N-1)` dilution correction, not sampling noise. These tests validate the
implementation and a controlled small model; they do not establish accuracy
for a particular production morphology or density.

Methodological basis: Cao, Gillespie and Petzold, *The slow-scale stochastic
simulation algorithm*, J. Chem. Phys. 122, 014116 (2005),
https://doi.org/10.1063/1.1824902. The independent spatial approximation here
is an additional assumption, not a claim of exact lattice averaging.

## Effective diffusion and output columns

For periodic boundaries in all three directions use the existing setting:

```python
periodic = True
diffusion_average = True
```

No boundary option or additional diffusion coefficient is required. For each
species and closed migration class, the code computes the long-time unwrapped
diffusion tensor from the original rates and edge vectors. It solves the
centered lattice corrector equation

`Q chi = v - b`, with `b_i = sum_j w_ij d_ij` and `v = sum_i pi_i b_i`,

then evaluates

`D = (1/2) sum_ij pi_i w_ij (d_ij + chi_j - chi_i)(d_ij + chi_j - chi_i)^T`.

Thus localized barrier regions, their connectivity and return correlations
enter the result. Nonzero macroscopic drift is supported and reported separately from the
centered diffusion tensor. Finite confined classes have zero
long-time diffusion in bounded directions; this does not mean zero local
mobility. Small graphs use sparse LU; large reversible graphs use a checked,
preconditioned iterative solve. A tensor is reused until the averaged model
changes.

For each individual free defect, accumulate `J += D * dt`, `M += v * dt`, and
`AveragedDiffusionTime += dt` before executing the next slow event or reaching
a deterministic time boundary. New defects begin with zero accumulation,
including defects born through FP dissociation. The original CSV log gains:

| Columns | Meaning |
| --- | --- |
| `AveragedDiffusionTime` | Time spent in averaged free-defect evolution |
| `D_est_xx`, `D_est_yy`, `D_est_zz`, `D_est_xy`, `D_est_xz`, `D_est_yz` | Components of `J / AveragedDiffusionTime` |
| `MSD_est_x`, `MSD_est_y`, `MSD_est_z` | Diagonal components of `2 J + M M^T` |
| `Cov_est_xy`, `Cov_est_xz`, `Cov_est_yz` | Off-diagonal components of centered covariance `2 J` |
| `Drift_est_x`, `Drift_est_y`, `Drift_est_z` | Lifetime-averaged velocity `M / AveragedDiffusionTime` |
| `MeanDelta_est_x`, `MeanDelta_est_y`, `MeanDelta_est_z` | Integrated mean displacement `M` |
| `Var_est_x`, `Var_est_y`, `Var_est_z` | Centered directional variance, diagonal components of `2 J` |

D units are lattice-coordinate units squared per simulation time unit; the
MSD and covariance units are lattice-coordinate units squared. For coordinates
in Angstrom and times in ps these are Angstrom squared/ps and Angstrom squared.
Particles never averaged, such as retained FP products, report NaN in these
columns. The explicit displacement columns remain NaN for averaged defects.
The added columns are emitted only with diffusion_average enabled. Appending
an incompatible older CSV schema is rejected instead of mixing row formats.

These are migration-model estimates, not reconstructed trajectories or an
independent measurement of D. Their conversion to an expected lifetime MSD
assumes the long-time diffusive regime and neglects reaction conditioning and
transients after environment changes. They do not specify the winding numbers
of the sampled reaction endpoint. For finite nonperiodic geometries they are
asymptotic transport covariances, not finite-time confined MSDs.

`tests/test_diffusion_tensor.py` checks homogeneous and heterogeneous periodic
rings against arithmetic/harmonic analytic results, a destination-dependent
localized barrier profile, off-diagonal components, a fully periodic cubic
lattice, interval accumulation with changed obstacles, reset for new defects,
zero diffusion for a confined chain, biased walks, and the large-graph
iterative solver. The CLI test reads the new CSV columns for a fully periodic
3D lattice, checks their numerical values and handles incompatible append
schemas. The existing direct-KMC CLI test checks that its original log schema
continues to work.

Transport reference: Chapman and Shabala, *Effective Transport Properties of
Lattices*, SIAM J. Appl. Math. 77 (2017),
https://doi.org/10.1137/16M1092039.


## Nonzero-drift failure in a grain-boundary run

The earlier evaluator rejected a stationary drift above a relative tolerance.
A traceback at that check identifies the failed assumption, not its physical
cause. Localized grain boundaries and periodic boundaries alone do not prove
zero drift. Drift can reflect nonreversible rate rules, a directed neighbor
graph, inconsistent periodic edge vectors, or numerical error in the stationary
solve. Symmetric rates on a reciprocal graph with antisymmetric edge vectors
should give zero drift. The actual input and morphology are needed to diagnose
a particular run.

The corrected evaluator centers the Poisson equation on the computed drift
and keeps the existing stationary and corrector residual checks. It does not
change the microscopic migration rates, force reversible edges, or alter
formation propensities. It reports the drift rather than folding ballistic
transport into D. To recover D from accumulated second moments use Var_est,
not MSD_est when MeanDelta_est is nonzero. Existing zero-drift results have
unchanged interpretation. The new log schema cannot be appended to an old
header; use a new output or append_output=False.

Additional checks cover the known biased-ring result
`v = (k_plus - k_minus) * step` and
`D = (k_plus + k_minus) * step^2 / 2`, a heterogeneous biased cycle checked
against independent tilted-generator eigenvalue derivatives, drift accumulation
over changing intervals, and a two-worker CLI run with contiguous material
regions and nonzero drift. These verify the extension; they do not reproduce
the user's unprovided grain-boundary morphology.

## Migration graph, angular balance, and periodic-image fixes

Migration now has a separate `migration` neighbor mode. Its geometric edge set
is the union of the material-resolved first-shell selections and their exact
reverse edges. Existing accepted migration edges are retained; missing reverse
geometric edges are added. Original base-rate dictionaries still determine
which channels have nonzero rates. Reaction and generation partner lists retain
their original per-origin first-shell rule and count each wrapped partner once.

Each physical periodic-image migration is a separate channel, even when two
images share a wrapped destination. Hops to a translated image of the current
site are also retained. Their propensities contribute to escape/event clocks;
the wrapped-site generator sums these rates, while transport uses the unmerged
physical jump vectors. The explicit Migration action samples the physical image
conditional on its selected destination using the same rates. Endpoint-only
minimum-image ties have a deterministic reciprocal convention. This avoids
both false drift and the loss of hop multiplicity in small periodic cells.

Scalar site barrier profiles now default to `directional=False`. Their existing
nondirectional rate formulas are preserved. In explicit directional destination
mode, let U_i be the profile, c_ij = abs(cos(theta))**angular_exponent, and define

`delta_ij = (U_j - U_i)/2 + (U_i + U_j)*c_ij/2`.

Only the symmetric bond part is angularly weighted. For reciprocal physical
edges and symmetric base rates, the reverse ratio is exp(-(U_j-U_i)/kBT),
so stationary occupancy is proportional to exp(-U_i/kBT) and equilibrium drift
is zero. Along the selected axis c=1 recovers the previous destination formula.
In `coupling="symmetric"` mode, the profile supplies only the unweighted half
site-energy difference; the directional flag does not multiply that difference.
This changes the previously nonconservative directional convention deliberately.

Image-specific barrier cache entries retain their source arrays and are
recovered from topology after worker pickling. Cached or fallback calculations
therefore use each physical image vector, not one representative vector for all
copies of a destination.

The regression suite reproduces the former missing-reverse-edge examples,
tests passive grain profiles on irregular periodic graphs, verifies detailed
balance on an oblique winding cycle, checks 1-site/2-site/4-site periodic cube
transport, exercises explicit image selection and reciprocal endpoint ties,
and checks the barrier cache after worker serialization. The original reaction
and generation tests continue to pass. Genuine bias explicitly present in the
base rates remains supported through the separate drift columns.
