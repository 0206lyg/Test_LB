# Validation of the PETSc particle solver update

Base: `0206lyg/Test_LB` commit `ad0befc0f75b40e16bbd81fc3f29c53e76a8b551`.

## Environment

- GCC/G++ 13.3, Linux, C++20.
- Actual PETSc 3.25.5, real double precision, MPIUNI/serial build.
- Full OpenLB source from the repository above.
- Slurm and the user's GCC 12.2/OpenMPI 4.1.4 cluster were not available here.

## Completed checks

| Check | Result |
|---|---|
| Full OpenLB application compilation and PETSc linking | Passed |
| Build source guard, release manifest and `build/slurry/current` publication | Passed (`BUILD COMPLETE`) |
| Executable `--build-info` | `particle_solver=petsc`, `petsc_version=3.25.5`, `mpi_enabled=false` |
| Reduced integration run | 4 particles, 80³ grid, 100/s, 7 LB steps, strain 0.00202072594, status `COMPLETED` |
| Common batch controller with the built executable | Passed (`BATCH COMPLETED`); diagnostic summary script copied into the result folder |
| Standalone physical solver regressions | 14/14 passed |
| Diagnostic summarizer regressions | 10/10 passed |
| Installer preservation/backup/input validation fixtures | 10/10 passed |
| Saved failure replay | Same failure and physical residuals reproduced |
| Failed fixture with Newton cap changed from 1 to 80 | Converged in 37 iterations; comparison fixture only |
| MPI compatibility preflight | Correctly rejected serial MPIUNI PETSc for an MPI application |

The 14 physical fixtures cover free flight, forced motion, normal equilibrium,
approaching contact, separation/history release, sliding and rolling history,
elastic/yielded rolling dynamics, rigid rotation of history, friction/rolling
yield, a three-body contact chain, failed-step rollback, recovered-retry output
suppression, and production oblate contact. The failure fixture also checks the
exact replay roundtrip for bodies, contact history, cache, loads, and settings.

The production contact fixture uses graphite oblate axes (1.65, 1.65, 0.20) µm,
the production inertia scaling and timestep, Hamaker 9.9e-20 J, water lubrication,
relative tolerance 0.01, absolute force tolerance 1e-13 N and torque tolerance
1.65e-19 N m. Three successive steps pass with adhesive rolling history.

The reduced integration run checks application wiring; its small particle set is
not a rheology benchmark. It does not reproduce the dense 108-particle failure.

The files defining RE² forces, lubrication, rough friction/adhesive rolling,
ellipsoid gap geometry and particle math were compared byte-for-byte with the
base source and remain unchanged. The installer preserves existing physical and
numerical settings in `pure_gr.json` apart from the explicit backend/diagnostic
fields and the new Krylov default if absent.

## Reproduce the standalone tests

On the cluster, use the same modules and PETSc installation as the build job:

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/contact_solver_tests.cpp --run
python3 -m unittest discover -s slurry/tests -p 'test_particle_diagnostics.py'
```

For an intentionally serial PETSc installation, the C++ test command needs
`--compiler g++ --serial`.

## Remaining verification

The user's multi-rank MPI application and full 108-particle production trajectory
have not been run here. Their convergence, long-run stress statistics and
performance need the rebuilt executable on the user's cluster. PETSc supplies
the nonlinear/linear solution framework; the contact residual and preconditioner
remain application code. This update does not establish convergence for every
contact-network transition or validate the physical model against experiment.

## Single build entry correction

The existing `build_slurry_cpu.sbatch` now invokes dependency preparation from
inside the common Python build. Only this one build job is submitted. The extra
root-level build and setup wrappers have been removed.

Ten focused checks with mocked compiler/build processes passed: shell dispatch,
root detection, saved environment reuse, explicit PETSc environment preservation,
MPI dependency setup ordering, job-count forwarding, cache reuse, and legacy/serial
mode routing. Fresh installation and replacement of the previous package also
passed, including obsolete-script backup/removal and physical-setting preservation.
No numerical solver code changed in this correction; C++ edits only correct an
error-message build command. The earlier real PETSc/full serial build checks above
are from the solver update; MPI compilation was not rerun for this script correction.

## MIT MPI preflight correction

User log `build_slurry_22944246.out` shows successful PETSc installation followed
by MPI_Init_thread failure when the probe was executed directly under Slurm.
The probe now uses `mpirun -np 1` (or `mpiexec`) from the selected compiler's
MPI installation. Compiler/link and runtime failures are reported separately.
The optional standalone test runner uses the same launch command.

Eight targeted regression checks passed, including launcher selection, Slurm
one-rank launch routing, stable build fingerprints, error classification, serial
mode, and reuse of a completed install after its previous check failed. MPI calls
in these regressions were mocked; actual MIT Slurm execution remains to be run.
The real serial PETSc 3.25.5 compile/link/runtime probe passed again.

The installed prefix is reused with the same compiler/MPI modules, even when the
previous attempt failed before writing `build/petsc/env.sh`. No force-law or solver
algorithm changes were made for this correction.

## Graphite local-gap adhesion

This change starts from `642ffc182930c9d9fe5114bf77ab99951dd30009` and explicitly
changes the graphite pair potential and the `pure_gr` physical settings. The
unchanged-force statements in the historical sections above describe those
earlier solver/build updates, not this adhesion update.

The production case now uses sigma = 0.4197 nm, D0 = 0.30 nm, local fraction = 1,
and tangential stiffness = 80 N/m. Hamaker, geometric roughness gap, sliding
friction, rolling length/angle, particle dimensions, timestep and solver
tolerances retain their previous values. The local potential smoothly replaces
the near-contact contribution, including both attraction and repulsion; forces
and torques differentiate the complete energy. See
[the model and run instructions](docs/graphite-local-adhesion.md).

Checks used GCC/G++ 13.3.0 and real, double-precision PETSc 3.25.5 with serial
MPIUNI. No MPI compiler was available in this environment.

| Check | Result |
| --- | --- |
| Standalone adhesion physics groups | 6/6 passed |
| PETSc contact regression fixtures | 25/25 passed |
| Local-gap solver domain and replay checks, linked with PETSc | Passed |
| Python configuration/build/diagnostic tests | 26/26 passed |
| Complete shared OpenLB executable, serial PETSc build | Passed |
| Executable capability metadata | `local_gap_adhesion: true`, PETSc 3.25.5 |
| Actual OpenLB integration: two FF particles, seven LB steps | Passed |

The physics groups check legacy alpha=0 results, contact-force magnitudes,
translation/rotation finite differences of energy, total angular momentum,
smooth local-switch release, parameter/domain rejection, and the independent
effects of tangential stiffness and sliding/rolling caps. The FF result is
936.214833 nN with the rounded input sigma, consistent with the small-gap
Derjaguin limit. EF and EE forces are 25.950663 nN and 13.748022 nN.

The new PETSc fixture advances the production oblate pair with water lubrication,
the production timestep/tolerances, and corrected adhesive rolling history.
The domain/replay checks exercise residual-domain rejection, unchanged committed
history after failed trials, corrected contact-birth strength, exact new-format
replay, and old-format replay with local adhesion disabled.

The actual OpenLB integration used the new `pure_gr` physical settings with a
temporary two-particle, 80-cubed-grid fixture (125 nm spacing), shear rate 100/s,
and end strain 0.002. The final strain was 0.00202073 after seven LB steps, all
accepted with one particle substep. All eight sampled history rows were finite.
The final pair force was 936.2185 nN and gap 1.99999999971 nm; the maximum gap
violation was 9.94e-14 m against the unchanged 1e-12 m tolerance. The maximum
accepted normalized force residual was 0.9784, below its acceptance limit of 1.
Output metadata confirmed all new interaction settings and retained contact
parameters. This fixture checks the actual application path; it is not a bulk
rheology result and does not change the committed production case geometry.

The full 108-particle trajectory and multi-rank MPI execution have not been run
for this change. These checks establish implementation and integration behavior;
they do not establish a 500 Pa bulk yield stress.

## MIT job 23075858: adhesive-network Newton stagnation

The supplied rank-0 and rank-10 replay snapshots are byte-identical. They record
the 10/s run at LB step 3264, with the terminal failure at subdivision count 256,
index 112 and substep duration 4.5105489780439514e-8 s. The previously accepted
LB state is at time 0.037689425572698769 s. No MPI rank-state disagreement is
present in those two snapshots.

The failure was reproduced with real, double-precision PETSc 3.25.5, serial
MPIUNI. The default linear GMRES relative tolerance of 0.05 accepted the first
Newton direction after two Krylov iterations, leaving scaled linear residual
16.9317. The resulting nonlinear path stalled at force/torque ratios
3.11888/4.93999, despite a passing complementarity ratio 0.812118 and no domain
rejections. Increasing the number of identical stalled nonlinear iterations
does not address that inaccurate direction.

Only the default linear GMRES relative tolerance was changed to 0.001. The first
direction then used six Krylov iterations, leaving scaled linear residual
0.836586, and the recorded substep converged after one Newton step. Force/torque
ratios were 0.2413223374/0.2025331374 and gap violation 8.4035373e-13 m. Both
uploaded snapshots passed. Alpha remained 1; pair/contact physics, force/torque/
gap acceptance tolerances, substep duration and all iteration limits remained
unchanged. The NGMRES and Newton line-search architecture was retained.

The exact state is retained as
`tests/petsc_contact/fixtures/adhesive_network_108.dat`, with reproduction
instructions alongside it. This verifies the formerly failing particle
substep; it does not rerun the entire MPI trajectory.

An additional continuation carried the accepted bodies, contact history and gap
cache through all 144 remaining particle substeps, indices 112 through 255,
using the original frozen hydrodynamic force/torque and substep duration. All
144 passed, reaching the original LB-step endpoint. The worst accepted force
ratio was 0.634466, torque ratio 0.759672 and gap violation 8.4035373e-13 m.
No substep used more than 19 Newton iterations against its unchanged limit 60.

After the fix, all 25 PETSc contact regression fixtures and all 26 Python tests
passed, and the full shared OpenLB executable rebuilt successfully with serial
PETSc. The diagnostics retry test's deliberately injected Newton budget was
changed from 3 to 2 because the more accurate solver now converges without a
retry at 3; this test-only setting continues to check actual retry/no-output
behavior and does not change the production configuration.

The attached history belongs to the 10/s run and ends at sampled stress
243.99 Pa; it does not contain the separately reported 100/s stress result.
It also records transient fluid Mach numbers up to 0.622 and fluid density
deviation up to 39.0%. Those are a separate limitation of the fluid time/velocity
scaling; the saved-step replay isolates the particle convergence fix and does
not demonstrate that these fluid transients have been resolved.

## pure_gr coupled checkpoint / restart

The restart update saves the complete serializable OpenLB lattice on every
rank, particle bodies, pair cache, persistent sliding/rolling contact state,
last angular acceleration, diagnostic state, accumulated timers and LB step.
Lees–Edwards time is recovered from the saved step and unchanged timestep.
Loading clears the restored particle auxiliary fields before rebuilding the
process-local mapping cache; omitting this would double the particle mask.
Completed checkpoints are published by directory rename after every rank has
closed its file. The latest two completed checkpoints are retained by default.

Checks used GCC/G++ 13.3 and the existing PETSc 3.25.5 MPIUNI installation.
A separately installed OpenMPI 4.1.4 toolchain was also located and used for
real two-rank compilation and execution of the checkpoint path with the
legacy particle backend. The MPI tests therefore exercise actual rank-local
lattice files and collective publication/restoration; they do not constitute
a PETSc-with-MPI or MIT Slurm production test.

| Check | Result |
| --- | --- |
| Final complete OpenLB/PETSc serial build | Passed |
| Complete OpenLB/OpenMPI build, legacy particle backend | Passed |
| New restart Python tests | 10/10 passed |
| Full Python configuration/build/diagnostic suite | 36/36 passed |
| Serial PETSc: continuous 16 steps vs. 8 steps + restart to 16 | All physical CSV fields exactly equal; final lattice bytes identical |
| Two-rank OpenMPI: continuous 16 steps vs. 8 steps + restart to 16 | All physical CSV fields exactly equal; both rank lattice files identical |
| SIGUSR1 through common controller and Python driver, serial and two ranks | Clean CHECKPOINTED exit after a complete LB step; subsequent restart passed |
| Deliberately corrupted state byte with unchanged file length, serial and two ranks | C++ checksum rejected restoration |

The equivalence fixture contains two touching production-size oblate particles
on an 80³ grid at 100/s, with the configured local adhesion and persistent
contact elasticity. Contact count and nonzero stored contact elastic energy
are asserted, so the test is not merely a free-particle restart. All physical
history columns and particle CSV rows are compared, including the copied
prefix; only elapsed/performance timing columns are excluded. The same common
controller invoked by `run_slurry_cpu.sbatch` performs the restart.

Unit checks cover incomplete staging directories, missing rank files, old
history-only runs, physical/timestep/rank compatibility, cumulative endpoints,
truncation to the saved CSV boundary, and mixed-rate batches that skip already
completed runs while preserving the original rate of each resumed run.

Reproduction and user commands are in [graphite-restart.md](docs/graphite-restart.md).
The old production runs have no coupled checkpoint and cannot be recovered
exactly from their diagnostic/history files. Native checkpoint version 1
requires the same OpenLB data layout, physical settings, timestep and MPI rank
count. Full 32-rank production execution and Slurm's actual timeout delivery
have not been tested here.
