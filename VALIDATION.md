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
