# Adhesive-network convergence regression

`adhesive_network_108.dat` is the unmodified rank-0 replay from MIT job
`slurry_23075858`, also byte-identical to the supplied rank-10 replay. It retains
108 particles, alpha=1, the original forces, contact history, timestep, physical
acceptance tolerances and iteration limits. It records the failed substep at
LB step 3264, subdivision count 256, index 112, shear rate 10/s.

From the repository root, with the PETSc environment loaded:

```bash
python3 slurry/tools/compile_petsc_test.py tests/petsc_contact/replay_particle_step.cpp --run -- tests/petsc_contact/fixtures/adhesive_network_108.dat
```

Add `--compiler g++ --serial` before `--run` for a serial MPIUNI PETSc install.
The command must exit successfully and report `REPLAY ACCEPTED`. It replays one
recorded particle substep, not the preceding LB trajectory or the remaining
part of that outer step.

With the corrected default linear relative tolerance 0.001, PETSc 3.25.5
converges in 1 Newton step / 6 Krylov iterations. Force and torque residual
ratios are approximately 0.2413 and 0.2025, both below the unchanged limit 1.
The prior 0.05 default stalls at ratios 3.1189 and 4.9400 after 60 Newton steps.
Iteration counts may vary with compiler and PETSc version; physical acceptance
is the regression criterion.
