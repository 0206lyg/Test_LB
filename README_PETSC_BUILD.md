# PETSc particle-solver build

The only build entry point is the existing `build_slurry_cpu.sbatch`. On MIT,
run it from `/home/lyjania/OpenLB` after applying the update:

```bash
sbatch build_slurry_cpu.sbatch
```

Wait for `BUILD COMPLETE` in `build_slurry_<jobid>.out`, then use the existing run command:

```bash
sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
```

The first job reuses a compatible PETSc installation if available. Otherwise it
downloads official PETSc 3.25.5 and builds it under
`$HOME/.local/slurry-petsc/3.25.5/<toolchain-id>/install` using the same MPI
compiler wrappers as OpenLB. It needs internet access for PETSc and the
`f2cblaslapack` dependency on this first installation. It does not need sudo and
does not install another MPI implementation. Later builds reuse the installation.

The existing GCC 12.2.0/OpenMPI 4.1.4 module defaults and `mit_normal` partition
are retained. On a different cluster, provide its actual module names through
`SLURRY_GCC_MODULE` / `SLURRY_MPI_MODULE` and the appropriate `sbatch --partition`
argument. Use the same modules for dependency installation, application build,
and simulation. The scripts do not guess a PETSc module name.

PETSc preparation is an internal stage of `slurry/tools/build_slurry.py`;
it is not a second build job. The one job prepares dependencies, compiles all
registered cases, and publishes the existing `build/slurry/current/slurry`.
If downloads are blocked, a compatible PETSc installation must be made available
through `PETSC_DIR` before submitting this same build job.

An existing PETSc prefix can be selected with `PETSC_DIR`; for an in-place PETSc
build also set `PETSC_ARCH`. Otherwise `pkg-config` discovers it. An incompatible
visible installation causes an explicit compiler/MPI preflight failure. The default official HTTPS archive's computed SHA256 is recorded in the setup
manifest; it is not represented as an independently authenticated expected hash.

Builds default to the PETSc backend and define `SLURRY_USE_PETSC`. A controlled
comparison build can use `sbatch build_slurry_cpu.sbatch --particle-solver legacy`.
The executable `--build-info` and `build_manifest.json` identify the backend.
The PETSc version, exact flags, compiler/wrapper identity, preflight output, and
configuration/library hashes enter the build fingerprint. Link commands include
runtime paths to the discovered PETSc libraries; MPI modules must still match.
Changed or missing PETSc dependencies are rejected by the run preflight.

Real double-precision PETSc 3.18 or later is required. The MPI build rejects
MPIUNI. Before compiling OpenLB, a small program checks PETSc initialization,
SNES creation, linking, and the MPI runtime with the selected C++ compiler.
The MPI runtime check runs through that MPI installation's `mpirun -np 1`
(or `mpiexec`), including inside Slurm. A completed PETSc installation is reused
if an earlier attempt failed only at this runtime check.
See the [Open MPI 4.x Slurm FAQ](https://www.open-mpi.org/faq/?category=slurm).

To compile an independent C++ solver check using the same dependency discovery:

```bash
source build/petsc/env.sh
python3 slurry/tools/compile_petsc_test.py path/to/test.cpp --run
```

This update requires a C++ rebuild. It does not modify the selected physical
force parameters or introduce checkpoint/restart behavior.

Official installation references:
[PETSc installation](https://petsc.org/release/install/install/) and
[PETSC_DIR/PETSC_ARCH](https://petsc.org/release/install/multibuild/).
