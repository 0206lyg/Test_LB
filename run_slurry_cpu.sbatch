#!/bin/bash
#SBATCH --job-name=slurry
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@180
#SBATCH --output=slurry_%j.out

set -Eeuo pipefail
module purge
module load "${SLURRY_GCC_MODULE:-gcc/12.2.0}"
module load "${SLURRY_MPI_MODULE:-openmpi/4.1.4}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" PYTHONUNBUFFERED=1
if [[ ! "$OMP_NUM_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'SLURM_CPUS_PER_TASK must be a positive integer.\n' >&2
    exit 1
fi
if [[ "$OMP_NUM_THREADS" != 1 ]]; then
    # Allocate every rank its requested cores, rather than binding all of its
    # OpenMP threads to one core. The default 32 ranks x 1 thread is unchanged.
    export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
    # Open MPI 4 uses ORTE; Open MPI 5 uses the renamed PRRTE parameters.
    export OMPI_MCA_rmaps_base_mapping_policy="slot:PE=$OMP_NUM_THREADS"
    export OMPI_MCA_hwloc_base_binding_policy=core
    export PRTE_MCA_rmaps_default_mapping_policy="slot:PE=$OMP_NUM_THREADS"
    export PRTE_MCA_hwloc_default_binding_policy=core
fi
SLURRY_BASE="${SLURRY_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "$SLURRY_BASE/slurry/tools/run_slurry.py" ]]; then
    printf 'Submit from your OpenLB base directory, or set SLURRY_ROOT.\n' >&2
    exit 1
fi
if [[ -z "${PETSC_DIR:-}" && -z "${PETSC_ARCH:-}" && -f "$SLURRY_BASE/build/petsc/env.sh" ]]; then
    source "$SLURRY_BASE/build/petsc/env.sh"
fi
exec python3 "$SLURRY_BASE/slurry/tools/run_slurry.py" "$@"
