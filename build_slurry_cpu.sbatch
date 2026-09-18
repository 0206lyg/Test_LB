#!/bin/bash
#SBATCH --job-name=build_slurry
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=build_slurry_%j.out

set -Eeuo pipefail
module purge
module load "${SLURRY_GCC_MODULE:-gcc/12.2.0}"
module load "${SLURRY_MPI_MODULE:-openmpi/4.1.4}"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
SLURRY_BASE="${SLURRY_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "$SLURRY_BASE/slurry/tools/build_slurry.py" ]]; then
    printf 'Submit from your OpenLB base directory, or set SLURRY_ROOT.\n' >&2
    exit 1
fi
if [[ -z "${PETSC_DIR:-}" && -z "${PETSC_ARCH:-}" && -f "$SLURRY_BASE/build/petsc/env.sh" ]]; then
    source "$SLURRY_BASE/build/petsc/env.sh"
fi
# The common build prepares PETSc and compiles the single all-case executable.
exec python3 "$SLURRY_BASE/slurry/tools/build_slurry.py" --mode mpi "$@"
