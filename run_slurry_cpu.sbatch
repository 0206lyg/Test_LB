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
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
SLURRY_BASE="${SLURRY_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [[ ! -f "$SLURRY_BASE/slurry/tools/run_slurry.py" ]]; then
    printf 'Submit from your OpenLB base directory, or set SLURRY_ROOT.\n' >&2
    exit 1
fi
exec python3 "$SLURRY_BASE/slurry/tools/run_slurry.py" "$@"
