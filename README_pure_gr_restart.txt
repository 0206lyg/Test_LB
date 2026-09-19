pure_gr checkpoint / restart update

Apply this cumulative source update over the existing OpenLB installation.
Existing slurry/cases/pure_gr.json and Slurm resource scripts are not included
and are preserved. Local adhesion and the previously delivered convergence
fix are included in the source. Checkpoint defaults apply without JSON edits.

1. cd /home/lyjania/OpenLB
2. unzip -o OpenLB_pure_gr_restart.zip
3. sbatch build_slurry_cpu.sbatch
4. Wait for BUILD COMPLETE or BUILD REUSED.
5. New run: sbatch run_slurry_cpu.sbatch --cases pure_gr --shear-rates 100
6. Resume:  sbatch run_slurry_cpu.sbatch --restart FOLDER_NAME_UNDER_RUNS

Only runs created with this update have restart checkpoints. Old history and
particle replay files do not contain the complete fluid state.

Restart reads current pure_gr.json, retains each saved shear rate and requires
the same physical settings, timestep and MPI rank count. Total end strain,
solver limits/tolerances and output intervals can change. Results go to a new
runs folder with history copied only through the saved checkpoint boundary.

See docs/graphite-restart.md for Korean instructions and VALIDATION.md for
serial PETSc and actual two-rank OpenMPI checkpoint verification.
