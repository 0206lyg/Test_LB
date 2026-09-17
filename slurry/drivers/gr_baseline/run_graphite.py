#!/usr/bin/env python3
"""Prepare and run one graphite Couette case; Python 3 standard library only.

Use --dry-run to inspect units and estimated steps without launching OpenLB.
Normal cluster runs belong inside a Slurm allocation. All inputs are recorded;
restarts continue into a new output directory, leaving the previous run intact.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1048576), b""):
            result.update(chunk)
    return result.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def finite_positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")
    return value


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))
    return value


def validate_restart_manifest(checkpoint, executable, ranks):
    """Require the original run manifest to guard binary/ABI and MPI continuity."""
    source = checkpoint.parent / "manifest.json"
    if not source.is_file():
        raise ValueError("Restart requires manifest.json in the checkpoint's parent run directory; retain the complete run provenance")
    previous = json.loads(source.read_text(encoding="utf-8"))
    if previous.get("sha256", {}).get("executable") != digest(executable):
        raise ValueError("Restart executable SHA256 differs from the original run; use the exact original binary")
    if previous.get("ranks") != ranks:
        raise ValueError("Restart requires the original MPI rank count")
    if not previous.get("sha256", {}).get("particles"):
        raise ValueError("Original run manifest has no initial-particle checksum")
    return previous


def resolve(config, shear_rate=None, max_steps=0):
    """Validate inputs and derive dimensional metadata without starting a run."""
    cfg = copy.deepcopy(config)
    if cfg.get("schema_version") != 1:
        raise ValueError("Only config schema_version=1 is supported")
    if shear_rate is not None:
        cfg["flow"]["shear_rate_s_inv"] = shear_rate
    f, g, p, flow, n, out = (cfg[key] for key in ("fluid", "geometry", "particles", "flow", "numerics", "output"))
    for group_name, group, names in (
        ("fluid", f, ("density_kg_m3", "dynamic_viscosity_Pa_s", "temperature_K")),
        ("geometry", g, ("box_length_m", "dx_m")),
        ("particles", p, ("diameter_m", "thickness_m", "density_kg_m3")),
        ("flow", flow, ("shear_rate_s_inv", "end_strain")),
        ("numerics", n, ("nu_lattice", "epsilon_cells", "contact_young_Pa")),
    ):
        for name in names:
            finite_positive(group[name], group_name + "." + name)
    integer(p["count"], "particles.count")
    integer(p["seed"], "particles.seed")
    integer(n["contact_resolution"], "numerics.contact_resolution", 4)
    integer(max_steps, "max_steps")
    for name in ("vtk_every_steps", "checkpoint_every_steps"):
        integer(out[name], "output." + name)
    integer(out["sample_every_steps"], "output.sample_every_steps", 1)
    for name, value in (("minimum_gap_m", p.get("minimum_gap_m", 0)),
                        ("target_re", n["target_re"]),
                        ("contact_enlargement_cells", n["contact_enlargement_cells"])):
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError(name + " must be finite and nonnegative")
    if not 0 <= n["contact_poisson"] < 0.5:
        raise ValueError("contact_poisson must lie in [0, 0.5)")
    if n["contact_restitution"] != 1:
        raise ValueError("This frictionless elastic baseline requires contact_restitution=1")
    dx, length, diameter, thickness, rate = g["dx_m"], g["box_length_m"], p["diameter_m"], p["thickness_m"], flow["shear_rate_s_inv"]
    intervals = round(length / dx)
    if intervals < 8 or not math.isclose(length / dx, intervals, rel_tol=0, abs_tol=1e-7):
        raise ValueError("box_length_m / dx_m must be an integer >= 8")
    if not 0 < thickness <= diameter < length:
        raise ValueError("Need 0 < thickness <= diameter < box length")
    if p["count"] > 1 and 2 * diameter >= length:
        raise ValueError("Multiple particles require 2*diameter < box length: native contacts retain one periodic image per particle pair")
    volume = math.pi * diameter * diameter * thickness / 6
    phi = p["count"] * volume / length ** 3
    if not 0 <= phi < 1:
        raise ValueError("Particle volume fraction must be < 1")
    physical_nu = f["dynamic_viscosity_Pa_s"] / f["density_kg_m3"]
    physical_re = rate * diameter ** 2 / physical_nu
    if n["target_re"] > 0:
        dt = n["target_re"] * n["nu_lattice"] * dx ** 2 / (rate * diameter ** 2)
    else:
        dt = n["nu_lattice"] * dx ** 2 / physical_nu
    rho_inertial = f["dynamic_viscosity_Pa_s"] * dt / (n["nu_lattice"] * dx ** 2)
    requested_steps = math.ceil(flow["end_strain"] / (rate * dt))
    stop_steps = min(requested_steps, max_steps) if max_steps else requested_steps
    mach = 0.5 * rate * length * dt / dx * math.sqrt(3)
    if mach >= 0.1:
        raise ValueError("Wall Mach >= 0.1; reduce target_re or time step")
    meta = {
        "shear_rate_s_inv": rate, "dt_s": dt, "tau_lattice": 0.5 + 3 * n["nu_lattice"],
        "physical_fluid_re_D": physical_re,
        "numerical_fluid_re_D": rate * diameter ** 2 * rho_inertial / f["dynamic_viscosity_Pa_s"],
        "physical_particle_St_D": rate * diameter ** 2 * p["density_kg_m3"] / f["dynamic_viscosity_Pa_s"],
        "numerical_particle_St_D": rate * diameter ** 2 * p["density_kg_m3"] * rho_inertial / f["density_kg_m3"] / f["dynamic_viscosity_Pa_s"],
        "fluid_inertial_density_kg_m3": rho_inertial,
        "particle_inertial_density_kg_m3": p["density_kg_m3"] * rho_inertial / f["density_kg_m3"],
        "inertial_density_scale": rho_inertial / f["density_kg_m3"],
        "particle_volume_m3": volume, "actual_volume_fraction": phi,
        "actual_mass_fraction": phi * p["density_kg_m3"] / (phi * p["density_kg_m3"] + (1 - phi) * f["density_kg_m3"]),
        "grid_intervals_per_direction": intervals, "nominal_bulk_cells": intervals ** 3,
        "d3q19_single_population_bytes": intervals ** 3 * 19 * 8,
        "thickness_cells": thickness / dx, "diameter_cells": diameter / dx,
        "wall_mach": mach, "steps_to_requested_end_strain": requested_steps,
        "this_run_max_step": stop_steps, "this_run_target_strain": stop_steps * rate * dt,
        "notes": [
            "Bulk cell count excludes boundary and MPI overlap cells; memory is population storage only.",
            "End strain is a run budget, not a claim of rheological steady state.",
            "target_re > 0 changes numerical inertia while preserving physical dynamic viscosity; validate inertia convergence.",
            "No attraction; frictionless elastic normal contact limits overlap and requires stiffness/resolution checks.",
        ],
    }
    return cfg, meta


def solver_values(cfg, output, particles, restart, max_steps):
    f, g, p, flow, n, out = (cfg[key] for key in ("fluid", "geometry", "particles", "flow", "numerics", "output"))
    values = {
        "shear_rate": flow["shear_rate_s_inv"], "box_x": g["box_length_m"],
        "box_y": g["box_length_m"], "box_z": g["box_length_m"], "dx": g["dx_m"],
        "diameter": p["diameter_m"], "thickness": p["thickness_m"],
        "rho_particle": p["density_kg_m3"], "rho_fluid": f["density_kg_m3"],
        "dynamic_viscosity": f["dynamic_viscosity_Pa_s"],
        "nu_lattice": n["nu_lattice"], "target_re": n["target_re"],
        "epsilon_cells": n["epsilon_cells"], "contact_enlargement_cells": n["contact_enlargement_cells"],
        "contact_young": n["contact_young_Pa"], "contact_poisson": n["contact_poisson"],
        "contact_restitution": n["contact_restitution"], "contact_resolution": n["contact_resolution"],
        "end_strain": flow["end_strain"], "max_steps": max_steps,
        "sample_every": out["sample_every_steps"], "vtk_every": out["vtk_every_steps"],
        "checkpoint_every": out["checkpoint_every_steps"], "output_dir": str(output),
        "particles_csv": str(particles), "restart_dir": str(restart) if restart else "",
    }
    return values


def main():
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--config", type=Path, default=Path(__file__).with_name("config.json"))
    command.add_argument("--shear-rate", type=float)
    command.add_argument("--dry-run", action="store_true", help="Print resolved parameters without creating files or launching a simulation")
    command.add_argument("--output", type=Path)
    command.add_argument("--executable", type=Path, default=Path(__file__).with_name("build") / "current" / "graphiteCouette3d")
    command.add_argument("--generator", type=Path, default=Path(__file__).with_name("generate_particles.py"))
    command.add_argument("--ranks", type=int, default=1)
    command.add_argument("--benchmark-steps", type=int, default=0, help="Absolute maximum step, including steps already saved in a restart")
    command.add_argument("--restart", type=Path, help="Completed checkpoint directory from an earlier run; same geometry, material, MPI layout")
    args = command.parse_args()
    try:
        cfg, meta = resolve(json.loads(args.config.read_text(encoding="utf-8")), args.shear_rate, args.benchmark_steps)
        integer(args.ranks, "ranks", 1)
        if args.dry_run:
            print(json.dumps({"config": cfg, "derived": meta}, indent=2, allow_nan=False))
            return 0
        if not os.environ.get("SLURM_JOB_ID") and not os.environ.get("GRAPHITE_ALLOW_LOCAL"):
            raise ValueError("Submit run_graphite_cpu.sbatch; local validation may explicitly set GRAPHITE_ALLOW_LOCAL=1")
        if args.output is None:
            raise ValueError("--output is required for a simulation")
        output, executable, generator = args.output.expanduser().resolve(), args.executable.expanduser().resolve(), args.generator.expanduser().resolve()
        restart = args.restart.expanduser().resolve() if args.restart else None
        if not executable.is_file() or not os.access(str(executable), os.X_OK):
            raise ValueError("Missing executable: " + str(executable))
        if restart and not restart.is_dir():
            raise ValueError("Restart directory not found: " + str(restart))
        previous = validate_restart_manifest(restart, executable, args.ranks) if restart else None
        if any((output / name).exists() for name in ("manifest.json", "resolved_run.cfg", "effective_config.json", "history.csv", "initial_particles.csv")):
            raise ValueError("Output already contains a run; use a new directory")
        if not generator.is_file():
            raise ValueError("Missing particle generator: " + str(generator))
        mpirun = shutil.which("mpirun") if args.ranks > 1 else None
        if args.ranks > 1 and not mpirun:
            raise ValueError("mpirun was not found")
    except (ValueError, OSError, KeyError, TypeError) as exc:
        command.error(str(exc))
    output.mkdir(parents=True, exist_ok=True)
    effective = output / "effective_config.json"
    atomic_json(effective, cfg)
    particles = output / "initial_particles.csv"
    # Retain the same deterministic initial particle inventory for restart validation.
    subprocess.run([sys.executable, str(generator), "--config", str(effective), "--output", str(particles)], check=True)
    if previous and digest(particles) != previous["sha256"]["particles"]:
        raise ValueError("Restart initial-particle inventory differs from the original run; preserve the original particle configuration and seed")
    values = solver_values(cfg, output, particles, restart, args.benchmark_steps)
    config_path = output / "resolved_run.cfg"
    with config_path.open("w", encoding="utf-8") as stream:
        for key, value in values.items():
            if "\n" in str(value) or "\r" in str(value):
                raise ValueError("Newlines cannot occur in resolved configuration values")
            stream.write(str(key) + "=" + str(value) + "\n")
    argv = ([mpirun, "-np", str(args.ranks)] if mpirun else []) + [str(executable), "--config", str(config_path)]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "config": cfg, "derived": meta, "ranks": args.ranks, "argv": argv,
        "restart_directory": str(restart) if restart else None,
        "parent_manifest_sha256": digest(restart.parent / "manifest.json") if restart else None,
        "sha256": {"executable": digest(executable), "driver": digest(__file__), "generator": digest(generator), "particles": digest(particles), "effective_config": digest(effective)},
    }
    atomic_json(output / "manifest.json", manifest)
    print(json.dumps(meta, indent=2), flush=True)
    print("Run directory: " + str(output), flush=True)
    # Signal handling asks all ranks to checkpoint at a synchronized time step.
    # SIGUSR1 forwarding by MPI launchers is implementation dependent, so a shared
    # STOP_REQUEST file is also written for the solver to poll on rank zero.
    requested = []
    process = None
    def request_stop(number, frame):
        requested.append(number)
        (output / "STOP_REQUEST").touch()
        print("Stop requested; waiting for a collective checkpoint.", flush=True)
    for number in (signal.SIGUSR1, signal.SIGTERM, signal.SIGINT):
        signal.signal(number, request_stop)
    start = time.monotonic()
    with (output / "solver.log").open("w", encoding="utf-8") as logfile:
        process = subprocess.Popen(argv, cwd=str(output), stdout=logfile, stderr=subprocess.STDOUT, start_new_session=True)
        rc = process.wait()
    solver_status_path = output / "status.json"
    solver_status = json.loads(solver_status_path.read_text(encoding="utf-8")) if solver_status_path.is_file() else {}
    if rc == 0 and solver_status.get("status") not in ("COMPLETED", "CHECKPOINTED"):
        print("Solver did not write a valid completion/checkpoint status.", file=sys.stderr)
        rc = 1
    status = {"exit_code": rc, "wall_seconds": time.monotonic() - start, "stop_requested": bool(requested) or (output / "STOP_REQUEST").exists(),
              "status": solver_status.get("status", "UNKNOWN") if rc == 0 else "FAILED",
              "solver_status": solver_status}
    atomic_json(output / "driver_status.json", status)
    print(json.dumps(status, indent=2), flush=True)
    print("Solver log: " + str(output / "solver.log"), flush=True)
    return rc if rc >= 0 else 128 - rc


if __name__ == "__main__":
    sys.exit(main())
