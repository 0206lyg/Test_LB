#!/usr/bin/env python3
"""Summarize measured viscosity history and benchmark runtime.

No packages are required for summary.json. If matplotlib is available, also
write viscosity_time.png. The averaging interval must be inspected physically;
this script never labels a transient trace as a converged steady viscosity.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics


def read_history(path):
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {"step", "time_s", "strain", "eta_bulk_Pa_s", "eta_relative", "wall_seconds"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("history.csv is missing required fields: " + ", ".join(sorted(required - set(reader.fieldnames or []))))
        rows = []
        for row in reader:
            values = {name: float(value) for name, value in row.items() if name and value is not None}
            if any(not math.isfinite(values[name]) for name in required):
                raise ValueError("Non-finite measured history value")
            if rows and values["step"] <= rows[-1]["step"]:
                raise ValueError("History steps must increase strictly")
            rows.append(values)
    if not rows:
        raise ValueError("No measurement rows in history.csv")
    return rows


def mean_over_time(rows, field):
    if len(rows) == 1:
        return rows[0][field]
    elapsed = rows[-1]["time_s"] - rows[0]["time_s"]
    if elapsed <= 0:
        raise ValueError("Selected averaging interval has no duration")
    integral = sum(0.5 * (a[field] + b[field]) * (b["time_s"] - a["time_s"]) for a, b in zip(rows, rows[1:]))
    return integral / elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--start-strain", type=float, default=0.0, help="Select an explicit averaging interval; default includes startup")
    parser.add_argument("--end-strain", type=float)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    run = args.run_directory.expanduser().resolve()
    try:
        history = read_history(run / "history.csv")
        selected = [row for row in history if row["strain"] >= args.start_strain and (args.end_strain is None or row["strain"] <= args.end_strain)]
        if not selected:
            raise ValueError("No samples in the requested strain interval")
        manifest_path = run / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
        intervals = [(b["step"] - a["step"]) / (b["wall_seconds"] - a["wall_seconds"]) for a, b in zip(history, history[1:]) if b["wall_seconds"] > a["wall_seconds"]]
        speed = statistics.median(intervals) if intervals else None
        target = manifest.get("derived", {}).get("steps_to_requested_end_strain")
        summary = {
            "run_directory": str(run), "number_of_samples": len(history), "selected_samples": len(selected),
            "selected_strain_interval": [selected[0]["strain"], selected[-1]["strain"]],
            "selected_time_interval_s": [selected[0]["time_s"], selected[-1]["time_s"]],
            "time_weighted_mean_eta_bulk_Pa_s": mean_over_time(selected, "eta_bulk_Pa_s"),
            "time_weighted_mean_relative_eta": mean_over_time(selected, "eta_relative"),
            "eta_bulk_min_Pa_s": min(row["eta_bulk_Pa_s"] for row in selected),
            "eta_bulk_max_Pa_s": max(row["eta_bulk_Pa_s"] for row in selected),
            "median_observed_steps_per_second": speed,
            "estimated_remaining_wall_seconds": max(0, target - history[-1]["step"]) / speed if speed and target else None,

            "notes": ["Mean describes only the stated sampled interval; default includes startup.",
                      "Runtime extrapolation uses measured intervals of this run and is not a hardware guarantee.",
                      "No experimental viscosity or shear-thinning fit is imposed."],
        }
        for name in ("stress_total_Pa", "stress_fluid_Pa", "stress_surface_Pa",
                     "stress_pair_attractive_Pa", "stress_pair_repulsive_Pa",
                     "stress_contact_normal_Pa", "stress_contact_tangential_Pa",
                     "stress_lubrication_Pa", "stress_acceleration_Pa",
                     "stress_fluid_reynolds_Pa", "stress_particle_reynolds_Pa",
                     "stress_noninertial_Pa", "stress_inertial_Pa"):
            if all(name in row for row in selected):
                summary["time_weighted_mean_" + name] = mean_over_time(selected, name)
        for name in ("max_mach", "max_fluid_mach", "density_drift", "fluid_density_drift",
                     "particle_substeps", "newton_iterations", "force_residual_ratio",
                     "torque_residual_ratio", "contact_gap_violation_m"):
            if all(name in row for row in history):
                summary["max_abs_" + name] = max(abs(row[name]) for row in history)
        (run / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed; history.csv and summary.json are available.")
        else:
            fig, axes = plt.subplots(2, 1, figsize=(7, 6), constrained_layout=True)
            axes[0].plot([r["time_s"] for r in history], [r["eta_bulk_Pa_s"] for r in history], linewidth=1.2)
            axes[0].set(xlabel="Time (s)", ylabel="Apparent viscosity (Pa s)")
            axes[1].plot([r["strain"] for r in history], [r["eta_relative"] for r in history], linewidth=1.2)
            axes[1].axhline(1, color="0.5", linewidth=0.8, linestyle="--")
            axes[1].set(xlabel="Accumulated strain", ylabel="Relative apparent viscosity")
            for ax in axes:
                ax.grid(alpha=0.2)
            fig.savefig(run / "viscosity_time.png", dpi=180)
            plt.close(fig)
            print("Plot: " + str(run / "viscosity_time.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
