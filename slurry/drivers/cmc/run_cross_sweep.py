#!/usr/bin/env python3
"""Run independent OpenLB Couette cases; Python standard library only.

The material table is the inverse fit to Figure S2's printed black curves.
Input shear rates are physical s^-1, viscosities Pa s, and all cases use a
fixed time step selected separately for each concentration. This driver does
not run OpenLB on import, install packages, or alter the OpenLB installation.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


def cross_eta(rate, row):
    """Dynamic viscosity in Pa s; tau_s is the Cross time, not LBM tau."""
    return row["eta_inf_Pa_s"] + (row["eta0_Pa_s"] - row["eta_inf_Pa_s"]) / (
        1.0 + (row["tau_s"] * rate) ** row["m"]
    )


def logspace(start, stop, count):
    if count == 1:
        return [start]
    a, b = math.log(start), math.log(stop)
    rates = [math.exp(a + (b - a) * i / (count - 1)) for i in range(count)]
    rates[0], rates[-1] = start, stop
    return rates


def positive_float(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return result


def positive_int(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def atomic_json(path, data):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def read_materials(path):
    required = (
        "concentration_g_L", "eta0_Pa_s", "eta_inf_Pa_s", "tau_s", "m",
        "digitized_gamma_min_s_inv", "digitized_gamma_max_s_inv",
    )
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing parameter columns: " + ", ".join(sorted(missing)))
        rows = [{key: float(row[key]) for key in required} for row in reader]
    if not rows:
        raise ValueError("Empty parameter table")
    seen = set()
    for row in rows:
        c = row["concentration_g_L"]
        if any(not math.isfinite(value) for value in row.values()):
            raise ValueError(f"Non-finite parameter at concentration {c}")
        if c <= 0 or c in seen:
            raise ValueError(f"Invalid or duplicate concentration {c}")
        seen.add(c)
        if not (row["eta0_Pa_s"] > 0 and 0 <= row["eta_inf_Pa_s"] <= row["eta0_Pa_s"]):
            raise ValueError(f"Invalid viscosity parameters at concentration {c}")
        if row["tau_s"] <= 0 or row["m"] <= 0:
            raise ValueError(f"Invalid Cross parameters at concentration {c}")
        if not (0 < row["digitized_gamma_min_s_inv"] <= row["digitized_gamma_max_s_inv"]):
            raise ValueError(f"Invalid digitized shear-rate interval at concentration {c}")
    return sorted(rows, key=lambda row: row["concentration_g_L"])


def choose_materials(rows, requested):
    if requested.strip().lower() == "all":
        return rows
    values = [positive_float(value.strip()) for value in requested.split(",")]
    selected = []
    for value in values:
        matches = [row for row in rows if math.isclose(row["concentration_g_L"], value, rel_tol=1e-10)]
        if len(matches) != 1:
            raise ValueError(f"Concentration {value:g} g/L not in the parameter table")
        if matches[0] not in selected:
            selected.append(matches[0])
    return sorted(selected, key=lambda row: row["concentration_g_L"])


def save_summary(output, manifest, records, run_status):
    base = [
        "case_id", "status", "concentration_g_L", "shear_rate_s_inv",
        "shear_rate_measured_s_inv", "eta_target_Pa_s", "eta_measured_Pa_s",
        "shear_stress_Pa", "relative_error", "converged", "steps", "dt_s",
        "max_mach", "tau_min", "tau_max", "profile_relative_error",
        "density_relative_drift", "stress_method", "outside_digitized_range",
        "exit_code", "case_directory", "error",
    ]
    extra = sorted(set().union(*(record.keys() for record in records)) - set(base)) if records else []
    temp = output / "results.csv.tmp"
    with temp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=base + extra)
        writer.writeheader()
        for record in records:
            writer.writerow({key: json.dumps(value) if isinstance(value, (dict, list)) else value
                             for key, value in record.items()})
    temp.replace(output / "results.csv")
    atomic_json(output / "results.json", {
        "run_status": run_status,
        "planned_cases": len(manifest["cases"]),
        "completed_cases": len(records),
        "pass_count": sum(record.get("status") == "PASS" for record in records),
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": "manifest.json",
        "results": records,
    })
    write_svg(output / "cross_comparison.svg", manifest, records)


def write_svg(path, manifest, records):
    """Standalone scientific plot without numpy/matplotlib dependency."""
    rows = manifest["materials"]
    cases = manifest["cases"]
    palette = ["#1764ab", "#e07b18", "#228443", "#c93236", "#8658a5", "#8c564b",
               "#cf5c9c", "#6b6b6b", "#929015", "#15929a", "#123b70", "#a94e00",
               "#105629", "#821e25", "#563777", "#583225", "#882b67"]
    colors = {row["concentration_g_L"]: palette[i % len(palette)] for i, row in enumerate(rows)}
    all_rates = [case["shear_rate_s_inv"] for case in cases]
    xmin, xmax = math.log10(min(all_rates)), math.log10(max(all_rates))
    if xmax - xmin < 0.2:
        xmin, xmax = xmin - 0.5, xmax + 0.5
    curves = []
    viscosities = []
    for row in rows:
        rates = [case["shear_rate_s_inv"] for case in cases if case["concentration_g_L"] == row["concentration_g_L"]]
        gamma = logspace(min(rates), max(rates), 240)
        curve = [(rate, 1000 * cross_eta(rate, row)) for rate in gamma]
        curves.append((row, curve))
        viscosities.extend(eta for _, eta in curve)
    for record in records:
        eta = record.get("eta_measured_Pa_s")
        if isinstance(eta, (int, float)) and math.isfinite(eta) and eta > 0:
            viscosities.append(1000 * eta)
    ymin, ymax = math.log10(min(viscosities)), math.log10(max(viscosities))
    padding = max(0.12, 0.08 * (ymax - ymin))
    ymin, ymax = ymin - padding, ymax + padding
    width, height = 1040, max(640, 235 + len(rows) * 24)
    left, right, top, bottom = 95, 785, 90, height - 95
    xp = lambda value: left + (math.log10(value) - xmin) / (xmax - xmin) * (right - left)
    yp = lambda value: bottom - (math.log10(value) - ymin) / (ymax - ymin) * (bottom - top)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<g font-family="Arial, Helvetica, sans-serif" font-size="15" fill="#202020">',
             '<text x="95" y="31" font-size="22" font-weight="bold">CMC 250k: Cross law and OpenLB Couette results</text>',
             '<text x="95" y="57" font-size="14">Lines: input Cross law. Circles: viscosity from LBM momentum-flux stress.</text>']
    for axis, low, high in (("x", xmin, xmax), ("y", ymin, ymax)):
        for exponent in range(math.ceil(low), math.floor(high) + 1):
            value = 10.0 ** exponent
            label = f"10^{exponent}"
            if axis == "x":
                x = xp(value)
                parts += [f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#e5e5e5"/>',
                          f'<text x="{x:.2f}" y="{bottom + 25}" text-anchor="middle">{label}</text>']
            else:
                y = yp(value)
                parts += [f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#e5e5e5"/>',
                          f'<text x="{left - 12}" y="{y + 5:.2f}" text-anchor="end">{label}</text>']
    # A narrow y range may contain no power-of-ten ticks.
    if math.ceil(ymin) > math.floor(ymax):
        for fraction in (0.15, 0.5, 0.85):
            value = 10 ** (ymin + fraction * (ymax - ymin))
            y = yp(value)
            parts += [f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#e5e5e5"/>',
                      f'<text x="{left - 12}" y="{y + 5:.2f}" text-anchor="end">{value:.3g}</text>']
    parts += [f'<rect x="{left}" y="{top}" width="{right-left}" height="{bottom-top}" fill="none" stroke="#333"/>',
              f'<text x="{(left+right)/2}" y="{bottom + 56}" text-anchor="middle">Shear rate (s^-1)</text>',
              f'<text transform="translate(25 {(top+bottom)/2}) rotate(-90)" text-anchor="middle">Dynamic viscosity (mPa s)</text>']
    for row, curve in curves:
        color = colors[row["concentration_g_L"]]
        points = " ".join(f"{xp(rate):.2f},{yp(eta):.2f}" for rate, eta in curve)
        parts.append(f'<polyline points="{points}" stroke="{color}" stroke-width="2" fill="none"/>')
    for record in records:
        eta = record.get("eta_measured_Pa_s")
        if not isinstance(eta, (int, float)) or not math.isfinite(eta) or eta <= 0:
            continue
        color = colors[record["concentration_g_L"]]
        fill = color if record.get("status") == "PASS" else "white"
        title = html.escape(f"c={record['concentration_g_L']:g} g/L, shear={record['shear_rate_s_inv']:.6g} s^-1, {record.get('status', '')}")
        parts.append(f'<circle cx="{xp(record["shear_rate_s_inv"]):.2f}" cy="{yp(1000*eta):.2f}" r="4.3" fill="{fill}" stroke="{color}" stroke-width="1.5"><title>{title}</title></circle>')
    parts.append('<text x="813" y="106" font-weight="bold">CMC (g/L)</text>')
    for i, row in enumerate(rows):
        y, color = 132 + i * 24, colors[row["concentration_g_L"]]
        parts += [f'<line x1="813" y1="{y}" x2="842" y2="{y}" stroke="{color}" stroke-width="2"/>',
                  f'<text x="852" y="{y+5}">{row["concentration_g_L"]:g}</text>']
    count = sum(record.get("status") == "PASS" for record in records)
    parts += [f'<text x="95" y="{height-13}" font-size="13">PASS: {count}/{len(cases)} planned cases. Open circles: failed case diagnostics, if available.</text>',
              '</g></svg>']
    temp = path.with_name(path.name + ".tmp")
    temp.write_text("\n".join(parts) + "\n", encoding="utf-8")
    temp.replace(path)


def parser():
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--executable", type=Path, required=True)
    command.add_argument("--parameters", type=Path, default=Path(__file__).with_name("cmc250k_cross_parameters.csv"))
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--ranks", type=positive_int, default=1)
    command.add_argument("--concentrations", default="16", help="Comma-separated g/L values, or all; default: 16")
    command.add_argument("--points", type=positive_int, default=9, help="Log-spaced rates per concentration")
    command.add_argument("--shear-rates", help="Explicit comma-separated positive rates in s^-1; overrides --points")
    command.add_argument("--resolution", type=positive_int, default=64)
    command.add_argument("--gap", type=positive_float, default=1e-4)
    command.add_argument("--density", type=positive_float, default=1000.0)
    command.add_argument("--width-cells", type=positive_int, default=16)
    command.add_argument("--max-steps", type=positive_int,
                         help="Default: max(400000, 100*resolution^2); stops early on convergence")
    command.add_argument("--check-every", type=positive_int, default=500)
    command.add_argument("--vtk", choices=("final", "off"), default="final")
    return command


def main():
    command = parser()
    args = command.parse_args()
    if args.max_steps is None:
        args.max_steps = max(400000, 100 * args.resolution ** 2)
    executable = args.executable.expanduser().resolve()
    parameters = args.parameters.expanduser().resolve()
    output = args.output.expanduser().resolve()
    try:
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f"Executable not found or not executable: {executable}")
        if args.resolution < 16 or args.width_cells < 8:
            raise ValueError("Use --resolution >= 16 and --width-cells >= 8")
        if args.ranks > args.resolution // 4:
            raise ValueError("Use --ranks <= resolution/4 for this small grid")
        if args.check_every * 5 > args.max_steps:
            raise ValueError("Allow at least five convergence checks with --max-steps")
        mpirun = shutil.which("mpirun") if args.ranks > 1 else None
        if args.ranks > 1 and not mpirun:
            raise ValueError("mpirun was not found; load the same MPI module used to build OpenLB")
        materials = choose_materials(read_materials(parameters), args.concentrations)
        explicit_rates = sorted(set(positive_float(value.strip()) for value in args.shear_rates.split(","))) if args.shear_rates else None
        if any((output / name).exists() for name in ("manifest.json", "cases", "results.csv", "results.json")):
            raise ValueError(f"Output already contains a sweep: {output}. Choose a new output directory.")
        dx = args.gap / args.resolution
        cases = []
        for row in materials:
            c = row["concentration_g_L"]
            rates = explicit_rates or logspace(row["digitized_gamma_min_s_inv"], row["digitized_gamma_max_s_inv"], args.points)
            nu_ref = math.sqrt(row["eta0_Pa_s"] * cross_eta(1000.0, row)) / args.density
            dt = 0.1 * dx * dx / nu_ref
            for i, rate in enumerate(rates):
                mach = math.sqrt(3) * rate * args.gap * dt / dx
                if mach >= 0.1:
                    raise ValueError(f"c={c:g}, shear={rate:g}: estimated wall Mach={mach:.4g} > 0.1; reduce --gap or increase --resolution")
                cases.append({
                    "case_id": f"c{c:g}_g{i:03d}", "concentration_g_L": c,
                    "shear_rate_s_inv": rate, "dt_s": dt,
                    "eta_target_Pa_s": cross_eta(rate, row),
                    "estimated_wall_mach": mach,
                    "outside_digitized_range": not (row["digitized_gamma_min_s_inv"] <= rate <= row["digitized_gamma_max_s_inv"]),
                    "material": row,
                })
        manifest = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "executable": str(executable), "executable_sha256": digest(executable),
            "parameters": str(parameters), "parameters_sha256": digest(parameters),
            "driver_sha256": digest(Path(__file__).resolve()),
            "settings": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "dx_m": dx, "dt_rule": "0.1 * dx^2 * density / sqrt(eta0 * eta(1000 s^-1))",
            "materials": materials, "cases": cases,
            "notes": ["eta_inf=0 is a fitted boundary value; no claim of zero physical infinite-shear viscosity.",
                      "Default rates lie inside each concentration's digitized source interval.",
                      "Cases start independently; no shear-history or concentration interpolation is imposed."],
        }
    except (ValueError, OSError, argparse.ArgumentTypeError) as exc:
        command.error(str(exc))
    output.mkdir(parents=True, exist_ok=True)
    (output / "cases").mkdir()
    atomic_json(output / "manifest.json", manifest)
    records = []
    save_summary(output, manifest, records, "RUNNING")
    print(f"Running {len(cases)} independent cases in {output}", flush=True)
    interrupted = False
    try:
        for i, case in enumerate(cases, 1):
            row = case["material"]
            directory = output / "cases" / case["case_id"]
            directory.mkdir()
            argv = ([mpirun, "--np", str(args.ranks)] if args.ranks > 1 else []) + [str(executable)]
            options = {
                "eta0": row["eta0_Pa_s"], "eta-inf": row["eta_inf_Pa_s"],
                "lambda": row["tau_s"], "m": row["m"], "concentration": row["concentration_g_L"],
                "shear-rate": case["shear_rate_s_inv"], "resolution": args.resolution,
                "gap": args.gap, "density": args.density, "dt": case["dt_s"],
                "width-cells": args.width_cells, "max-steps": args.max_steps,
                "check-every": args.check_every, "vtk-every": 0,
            }
            for key, value in options.items():
                argv.extend(["--" + key, str(value)])
            if args.vtk == "off":
                argv.append("--no-vtk")
            atomic_json(directory / "command.json", {"argv": argv, "cwd": str(directory), "input": case})
            record = {key: value for key, value in case.items() if key != "material"}
            record["case_directory"] = str(directory.relative_to(output))
            print(f"[{i}/{len(cases)}] c={row['concentration_g_L']:g} g/L, shear={case['shear_rate_s_inv']:.6g} s^-1", flush=True)
            try:
                with (directory / "run.log").open("w", encoding="utf-8") as logfile:
                    process = subprocess.run(argv, cwd=directory, stdout=logfile, stderr=subprocess.STDOUT, check=False)
                record["exit_code"] = process.returncode
                result_path = directory / "result.json"
                if result_path.is_file():
                    result = json.loads(result_path.read_text(encoding="utf-8"), parse_constant=lambda value: None)
                    if not isinstance(result, dict):
                        raise ValueError("result.json must contain one JSON object")
                    # Preserve actual measured values; prevent result metadata from changing the requested case.
                    result.pop("case_id", None)
                    result.pop("case_directory", None)
                    record.update(result)
                    for name in ("concentration_g_L", "shear_rate_s_inv", "dt_s"):
                        if not isinstance(record.get(name), (int, float)) or not math.isclose(record[name], case[name], rel_tol=1e-7, abs_tol=1e-15):
                            raise ValueError(f"Result {name} does not match the requested case")
                    eta = record.get("eta_measured_Pa_s")
                    if not isinstance(eta, (int, float)) or not math.isfinite(eta) or eta <= 0:
                        raise ValueError("Missing or invalid eta_measured_Pa_s")
                    record["eta_target_Pa_s"] = case["eta_target_Pa_s"]
                    record["relative_error"] = abs(eta / case["eta_target_Pa_s"] - 1)
                    if record.get("status") != "PASS" or record.get("converged") is not True:
                        record["status"] = "FAIL"
                else:
                    raise ValueError("result.json was not written; inspect run.log")
                if process.returncode != 0:
                    record["status"] = "FAIL"
                    record["error"] = f"OpenLB/launcher exited with code {process.returncode}"
            except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
                record["status"] = "ERROR"
                record["error"] = str(exc)
            records.append(record)
            save_summary(output, manifest, records, "RUNNING")
            print(f"  {record.get('status', 'ERROR')}: {directory / 'run.log'}", flush=True)
    except KeyboardInterrupt:
        interrupted = True
        print("Interrupted; completed case results are preserved.", file=sys.stderr, flush=True)
    all_passed = len(records) == len(cases) and all(record.get("status") == "PASS" for record in records)
    status = "INTERRUPTED" if interrupted else ("PASS" if all_passed else "FAIL")
    save_summary(output, manifest, records, status)
    print(f"Sweep {status}: {sum(record.get('status') == 'PASS' for record in records)}/{len(cases)} PASS", flush=True)
    print(f"Summary: {output / 'results.csv'}\nPlot: {output / 'cross_comparison.svg'}", flush=True)
    return 130 if interrupted else (0 if all_passed else 1)


if __name__ == "__main__":
    sys.exit(main())
