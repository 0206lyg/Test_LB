#!/usr/bin/env python3
"""Summarize particle-solver failure traces using only Python's standard library.

Run from a graphite result directory: python3 summarize_particle_solver.py
Or point at a result directory/trace: python3 summarize_particle_solver.py PATH
No simulation data or diagnostics are modified.
"""
import argparse
from collections import OrderedDict, deque
import csv
import json
import math
from pathlib import Path
import sys


ATTEMPT_FIELDS = ("outer_time_s", "outer_dt_s", "subdivision_count", "substep_index")
TERMINAL_ACCEPTED = {"outer_accepted", "outer_success", "accepted_outer", "recovered"}
TERMINAL_FAILED = {"outer_failed", "outer_failure", "failed_outer", "failed_outer_step", "fatal"}


def number(value):
    """Return a finite number or None; an interrupted last CSV line is allowed."""
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) else None


def display(value):
    value = number(value)
    return "-" if value is None else format(value, ".6g")


def event_of(row):
    return str(row.get("event") or row.get("status") or "").lower()


def discover(paths, recursive=False):
    files = set()
    for path in paths:
        path = path.expanduser().resolve()
        if path.is_file():
            if path.suffix.lower() != ".csv":
                raise ValueError("Expected a CSV trace: " + str(path))
            if path.name.endswith("_outcomes.csv"):
                path = path.with_name(path.name[:-len("_outcomes.csv")] + "_trace.csv")
                if not path.is_file():
                    raise ValueError("No iteration trace accompanies the outcomes file")
            files.add(path)
        elif path.is_dir():
            iterator = path.rglob("particle_solver*.csv") if recursive else path.glob("particle_solver*.csv")
            files.update(item.resolve() for item in iterator
                         if item.is_file() and not item.name.endswith("_outcomes.csv"))
        else:
            raise ValueError("Path does not exist: " + str(path))
    return sorted(files)


def read_trace(path, tail=10):
    """Keep bounded tails for up to 64 recent attempts, including retry outcomes."""
    attempts = OrderedDict()
    last_event = None
    row_count = 0
    ignored_rows = 0
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        header = set(reader.fieldnames or ())
        required = set(ATTEMPT_FIELDS) | {"newton_iteration", "force_ratio", "torque_ratio"}
        if not required <= header:
            missing = ", ".join(sorted(required - header))
            raise ValueError(str(path) + ": missing diagnostic columns: " + missing)
        for row in reader:
            if None in row or any(row.get(key) in (None, "") for key in ATTEMPT_FIELDS):
                ignored_rows += 1
                continue
            row_count += 1
            event = event_of(row)
            if event in TERMINAL_ACCEPTED | TERMINAL_FAILED:
                last_event = dict(row)
                continue
            key = tuple(row[key] for key in ATTEMPT_FIELDS)
            if key not in attempts:
                attempts[key] = {"first": None, "last": None, "footer": None,
                                 "metadata": dict(row), "tail": deque(maxlen=tail), "rows": 0}
                if len(attempts) > 64:
                    attempts.popitem(last=False)
            attempt = attempts[key]
            if event == "failed_attempt":
                attempt["footer"] = dict(row)
                continue
            if attempt["first"] is None:
                attempt["first"] = dict(row)
            attempt["last"] = dict(row)
            attempt["tail"].append(dict(row))
            attempt["rows"] += 1
    prefix = path.stem.split("_iterations")[0].split("_trace")[0]
    outcome_path = path.with_name(prefix + "_outcomes.csv")
    artifacts = [str(item.resolve()) for item in sorted(path.parent.glob(prefix + "*"))
                 if item.is_file() and item != path and item.suffix.lower() != ".csv"]
    recent = list(attempts.values())
    if not recent:
        return {"file": str(path), "row_count": row_count,
                "ignored_rows": ignored_rows, "last_attempt": None,
                "terminal_event": last_event, "retry_counts": [], "artifacts": artifacts}
    last = recent[-1]
    outer_key = tuple(last["metadata"][key] for key in ATTEMPT_FIELDS[:2])
    retries = []
    for attempt in recent:
        if tuple(attempt["metadata"][key] for key in ATTEMPT_FIELDS[:2]) == outer_key:
            count = attempt["metadata"]["subdivision_count"]
            if count not in retries:
                retries.append(count)
    # An old terminal outcome does not describe a newer attempted outer step.
    if last_event and tuple(last_event[key] for key in ATTEMPT_FIELDS[:2]) != outer_key:
        last_event = None
    if outcome_path.is_file():
        with outcome_path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if tuple(row.get(key) for key in ATTEMPT_FIELDS[:2]) == outer_key:
                    if event_of(row) in TERMINAL_ACCEPTED | TERMINAL_FAILED:
                        last_event = dict(row)
        artifacts.append(str(outcome_path.resolve()))
    return {"file": str(path), "row_count": row_count,
            "ignored_rows": ignored_rows,
            "last_attempt": {"first": last["first"], "last": last["last"],
                             "metadata": last["metadata"], "footer": last["footer"],
                             "tail": list(last["tail"]), "rows": last["rows"]},
            "terminal_event": last_event, "retry_counts": retries, "artifacts": artifacts}


def render(report):
    lines = ["Trace: " + report["file"]]
    attempt = report["last_attempt"]
    if attempt is None:
        terminal = report["terminal_event"]
        if terminal and event_of(terminal) in TERMINAL_FAILED:
            lines.append("Outcome: outer step failed before any iteration trace was recorded")
            lines.append("Recorded reason: " + str(terminal.get("message") or "no message"))
        else:
            lines.append("No failed-attempt iterations recorded.")
        if report["artifacts"]:
            lines.append("Related diagnostics:")
            lines.extend("  " + path for path in report["artifacts"])
        return "\n".join(lines)
    first, last = attempt["first"], attempt["last"]
    metadata = attempt["metadata"]
    terminal = report["terminal_event"]
    if terminal and event_of(terminal) in TERMINAL_FAILED:
        outcome = "outer step failed"
    elif terminal and event_of(terminal) in TERMINAL_ACCEPTED:
        outcome = "retry recovered; outer step accepted"
    else:
        outcome = "failed attempt recorded; final outer-step outcome is not recorded here"
    lines.append("Outcome: " + outcome)
    lines.append("Outer time=" + display(metadata["outer_time_s"]) + " s; dt="
                 + display(metadata["outer_dt_s"]) + " s")
    lines.append("Recorded subdivision attempts: " + ", ".join(report["retry_counts"]))
    lines.append("Last attempt: count=" + str(metadata["subdivision_count"])
                 + "; substep_index=" + str(metadata["substep_index"])
                 + "; subdt=" + display(metadata.get("subdt_s")) + " s")
    reason_row = attempt["footer"] or last or metadata
    reason = str((terminal or {}).get("message") or reason_row.get("message") or "").strip()
    code = str(reason_row.get("snes_reason") or "").strip()
    lines.append("Recorded reason: " + (reason or "no message")
                 + (" [SNES reason " + code + "]" if code else ""))
    if first is None or last is None:
        lines.append("No nonlinear iteration values were recorded for this attempt.")
        if report["artifacts"]:
            lines.append("Related diagnostics:")
            lines.extend("  " + path for path in report["artifacts"])
        return "\n".join(lines)
    lines.append("Force ratio: " + display(first.get("force_ratio")) + " -> "
                 + display(last.get("force_ratio")) + "; torque ratio: "
                 + display(first.get("torque_ratio")) + " -> " + display(last.get("torque_ratio")))
    lines.append("Gap violation=" + display(last.get("gap_violation_m"))
                 + " m; complementarity ratio=" + display(last.get("complementarity_ratio")))
    has_contacts = "active_contacts" in last
    if has_contacts:
        lines.append("Contacts: active=" + display(last.get("active_contacts"))
                     + "; sliding=" + display(last.get("sliding_contacts"))
                     + "; rolling=" + display(last.get("rolling_contacts"))
                     + "; rejected domain evaluations=" + display(last.get("domain_errors"))
                     + "; KSP reason=" + display(last.get("ksp_reason")))
    if "friction_branch_attempts" in last:
        lines.append("Friction boundary corrections: attempts="
                     + display(last.get("friction_branch_attempts"))
                     + "; descent directions=" + display(last.get("friction_branch_corrections")))
    if "contact_state_updates" in last:
        lines.append("Contact state updates: " + display(last.get("contact_state_updates"))
                     + "; activated=" + display(last.get("contact_activations"))
                     + "; released=" + display(last.get("contact_releases")))
    has_ngmres = "nonlinear_iteration" in last
    if has_ngmres:
        lines.append("Nonlinear acceleration: NGMRES iteration="
                     + display(last.get("nonlinear_iteration"))
                     + "; actual Newton attempts=" + display(last.get("newton_iteration"))
                     + "; cumulative Krylov iterations=" + display(last.get("total_ksp_iterations"))
                     + "; residual evaluations=" + display(last.get("residual_evaluations")))
        lines.append("Newton candidate SNES reason=" + display(last.get("npc_snes_reason"))
                     + "; +5/-5 can mark one-step completion, not physical convergence.")
        lines.append("Step fraction below is the Newton candidate fraction, before NGMRES selection.")
    lines.append("Final iteration rows (force/torque/complementarity ratios pass at <=1):")
    lines.append("  iter       force      torque     contact   step_frac    KSP_residual"
                 + ("   active  added removed expand" if has_contacts else ""))
    for row in attempt["tail"]:
        values = [row.get(key) for key in ("newton_iteration", "force_ratio", "torque_ratio",
                                         "complementarity_ratio",
                                         "newton_step_fraction" if has_ngmres else "step_fraction",
                                         "ksp_residual_norm")]
        line = "  " + " ".join(display(value).rjust(width)
                               for value, width in zip(values, (4, 11, 11, 11, 11, 15)))
        if has_contacts:
            line += " " + " ".join(display(row.get(key)).rjust(6)
                                    for key in ("active_contacts", "activated_contacts", "released_contacts", "candidate_expansion"))
        lines.append(line)
    if report["ignored_rows"]:
        lines.append("Incomplete/malformed rows skipped: " + str(report["ignored_rows"]))
    if report["artifacts"]:
        lines.append("Related diagnostics:")
        lines.extend("  " + path for path in report["artifacts"])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, default=[Path.cwd()],
                        help="Result directories or diagnostic CSV files (default: current directory)")
    parser.add_argument("--recursive", action="store_true", help="Search nested run directories")
    parser.add_argument("--tail", type=int, default=10, help="Iteration rows to show for the last attempt")
    parser.add_argument("--json", action="store_true", help="Print machine-readable summaries")
    args = parser.parse_args(argv)
    if args.tail < 1:
        parser.error("--tail must be positive")
    try:
        files = discover(args.paths or [Path.cwd()], args.recursive)
        if not files:
            print("No particle_solver*.csv traces found in the selected directory.")
            return 0
        reports = [read_trace(path, args.tail) for path in files]
    except (OSError, ValueError, csv.Error) as error:
        print("DIAGNOSTIC READ FAILED: " + str(error), file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(reports, indent=2, allow_nan=False))
    else:
        print("\n\n".join(render(report) for report in reports))
    return 0


if __name__ == "__main__":
    sys.exit(main())
