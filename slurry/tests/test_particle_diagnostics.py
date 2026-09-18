#!/usr/bin/env python3
"""Regression checks for recorded failed retries, recovery and interrupted CSV."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / "tools/summarize_particle_solver.py"
SPEC = importlib.util.spec_from_file_location("particle_summary", MODULE)
summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary)

FIELDS = ["outer_time_s", "outer_dt_s", "subdivision_count", "substep_index",
          "subdt_s", "newton_iteration", "force_ratio", "torque_ratio",
          "gap_violation_m", "complementarity_ratio", "residual_norm",
          "step_fraction", "ksp_iterations", "ksp_residual_norm", "snes_reason",
          "message", "event"]


def row(**kwargs):
    result = {name: "" for name in FIELDS}
    result.update(outer_time_s="0.033", outer_dt_s="1e-6", subdivision_count="1",
                  substep_index="0", subdt_s="1e-6", newton_iteration="0",
                  force_ratio="100", torque_ratio="10", gap_violation_m="0",
                  complementarity_ratio="0", residual_norm="100",
                  step_fraction="1", ksp_iterations="4", ksp_residual_norm="1e-8",
                  snes_reason="0", message="", event="iteration")
    result.update({key: str(value) for key, value in kwargs.items()})
    return result


class DiagnosticSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.path = self.directory / "particle_solver_rank0_iterations.csv"

    def write(self, rows):
        with self.path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)

    def test_failed_retry_is_not_reported_as_aborted_run(self):
        self.write([row(), row(newton_iteration=60, force_ratio=40, torque_ratio=25,
                               snes_reason=-5, message="iteration limit"),
                    row(subdivision_count=2, event="outer_accepted", snes_reason=2,
                        message="accepted after retry")])
        result = summary.read_trace(self.path)
        self.assertIn("retry recovered; outer step accepted", summary.render(result))
        self.assertEqual(result["last_attempt"]["last"]["force_ratio"], "40")

    def test_aborted_outer_step_shows_terminal_reason_and_last_retry(self):
        self.write([row(), row(subdivision_count=2, force_ratio=70),
                    row(subdivision_count=256, substep_index=13, subdt_s="3.90625e-9"),
                    row(subdivision_count=256, substep_index=13, newton_iteration=60,
                        force_ratio=45.63, torque_ratio=29.36, snes_reason=-5,
                        message="nonlinear iterations exhausted"),
                    row(subdivision_count=256, substep_index=13, event="outer_failed",
                        message="all subdivisions exhausted")])
        result = summary.read_trace(self.path, tail=1)
        output = summary.render(result)
        self.assertEqual(result["retry_counts"], ["1", "2", "256"])
        self.assertEqual(len(result["last_attempt"]["tail"]), 1)
        self.assertIn("outer step failed", output)
        self.assertIn("all subdivisions exhausted", output)
        self.assertIn("45.63", output)
        self.assertIn("substep_index=13", output)

    def test_prior_recovery_does_not_hide_newer_unfinished_failure(self):
        self.write([row(event="outer_accepted"), row(outer_time_s="0.034", snes_reason=-6)])
        result = summary.read_trace(self.path)
        self.assertIsNone(result["terminal_event"])
        self.assertIn("final outer-step outcome is not recorded", summary.render(result))

    def test_quoted_message_and_truncated_record_are_readable(self):
        self.write([row(newton_iteration=4, message="line search, residual \"stalled\"", snes_reason=-6)])
        with self.path.open("a", encoding="utf-8") as stream:
            stream.write("0.033,1e-6\n")
        result = summary.read_trace(self.path)
        self.assertEqual(result["ignored_rows"], 1)
        self.assertIn('line search, residual "stalled"', summary.render(result))

    def test_related_snapshot_and_non_diagnostic_history(self):
        self.write([row()])
        (self.directory / "history.csv").write_text("time,stress\n0,0\n")
        snapshot = self.directory / "particle_solver_rank0_failure.txt"
        snapshot.write_text("failing state\n")
        self.assertEqual(summary.discover([self.directory]), [self.path.resolve()])
        self.assertEqual(summary.read_trace(self.path)["artifacts"], [str(snapshot.resolve())])

    def test_empty_header_is_not_interpreted_as_failure(self):
        self.write([])
        self.assertIn("No failed-attempt iterations", summary.render(summary.read_trace(self.path)))

    def test_pre_iteration_failure_retains_reason(self):
        self.write([row(event="outer_failed", message="initial state violates gap")])
        output = summary.render(summary.read_trace(self.path))
        self.assertIn("outer step failed before any iteration", output)
        self.assertIn("initial state violates gap", output)

    def write_outcome(self, status):
        path = self.directory / "particle_solver_rank0_outcomes.csv"
        path.write_text("outer_time_s,outer_dt_s,subdivision_count,status\n"
                        + "0.033,1e-6,2," + status + "\n")
        return path

    def test_actual_backend_footer_does_not_replace_last_residual(self):
        self.write([row(event="failed_iteration", force_ratio=45.63, torque_ratio=29.36),
                    row(event="failed_attempt", newton_iteration=-1, force_ratio=0,
                        torque_ratio=0, message="budget exhausted", snes_reason=-5)])
        outcome = self.write_outcome("failed_outer_step")
        result = summary.read_trace(self.path)
        self.assertIn("45.63 -> 45.63", summary.render(result))
        self.assertIn("outer step failed", summary.render(result))
        self.assertIn("budget exhausted", summary.render(result))
        self.assertEqual(summary.discover([self.directory]), [self.path.resolve()])
        self.assertIn(str(outcome.resolve()), result["artifacts"])

    def test_actual_backend_recovered_outcome(self):
        self.write([row(event="failed_iteration"), row(event="failed_attempt", newton_iteration=-1)])
        self.write_outcome("recovered")
        self.assertIn("retry recovered; outer step accepted", summary.render(summary.read_trace(self.path)))

    def test_actual_backend_footer_only_does_not_invent_zero_residuals(self):
        self.write([row(event="failed_attempt", newton_iteration=-1, force_ratio=0,
                        torque_ratio=0, message="initial gap invalid")])
        self.write_outcome("failed_outer_step")
        output = summary.render(summary.read_trace(self.path))
        self.assertIn("initial gap invalid", output)
        self.assertIn("No nonlinear iteration values", output)
        self.assertNotIn("Force ratio:", output)


if __name__ == "__main__":
    unittest.main()
