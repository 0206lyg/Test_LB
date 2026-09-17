#!/usr/bin/env python3
"""Optional PNG/PDF rendering of a completed or partial OpenLB sweep.

Usage: python3 plot_cross_results.py /path/to/run_directory
Requires matplotlib on the computer where this plotting script is run.
The cluster runner always writes a standalone SVG without matplotlib.
"""
import argparse
import json
import math
from pathlib import Path

from run_cross_sweep import cross_eta, logspace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--output-prefix", type=Path, help="Default: RUN_DIRECTORY/cross_comparison")
    args = parser.parse_args()
    directory = args.run_directory.expanduser().resolve()
    prefix = args.output_prefix.expanduser().resolve() if args.output_prefix else directory / "cross_comparison"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        parser.error("matplotlib is required only for PNG/PDF export. Open cross_comparison.svg directly, or run this script on a computer with matplotlib.")
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    summary = json.loads((directory / "results.json").read_text(encoding="utf-8"))
    records = summary["results"]
    materials = manifest["materials"]
    fig, ax = plt.subplots(figsize=(10.2, max(5.8, 2 + 0.23 * len(materials))))
    cmap = plt.get_cmap("turbo")
    handles = []
    for i, row in enumerate(materials):
        concentration = row["concentration_g_L"]
        color = cmap(0.08 + 0.84 * i / max(1, len(materials) - 1))
        rates = [case["shear_rate_s_inv"] for case in manifest["cases"] if case["concentration_g_L"] == concentration]
        gamma = logspace(min(rates), max(rates), 300)
        ax.plot(gamma, [1000 * cross_eta(rate, row) for rate in gamma], color=color, lw=1.8)
        handles.append(Line2D([], [], color=color, lw=1.8, label=f"{concentration:g}"))
        for passed in (True, False):
            points = [record for record in records
                      if record["concentration_g_L"] == concentration
                      and (record.get("status") == "PASS") == passed
                      and isinstance(record.get("eta_measured_Pa_s"), (int, float))
                      and math.isfinite(record["eta_measured_Pa_s"])
                      and record["eta_measured_Pa_s"] > 0]
            ax.scatter([record["shear_rate_s_inv"] for record in points],
                       [1000 * record["eta_measured_Pa_s"] for record in points],
                       s=33, facecolors=color if passed else "white", edgecolors=color,
                       linewidths=1.0, zorder=3)
    ax.set(xscale="log", yscale="log", xlabel=r"Shear rate $\dot\gamma$ (s$^{-1}$)",
           ylabel=r"Dynamic viscosity $\eta$ (mPa$\,\cdot\,$s)")
    ax.grid(which="major", alpha=0.2)
    legend = ax.legend(handles=handles, title="CMC (g/L)", loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.add_artist(legend)
    ax.legend(handles=[Line2D([], [], color="#444", lw=1.8, label="Input Cross law"),
                       Line2D([], [], color="#444", marker="o", ls="", label="LBM stress: PASS"),
                       Line2D([], [], color="#444", marker="o", markerfacecolor="white", ls="", label="LBM stress: failed")],
              loc="lower left", bbox_to_anchor=(1.01, 0), frameon=False, fontsize=9)
    fig.suptitle("CMC 250k: Cross law and OpenLB Couette results", y=0.985)
    fig.text(0.10, 0.015,
             f"Sweep {summary['run_status']}; {summary['pass_count']}/{summary['planned_cases']} PASS. "
             "LBM viscosity is inferred from momentum-flux stress.", fontsize=9)
    fig.subplots_adjust(left=0.10, right=0.73, bottom=0.13, top=0.91)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        path = prefix.with_suffix("." + extension)
        fig.savefig(path, dpi=220)
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
