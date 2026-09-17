#!/usr/bin/env python3
"""Reconstruct the CMC 250k black Cross curves in SI Figure S2.

python3 fit_cmc250k.py
python3 fit_cmc250k.py --pdf path/to/nn6c10201_si_001.pdf --outdir reproduced

Without --pdf, the adjacent digitized CSV is used. Dependencies: numpy,
scipy, matplotlib; PyMuPDF is needed only for extraction from the source PDF.
All fitted viscosity values use mPa s; Pa s values are also exported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares, linear_sum_assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


# These are source-PDF drawing indices, verified against the legend RGB and
# experimental marker groups. They are NOT assigned by vertical ordering.
PATHS = {
    0.6: [944], 1.5: [945], 3: [946], 4: [947], 5: [948], 6: [949],
    7: [958], 8: [959], 9: [950], 10: [951], 11: [962], 12: [952],
    15: [953], 16: [954], 20: [955], 31: [956, 957], 42: [960, 961],
}
CONCENTRATIONS = list(PATHS)
PANEL_CLIP = [228, 82, 364, 237]
PREFIX = "cmc250k"


def cross_viscosity(gamma_dot, eta0, eta_inf, tau, m):
    """Cross viscosity; gamma_dot in s^-1, tau in s, eta values in same units.

    Physical constraints: eta0 >= eta_inf >= 0, tau > 0, m > 0.
    The output is in the same viscosity unit as eta0 and eta_inf.
    """
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    if np.any(~np.isfinite(gamma_dot)) or np.any(gamma_dot < 0):
        raise ValueError("Shear-rate magnitude must be finite and nonnegative")
    if not (eta0 >= eta_inf >= 0 and tau > 0 and m > 0):
        raise ValueError("Invalid Cross parameters")
    return eta_inf + (eta0 - eta_inf) / (1 + (tau * gamma_dot) ** m)


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def extract_pdf(pdf, outdir):
    import fitz
    doc = fitz.open(pdf)
    page = doc[2]
    drawings = page.get_drawings()
    if "250k" not in page.get_text() or len(drawings) != 1918:
        raise ValueError("This extractor is specific to the supplied SI Figure S2 PDF")
    # The top frame lies above the 1e5 tick. Calibrate from SIX MAJOR TICKS.
    x_ticks = np.array([drawings[i]["items"][0][1].x for i in range(41, 47)])
    y_ticks = np.array([drawings[i]["items"][0][1].y for i in range(184, 190)])
    xcoef = np.polyfit(np.arange(-2, 4), x_ticks, 1)
    ycoef = np.polyfit(np.arange(6), y_ticks, 1)
    if not (19 < xcoef[0] < 20 and -24 < ycoef[0] < -22):
        raise ValueError("Unexpected axis calibration")
    curves = []
    for concentration, indices in PATHS.items():
        points = []
        for index in indices:
            drawing = drawings[index]
            if drawing["color"] != (0.0, 0.0, 0.0):
                raise ValueError(f"Path {index} is not a black curve")
            for item in drawing["items"]:
                if item[0] != "l":
                    raise ValueError("Unexpected curve segment type")
                points.extend([tuple(item[1]), tuple(item[2])])
        points = np.asarray(points)
        xp = np.unique(points[:, 0])
        # Median deals with repeated vertices and vertical quantization steps.
        yp = np.array([np.median(points[points[:, 0] == x, 1]) for x in xp])
        curves.append({
            "c": concentration, "x_pdf": xp, "y_pdf": yp,
            "gamma": 10 ** ((xp - xcoef[1]) / xcoef[0]),
            "eta_mPa_s": 10 ** ((yp - ycoef[1]) / ycoef[0]),
        })

    # Independently validate identities: exact legend color -> marker group
    # -> closest printed black curve. Concentrations 7-12 are not height-sorted.
    costs, marker_counts = [], []
    for concentration in CONCENTRATIONS:
        legend_index = 1866 + list(reversed(CONCENTRATIONS)).index(concentration)
        color = drawings[legend_index]["fill"]
        markers = []
        for drawing in drawings[:944]:
            if drawing["fill"] != color or len(drawing["items"]) != 4:
                continue
            if not all(item[0] == "c" for item in drawing["items"]):
                continue
            rect = drawing["rect"]
            x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            if x_ticks[0] - .1 <= x <= x_ticks[-1] + .1 and 87 < y < 209:
                markers.append((x, y))
        markers = np.asarray(markers)
        if len(markers) < 3:
            raise ValueError(f"Cannot verify marker identity for {concentration} g/L")
        marker_counts.append(len(markers))
        row = []
        for curve in curves:
            keep = ((markers[:, 0] >= curve["x_pdf"][0] - .1)
                    & (markers[:, 0] <= curve["x_pdf"][-1] + .1))
            if keep.sum() < 3:
                row.append(1e6)
            else:
                yline = np.interp(markers[keep, 0], curve["x_pdf"], curve["y_pdf"])
                row.append(float(np.median(np.abs(markers[keep, 1] - yline))))
        costs.append(row)
    row, col = linear_sum_assignment(costs)
    if not np.array_equal(row, col):
        raise ValueError("Legend-marker matching disagrees with concentration labels")

    metadata = {
        "source": "Gwag et al., Supporting Information, Figure S2, CMC 250k middle panel",
        "source_pdf_sha256": hashlib.sha256(Path(pdf).read_bytes()).hexdigest(),
        "pdf_page_one_based": 3, "printed_page": "S3", "pdf_drawing_count": len(drawings),
        "method": "PDF vector black-line extraction, not raster clicks or raw experimental data",
        "x_major_ticks_pdf_pt": x_ticks.tolist(),
        "x_major_ticks_log10_shear_rate": list(range(-2, 4)),
        "y_major_ticks_pdf_pt": y_ticks.tolist(),
        "y_major_ticks_log10_viscosity_mPa_s": list(range(6)),
        "x_pdf_equals_a_log10_gamma_plus_b": xcoef.tolist(),
        "y_pdf_equals_a_log10_eta_mPa_s_plus_b": ycoef.tolist(),
        "panel_clip_pdf_pt": PANEL_CLIP,
        "approx_y_graphical_quantization_pt": .10316,
        "approx_y_graphical_quantization_decades": .10316 / abs(float(ycoef[0])),
        "curve_identification": [
            {"concentration_g_L": c, "pdf_drawing_indices": PATHS[c],
             "marker_count": marker_counts[i],
             "median_marker_distance_to_black_line_pt": costs[i][i]}
            for i, c in enumerate(CONCENTRATIONS)
        ],
    }
    page.get_pixmap(matrix=fitz.Matrix(10, 10), clip=fitz.Rect(*PANEL_CLIP)).save(
        str(outdir / f"{PREFIX}_source_panel.png"))
    doc.close()
    return curves, metadata


def read_digitized(path):
    groups = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            c = float(row["concentration_g_L"])
            groups.setdefault(c, []).append([
                float(row["shear_rate_s_inv"]), float(row["eta_mPa_s"]),
                float(row["x_pdf_pt"]), float(row["y_pdf_pt"]),
            ])
    if sorted(groups) != CONCENTRATIONS:
        raise ValueError("Expected the 17 CMC 250k concentrations")
    curves = []
    for c in CONCENTRATIONS:
        a = np.asarray(sorted(groups[c]))
        curves.append({"c": c, "gamma": a[:, 0], "eta_mPa_s": a[:, 1],
                       "x_pdf": a[:, 2], "y_pdf": a[:, 3]})
    return curves


def fit_curve(curve, rng):
    gamma, eta = curve["gamma"], curve["eta_mPa_s"]
    scale, log_eta = float(eta[-1]), np.log10(eta)

    def prediction(q):
        delta, tau, m, fraction = 10 ** q[0], 10 ** q[1], q[2], q[3]
        return scale * fraction + delta / (1 + (tau * gamma) ** m)

    def residual(q):
        return np.log10(prediction(q)) - log_eta

    solutions = []
    for _ in range(15):
        guess = [np.log10(eta[0]), rng.uniform(-5, 1),
                 rng.uniform(.1, 1.2), rng.uniform(0, .8)]
        solutions.append(least_squares(
            residual, guess, bounds=([-8, -12, .02, 0], [9, 8, 4, 1]),
            max_nfev=3000, ftol=1e-12, xtol=1e-12, gtol=1e-12))
    best = min(solutions, key=lambda s: float(np.sum(s.fun ** 2)))
    if not best.success:
        raise RuntimeError(f"Optimization did not converge for {curve['c']} g/L")
    q = best.x
    eta_inf = scale * q[3]
    at_bound = bool(eta_inf < 1e-7 * scale)
    if at_bound:
        eta_inf = 0.0
    eta0, tau, m = float(eta_inf + 10 ** q[0]), float(10 ** q[1]), float(q[2])
    predicted = cross_viscosity(gamma, eta0, eta_inf, tau, m)
    relative = predicted / eta - 1
    return {
        "concentration_g_L": curve["c"],
        "eta0_mPa_s": eta0, "eta_inf_mPa_s": float(eta_inf),
        "tau_s": tau, "m": m,
        "eta0_Pa_s": eta0 * 1e-3, "eta_inf_Pa_s": float(eta_inf) * 1e-3,
        "digitized_gamma_min_s_inv": float(gamma[0]),
        "digitized_gamma_max_s_inv": float(gamma[-1]),
        "rms_relative_error_percent": float(100 * np.sqrt(np.mean(relative ** 2))),
        "max_relative_error_percent": float(100 * np.max(np.abs(relative))),
        "graphical_points_count": len(gamma),
        "eta_inf_at_zero_bound": at_bound,
        "source_pdf_paths": ";".join(map(str, PATHS[curve["c"]])),
    }


def evaluate(fit, gamma):
    return cross_viscosity(gamma, fit["eta0_mPa_s"], fit["eta_inf_mPa_s"],
                           fit["tau_s"], fit["m"])


def make_figures(curves, fits, outdir, metadata, source_panel):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 12, "axes.titlesize": 12,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.dpi": 210})
    palette = plt.get_cmap("turbo")(np.linspace(.06, .93, len(curves)))
    colors = {c: palette[i] for i, c in enumerate(CONCENTRATIONS)}
    lookup = {f["concentration_g_L"]: f for f in fits}

    def draw(ax, selection, markers=True):
        handles = []
        for curve in curves:
            c = curve["c"]
            if c not in selection:
                continue
            gamma = curve["gamma"]
            dense = np.geomspace(gamma[0], gamma[-1], 350)
            line, = ax.loglog(dense, evaluate(lookup[c], dense), color=colors[c],
                              lw=1.8, label=f"{c:g}")
            handles.append(line)
            if markers:
                take = np.unique(np.linspace(0, len(gamma) - 1, 22).astype(int))
                ax.plot(gamma[take], curve["eta_mPa_s"][take], "o", ms=2.8,
                        mfc="white", mec=colors[c], mew=.75, zorder=3)
        ax.grid(which="major", color="#e3e7eb", lw=.7)
        ax.set_xlabel(r"Shear rate, $\dot\gamma$ (s$^{-1}$)")
        ax.set_ylabel(r"Viscosity, $\eta$ (mPa s)")
        return handles

    overview, ax = plt.subplots(figsize=(10.7, 7.3))
    overview.subplots_adjust(left=.105, right=.78, top=.87, bottom=.18)
    handles = draw(ax, CONCENTRATIONS)
    ax.set(xlim=(1e-2, 1e3), ylim=(1, 1e4))
    ax.legend(handles=handles, title="CMC 250k\nConcentration (g/L)",
              loc="upper left", bbox_to_anchor=(1.025, 1), frameon=False,
              borderaxespad=0, labelspacing=.47)
    overview.suptitle("CMC 250k | Cross-model reconstruction of Figure S2",
                      x=.105, ha="left", y=.965, fontsize=15, weight="bold")
    overview.text(.105, .916, "Lines: reconstructed Cross model    Circles: digitized black-line coordinates",
                  color="#495366", fontsize=10)
    overview.text(.105, .032, "Fits use only the shear-rate interval of each printed black curve.\n"
                  "Graphical reconstruction error is not an experimental uncertainty estimate.",
                  color="#495366", fontsize=9)
    overview.savefig(outdir / f"{PREFIX}_cross_reconstruction.png")

    zoom, axes = plt.subplots(1, 3, figsize=(14, 5.6))
    zoom.subplots_adjust(left=.065, right=.985, bottom=.27, top=.80, wspace=.30)
    selections = [CONCENTRATIONS[:6], CONCENTRATIONS[6:15], CONCENTRATIONS[15:]]
    for ax, selected, title in zip(axes, selections,
                                  ["0.6-6 g/L", "7-20 g/L", "31-42 g/L"]):
        handles = draw(ax, selected)
        ax.set_title(title, loc="left", weight="bold")
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.23),
                  ncol=3, frameon=False, title="Concentration (g/L)", fontsize=9)
    axes[0].set(xlim=(1, 1e3), ylim=(2, 16))
    axes[1].set(xlim=(1, 1e3), ylim=(12, 400))
    axes[2].set(xlim=(1e-2, 1e3), ylim=(250, 1e4))
    for ax, ticks in zip(axes, [[2, 3, 5, 10, 15], [20, 50, 100, 200], [300, 1000, 3000, 10000]]):
        ax.set_yticks(ticks, [f"{x:g}" for x in ticks])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    zoom.suptitle("Closer view | the same 17 fitted curves", x=.065, ha="left",
                   y=.95, fontsize=15, weight="bold")
    zoom.text(.065, .87, "Colored lines: Cross fits    Open circles: coordinates extracted from the published black lines",
              fontsize=10, color="#495366")
    zoom.savefig(outdir / f"{PREFIX}_cross_zoom.png")

    figures = [overview, zoom]
    if metadata is not None and source_panel is not None and source_panel.exists():
        overlay, ax = plt.subplots(figsize=(8, 9))
        overlay.subplots_adjust(left=.04, right=.96, top=.845, bottom=.07)
        clip = metadata["panel_clip_pdf_pt"]
        ax.imshow(plt.imread(source_panel), extent=(clip[0], clip[2], clip[3], clip[1]))
        xa, xb = metadata["x_pdf_equals_a_log10_gamma_plus_b"]
        ya, yb = metadata["y_pdf_equals_a_log10_eta_mPa_s_plus_b"]
        for curve in curves:
            gamma = np.geomspace(curve["gamma"][0], curve["gamma"][-1], 500)
            ax.plot(xa * np.log10(gamma) + xb,
                    ya * np.log10(evaluate(lookup[curve["c"]], gamma)) + yb,
                    color="#e85b12", lw=.95, dashes=(3, 3))
        ax.set(xlim=(clip[0], clip[2]), ylim=(clip[3], clip[1]))
        ax.axis("off")
        overlay.suptitle("Direct overlay on the original CMC 250k panel",
                        y=.98, fontsize=14, weight="bold")
        overlay.text(.5, .925, "Orange dashed: reconstructed Cross curves\n"
                     "Black lines and blue experimental symbols: original Figure S2",
                     ha="center", va="top", color="#495366", fontsize=10)
        overlay.text(.5, .035, "Source: supplied Supporting Information, Figure S2, page S3",
                     ha="center", fontsize=9, color="#495366")
        overlay.savefig(outdir / f"{PREFIX}_cross_overlay.png")
        figures.append(overlay)

    with PdfPages(outdir / f"{PREFIX}_cross_comparison.pdf") as pdf:
        pdf.infodict().update(Title="CMC 250k: digitized Cross-model reconstruction",
                              Subject="SI Figure S2; numerical reconstruction of printed black lines")
        for figure in figures:
            pdf.savefig(figure)
    for figure in figures:
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).with_name(f"{PREFIX}_digitized_black_lines.csv"))
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).parent / "reproduced_cmc250k")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata, source_panel = None, None
    if args.pdf:
        curves, metadata = extract_pdf(args.pdf, args.outdir)
        source_panel = args.outdir / f"{PREFIX}_source_panel.png"
    else:
        curves = read_digitized(args.data)
        meta_path = args.data.with_name(f"{PREFIX}_digitization_metadata.json")
        if meta_path.exists():
            metadata = json.loads(meta_path.read_text())
        source_panel = args.data.with_name(f"{PREFIX}_source_panel.png")

    rng = np.random.default_rng(41)
    fits = [fit_curve(curve, rng) for curve in curves]
    for fit in fits:
        print(f"c={fit['concentration_g_L']:g} g/L: "
              f"eta0={fit['eta0_mPa_s']:.7g}, eta_inf={fit['eta_inf_mPa_s']:.7g} mPa s; "
              f"tau={fit['tau_s']:.7g} s, m={fit['m']:.7g}; "
              f"RMS={fit['rms_relative_error_percent']:.4f}%", flush=True)
    write_csv(args.outdir / f"{PREFIX}_cross_parameters.csv", fits)
    specification = {
        "equation": "eta = eta_inf + (eta0 - eta_inf)/(1 + (tau*gamma_dot)**m)",
        "target": "Printed black fitted curves, not the raw experimental scatter",
        "fit_objective": "Unweighted least squares in log10(viscosity) at unique vector x coordinates",
        "constraints": "eta0 >= eta_inf >= 0; tau in [1e-12,1e8] s; m in [0.02,4]",
        "viscosity_unit_primary": "mPa s", "viscosity_unit_for_LBM_export": "Pa s",
        "shear_rate_unit": "s^-1", "tau_unit": "s", "m_unit": "dimensionless",
        "notes": [
            "Zero eta_inf is a fitted constraint boundary, not a measured zero solvent viscosity.",
            "Individual parameters, especially eta_inf, can be weakly identified from the displayed interval.",
            "The graphical points are not independent experimental observations; no experimental confidence intervals are claimed.",
            "Use each row's digitized shear-rate limits to distinguish reconstruction from extrapolation.",
            "The concentration ordering reversals at 8/9 and 11/12 g/L are present in the source and preserved.",
            "No continuous concentration interpolation has been assumed or fitted.",
        ],
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "parameters": fits,
    }
    (args.outdir / f"{PREFIX}_cross_parameters.json").write_text(
        json.dumps(specification, indent=2), encoding="utf-8")
    rows = []
    dense_rows = []
    for curve, fit in zip(curves, fits):
        for gamma, eta, x, y in zip(curve["gamma"], curve["eta_mPa_s"], curve["x_pdf"], curve["y_pdf"]):
            rows.append({"concentration_g_L": curve["c"], "shear_rate_s_inv": gamma,
                         "eta_mPa_s": eta, "eta_Pa_s": eta * 1e-3,
                         "x_pdf_pt": x, "y_pdf_pt": y})
        for gamma in np.geomspace(curve["gamma"][0], curve["gamma"][-1], 300):
            eta = float(evaluate(fit, gamma))
            dense_rows.append({"concentration_g_L": curve["c"], "shear_rate_s_inv": gamma,
                               "eta_mPa_s": eta, "eta_Pa_s": eta * 1e-3})
    write_csv(args.outdir / f"{PREFIX}_digitized_black_lines.csv", rows)
    write_csv(args.outdir / f"{PREFIX}_reconstructed_curves.csv", dense_rows)
    if metadata:
        (args.outdir / f"{PREFIX}_digitization_metadata.json").write_text(
            json.dumps(metadata, indent=2), encoding="utf-8")
    make_figures(curves, fits, args.outdir, metadata, source_panel)


if __name__ == "__main__":
    main()
