#!/usr/bin/env python3
"""Check graphical sensitivity to eta_inf, using the bundled digitized data.

python3 check_eta_inf_sensitivity.py

Output ranges are COARSE GRAPHICAL tolerance ranges, NOT statistical
confidence intervals or uncertainty estimates for experimental observations.
The fixed value 1 mPa s is illustrative, not a measured solvent viscosity.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
from fit_cmc250k import read_digitized, write_csv


def refit(gamma, eta, eta_inf, fit, warm=None):
    log_eta, log_gamma = np.log(eta), np.log(gamma)

    def calc(q):
        return eta_inf + np.exp(q[0]) * expit(-np.exp(q[2]) * (q[1] + log_gamma))

    def residual(q):
        return np.log(calc(q)) - log_eta

    low, high = [-20, -27.6, np.log(.01)], [30, 27.6, np.log(3)]
    starts = [np.log([max(fit['eta0_mPa_s'] - eta_inf, .01), fit['tau_s'], fit['m']]),
              np.log([max(eta.max() * 1.1 - eta_inf, .01), 1 / gamma.max(), .7])]
    if warm is not None:
        starts.append(warm)
    results = [least_squares(residual, np.clip(q, low, high), bounds=(low, high),
                             max_nfev=2000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
               for q in starts]
    result = min(results, key=lambda r: np.sum(r.fun ** 2))
    error = calc(result.x) / eta - 1
    return {
        'eta0_mPa_s': float(eta_inf + np.exp(result.x[0])),
        'eta_inf_mPa_s': float(eta_inf), 'tau_s': float(np.exp(result.x[1])),
        'm': float(np.exp(result.x[2])),
        'rms_relative_error_percent': float(100 * np.sqrt(np.mean(error ** 2))),
        'max_relative_error_percent': float(100 * np.max(np.abs(error))),
    }, result.x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).parent)
    parser.add_argument('--outdir', type=Path, default=Path(__file__).parent / 'sensitivity_reproduced')
    args = parser.parse_args()
    curves = read_digitized(args.data_dir / 'cmc250k_digitized_black_lines.csv')
    params = json.loads((args.data_dir / 'cmc250k_cross_parameters.json').read_text())['parameters']
    fits = {f['concentration_g_L']: f for f in params}
    full, summary = [], []
    for curve in curves:
        c, gamma, eta = curve['c'], curve['gamma'], curve['eta_mPa_s']
        fit = fits[c]
        fixed, _ = refit(gamma, eta, 1., fit)
        grid, warm = [], None
        for eta_inf in np.unique(np.r_[0, 1, fit['eta_inf_mPa_s'], np.linspace(0, .99 * eta.min(), 41)]):
            row, warm = refit(gamma, eta, eta_inf, fit, warm)
            grid.append(row)
        acceptable = [v['eta_inf_mPa_s'] for v in grid if v['max_relative_error_percent'] <= 1.]
        summary.append({
            'concentration_g_L': c,
            'free_rms_error_percent': fit['rms_relative_error_percent'],
            'free_max_error_percent': fit['max_relative_error_percent'],
            'fixed1_rms_error_percent': fixed['rms_relative_error_percent'],
            'fixed1_max_error_percent': fixed['max_relative_error_percent'],
            'fixed1_eta0_mPa_s': fixed['eta0_mPa_s'],
            'fixed1_tau_s': fixed['tau_s'], 'fixed1_m': fixed['m'],
            'coarse_eta_inf_min_mPa_s_at_1percent_max_error': min(acceptable) if acceptable else None,
            'coarse_eta_inf_max_mPa_s_at_1percent_max_error': max(acceptable) if acceptable else None,
        })
        full.append({'concentration_g_L': c, 'fixed_eta_inf_1_mPa_s': fixed, 'coarse_grid': grid})
        print(f'{c:g} g/L: fixed eta_inf=1 gives max graphical error '
              f'{fixed["max_relative_error_percent"]:.3f}%', flush=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / 'cmc250k_eta_inf_sensitivity.csv', summary)
    (args.outdir / 'cmc250k_eta_inf_sensitivity.json').write_text(json.dumps({
        'notes': ['Fixed eta_inf=1 mPa s is illustrative, not measured.',
                  'Ranges are sampled-grid results after log-least-squares refits, not confidence intervals.',
                  'Acceptance criterion: <=1% maximum absolute relative error to digitized black-line coordinates.'],
        'results': full,
    }, indent=2))


if __name__ == '__main__':
    main()
