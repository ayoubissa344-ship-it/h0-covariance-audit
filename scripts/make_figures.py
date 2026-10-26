#!/usr/bin/env python3
"""Regenerate all manuscript figures from packaged CSVs and closed-form results.
Usage: python scripts/make_figures.py [--out paper/figures]"""
import argparse, csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "src"))
from h0audit import core  # noqa: E402

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "figure.dpi": 150})
C1, C2, C3 = "#1f4e79", "#c0504d", "#4f7f4f"


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save(fig, out, name):
    fig.savefig(out / f"{name}.pdf"); fig.savefig(out / f"{name}.png", dpi=200); plt.close(fig)


def fig_scenarios(out):
    rows = read(ROOT / "results/legacy_v1/hubble_tension_phase7_scenarios.csv")
    lab = {"all_88_rows": "All rows", "exclude_preliminary_Pogosian_69.37": "Excl. Pogosian\n69.37", "leave_out_three_DESI_DR2_family_rows": "Excl. 3 DESI\nDR2 rows"}
    names = list(lab); ivw_ = [float(r["nominal_Z"]) for r in rows if r["estimator"].startswith("inverse")]
    gls_ = [float(r["nominal_Z"]) for r in rows if r["estimator"].startswith("GLS")]
    x = np.arange(3); fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.bar(x - .2, ivw_, .38, color=C1, label="IVW (independence)"); ax.bar(x + .2, gls_, .38, color=C2, label=r"GLS ($\rho_L{=}0.40,\ \rho_R{=}0.15$)")
    for xi, v in zip(x - .2, ivw_): ax.text(xi, v + .1, f"{v:.2f}", ha="center", fontsize=8)
    for xi, v in zip(x + .2, gls_): ax.text(xi, v + .1, f"{v:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x, [lab[n] for n in names]); ax.set_ylabel(r"nominal $Z$"); ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=8); ax.set_ylim(0, 8.6)
    save(fig, out, "fig_scenarios")


def fig_kappa(out):
    L = core.Combined(72.3020299956, 0.5712933539, np.ones(1), 0); R = core.Combined(69.5205256510, 0.4203840062, np.ones(1), 0)
    k = np.linspace(-1, 0.9, 400); z = [core.contrast_kappa(L, R, x).z for x in k]
    fig, ax = plt.subplots(figsize=(5.6, 3.2)); ax.plot(k, z, color=C1, lw=2)
    ax.axvline(0, color="k", lw=.8, ls=":"); ax.scatter([0, -1], [core.contrast_kappa(L, R, 0).z, core.contrast_kappa(L, R, -1).z], color=C2, zorder=5)
    ax.annotate(r"$\kappa{=}0$: $Z{=}3.92$", (0, 3.92), (0.12, 3.1), arrowprops=dict(arrowstyle="-", lw=.6), fontsize=8)
    ax.annotate(r"$\kappa{=}-1$: $Z{=}2.80$", (-1, 2.80), (-0.85, 4.9), arrowprops=dict(arrowstyle="-", lw=.6), fontsize=8)
    ax.set_xlabel(r"$\kappa=\mathrm{Corr}(\widehat H_L,\widehat H_R)$"); ax.set_ylabel(r"nominal $Z(\kappa)$"); ax.set_ylim(2, 12)
    save(fig, out, "fig_kappa")


def fig_cap(out):
    rows = [r for r in read(ROOT / "results/legacy_v1/hubble_tension_robust_bounds.csv") if r["analysis"] == "absolute_pairwise_cap"]
    r_ = np.array([float(r["rho_abs_cap"]) for r in rows]); z = np.array([float(r["z_lower_nominal"]) for r in rows])
    fig, ax = plt.subplots(figsize=(5.6, 3.2)); ax.plot(r_, z, "o-", color=C1, ms=4)
    for t in (3, 5): ax.axhline(t, color="grey", lw=.7, ls="--"); ax.text(1.0, t + .08, rf"${t}\sigma$", ha="right", fontsize=8)
    ax.set_xscale("symlog", linthresh=0.005); ax.set_xlim(0, 1.15); ax.set_xticks([0, 0.01, 0.1, 0.5, 1.0]); ax.set_xticklabels(["0", "0.01", "0.1", "0.5", "1"]); ax.minorticks_off(); ax.set_xlabel(r"entrywise correlation cap $r$ ($|R_{ij}|\leq r$)"); ax.set_ylabel(r"worst-case $Z_{\min}(r)$")
    save(fig, out, "fig_cap")


def fig_deff(out):
    rho = np.linspace(0, 0.6, 300); fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for n, c in ((30, C1), (58, C2)): ax.plot(rho, n / (1 + (n - 1) * rho), color=c, label=rf"$n={n}$")
    ax.plot(rho[1:], 1 / rho[1:], color="k", ls=":", lw=1, label=r"ceiling $1/\rho$")
    ax.scatter([0.40], [30 / (1 + 29 * .4)], color=C1, zorder=5); ax.scatter([0.15], [58 / (1 + 57 * .15)], color=C2, zorder=5)
    ax.set_ylim(0, 30); ax.set_xlabel(r"common correlation $\rho$"); ax.set_ylabel(r"$n_{\rm eff}=n/(1+(n-1)\rho)$"); ax.legend(frameon=False)
    save(fig, out, "fig_deff")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(ROOT / "paper/figures")); a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for f in (fig_scenarios, fig_kappa, fig_cap, fig_deff): f(out)
    print("figures written to", out)
