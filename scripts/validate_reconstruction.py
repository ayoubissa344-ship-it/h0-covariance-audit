#!/usr/bin/env python3
"""Compare the reconstructed table with the published/packaged group-level results."""
import csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "python" / "src"))
from h0audit import core
from h0audit.io import read_dataset
d = read_dataset(ROOT / "data" / "h0_compilation_reconstructed.csv")
Li, Ri = core.ivw(d.yL, d.sL), core.ivw(d.yR, d.sR)
Lg, Rg = core.gls_equicorr(d.yL, d.sL, .40), core.gls_equicorr(d.yR, d.sR, .15)
ref = {  # packaged (results/legacy_v1) / published values
 "IVW ladder mean": (Li.mean, 72.7324645383), "IVW ladder SE": (Li.se, 0.3852470770),
 "IVW comparison mean": (Ri.mean, 69.3722261931), "IVW comparison SE": (Ri.se, 0.3448874494),
 "GLS ladder mean": (Lg.mean, 72.3020299956), "GLS ladder SE": (Lg.se, 0.5712933539),
 "GLS comparison mean": (Rg.mean, 69.5205256510), "GLS comparison SE": (Rg.se, 0.4203840062),
 "negative weights ladder": (Lg.n_negative, 16), "negative weights comparison": (Rg.n_negative, 38)}
out = ROOT / "results" / "reconstructed" / "validation.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["quantity", "reconstructed", "packaged", "abs_diff"])
    for k, (a, b) in ref.items():
        w.writerow([k, f"{a:.6f}", f"{b:.6f}", f"{abs(a-b):.6f}"]); print(f"{k:30s} {a:12.6f} {b:12.6f} {abs(a-b):10.6f}")
