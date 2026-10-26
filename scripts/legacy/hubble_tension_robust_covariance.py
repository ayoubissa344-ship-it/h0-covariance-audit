#!/usr/bin/env python3
"""Sharp covariance-uncertainty bounds for the Hubble-tension compilation.

Python 3.10+, standard library only. Designed for Termux/Android too.

The script fetches a PINNED version of the public Pantos & Perivolaropoulos
analysis source, extracts the 30 distance-ladder and 58 comparison estimates,
and computes sensitivity bounds for a fixed inverse-variance weighted contrast.

Important: these are bounds on a *summary-statistic contrast under the stated
covariance sets*, not a replacement for the primary likelihoods, not calibrated
p-values under arbitrary non-Gaussian models, and not proof that any particular
correlation structure is physically plausible.
"""
from __future__ import annotations

import csv
import math
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

# Pin to an immutable commit so later upstream changes do not silently change
# what this script reproduces.
SOURCE_URL = (
    "https://raw.githubusercontent.com/ipantos/Dissecting-the-Hubble-Tension/"
    "5a0ddb5418a58796b39d6c1783df99ab362a9839/scripts/Mean_Gaussian(v2).py"
)
OUTPUT_NAME = "hubble_tension_robust_bounds.csv"

@dataclass(frozen=True)
class Estimate:
    h0: float
    sigma: float
    author: str
    year: str
    method: str


def fetch_source() -> str:
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "H0RobustCovariance/1.0"})
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.read().decode("utf-8")


def extract_table(source: str, name: str) -> str:
    m = re.search(rf"latex_data_{re.escape(name)}\s*=\s*r?\"\"\"([\s\S]*?)\"\"\"", source)
    if not m:
        raise ValueError(f"Table {name!r} not found in pinned source")
    return m.group(1)


def parse_sigma(error_text: str) -> float:
    # Mirrors the upstream script's approximation for asymmetric uncertainties.
    m = re.search(r"\^\{\+?([\d.]+)\}_\{-?([\d.]+)\}", error_text)
    if m:
        up, down = float(m.group(1)), float(m.group(2))
        return math.sqrt((up * up + down * down) / 2.0)
    m = re.search(r"[\d.]+", error_text)
    if not m:
        raise ValueError(f"No uncertainty found in {error_text!r}")
    return float(m.group(0))


def parse_table(raw: str) -> list[Estimate]:
    result: list[Estimate] = []
    for line in raw.strip().splitlines():
        p = [x.strip() for x in line.split("&")]
        if len(p) < 5:
            continue
        try:
            h0, sigma = float(p[3]), parse_sigma(p[4])
        except (ValueError, IndexError):
            continue
        if math.isfinite(h0) and math.isfinite(sigma) and sigma > 0:
            result.append(Estimate(h0, sigma, p[1], p[2], " & ".join(p[5:])))
    return result


def ivw(data: list[Estimate]) -> dict[str, object]:
    raw = [1.0 / (e.sigma * e.sigma) for e in data]
    total = sum(raw)
    weights = [w / total for w in raw]
    mean = sum(a * e.h0 for a, e in zip(weights, data))
    # q_i = a_i sigma_i; q^T q is the variance for independent errors.
    q = [a * e.sigma for a, e in zip(weights, data)]
    return {
        "n": len(data), "mean": mean, "weights": weights,
        "independent_variance": sum(x*x for x in q),
        "sum_abs_q": sum(abs(x) for x in q),
        "sum_q": q,
        "diagonal_se": math.sqrt(sum(x*x for x in q)),
    }


def unrestricted_correlation_cap(left: dict[str, object], right: dict[str, object], delta: float, rho_max: float) -> dict[str, float]:
    """Sharp Vmax when all individual pairwise |rho_ij| <= rho_max.

    For q = (a_i sigma_i on L, -b_j sigma_j on R),
    max_R q'R q = (1-rho)*sum(q_i^2) + rho*(sum(|q_i|))^2.
    The bound is attainable by R=(1-rho)I + rho*s*s', with s_i=sign(q_i).
    This allows negative cross-ensemble correlations; it is a deliberately broad
    worst case, not an astrophysical model.
    """
    if not 0.0 <= rho_max <= 1.0:
        raise ValueError("rho_max must be in [0,1]")
    q2 = float(left["independent_variance"]) + float(right["independent_variance"])
    q1 = float(left["sum_abs_q"]) + float(right["sum_abs_q"])
    vmax = (1.0 - rho_max) * q2 + rho_max * q1 * q1
    return {"rho_abs_cap": rho_max, "delta_H0": delta, "Vmax": vmax,
            "SE_max": math.sqrt(vmax), "Z_min_nominal": abs(delta) / math.sqrt(vmax)}


def nonnegative_within_group_cap(left: dict[str, object], right: dict[str, object], delta: float, rho_l: float, rho_r: float) -> dict[str, float]:
    """Stress test: 0<=rho_ij<=rho_group within each group, cross-group cov=0.

    For fixed positive IVW weights, maximum variance in a group occurs at the
    equicorrelation boundary rho_ij=rho_group. Cross-group covariance is set to
    zero, which is worst-case for the contrast if only nonnegative correlations
    are permitted and the weights in both groups are nonnegative.
    """
    if not 0 <= rho_l <= 1 or not 0 <= rho_r <= 1:
        raise ValueError("group rho caps must be in [0,1]")
    vl0, vr0 = float(left["independent_variance"]), float(right["independent_variance"])
    sl, sr = float(left["sum_abs_q"]), float(right["sum_abs_q"])
    vl = (1-rho_l)*vl0 + rho_l*sl*sl
    vr = (1-rho_r)*vr0 + rho_r*sr*sr
    variance = vl + vr
    return {"rho_l": rho_l, "rho_r": rho_r, "delta_H0": delta,
            "V_ladder": vl, "V_comparison": vr,
            "SE_ladder": math.sqrt(vl), "SE_comparison": math.sqrt(vr),
            "SE_delta": math.sqrt(variance), "Z_nominal": abs(delta)/math.sqrt(variance)}


def self_test() -> None:
    # Rank-one correlation construction reaches the claimed upper bound.
    q = [1.0, -2.0, 0.5]
    for rho in (0.0, 0.1, 0.5, 1.0):
        s = [1.0 if x >= 0 else -1.0 for x in q]
        R = [[(1.0 if i == j else 0.0)*(1-rho) + rho*s[i]*s[j] for j in range(3)] for i in range(3)]
        direct = sum(q[i]*R[i][j]*q[j] for i in range(3) for j in range(3))
        formula = (1-rho)*sum(x*x for x in q) + rho*(sum(abs(x) for x in q)**2)
        assert abs(direct-formula) < 1e-10, (rho, direct, formula)
    print("SELF-TEST PASS: sharp covariance upper-bound identity verified on synthetic cases.")


def main() -> int:
    if "--self-test" in sys.argv:
        self_test()
        return 0
    try:
        source = fetch_source()
        ladder = parse_table(extract_table(source, "distance_ladder"))
        other = parse_table(extract_table(source, "all"))
    except Exception as exc:
        print(f"ERROR loading the pinned source: {exc}", file=sys.stderr)
        print("Check internet access. Self-test works offline: python hubble_tension_robust_covariance.py --self-test", file=sys.stderr)
        return 2
    if (len(ladder), len(other)) != (30, 58):
        raise RuntimeError(f"Unexpected source table sizes: ladder={len(ladder)}, comparison={len(other)}")

    L, R = ivw(ladder), ivw(other)
    delta = float(L["mean"]) - float(R["mean"])
    z_ind = abs(delta) / math.sqrt(float(L["independent_variance"]) + float(R["independent_variance"]))

    rows: list[dict[str, object]] = []
    for rho_max in (0, 0.005, 0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0):
        row = unrestricted_correlation_cap(L, R, delta, rho_max)
        row.update({"analysis": "absolute_pairwise_cap", "n_ladder": len(ladder), "n_comparison": len(other),
                    "assumption": "all pairwise |correlation|<=rho_max; cross-group anticorrelation allowed"})
        rows.append(row)

    for rho_l in (0, 0.10, 0.20, 0.30, 0.40, 0.50):
        for rho_r in (0, 0.05, 0.10, 0.15, 0.20, 0.25):
            row = nonnegative_within_group_cap(L, R, delta, rho_l, rho_r)
            row.update({"analysis": "nonnegative_within_group_stress", "n_ladder": len(ladder), "n_comparison": len(other),
                        "assumption": "within-group equicorrelation at stated rho; cross-group covariance=0"})
            rows.append(row)

    output = Path.cwd() / OUTPUT_NAME
    keys = sorted({key for row in rows for key in row})
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    print("Pinned source:", SOURCE_URL)
    print(f"Parsed counts: ladder={len(ladder)}, comparison={len(other)}")
    print(f"Fixed IVW mean L={float(L['mean']):.6f} +/- {math.sqrt(float(L['independent_variance'])):.6f}")
    print(f"Fixed IVW mean R={float(R['mean']):.6f} +/- {math.sqrt(float(R['independent_variance'])):.6f}")
    print(f"Delta H0={delta:.6f}; independent nominal Z={z_ind:.4f}")
    print("\nSharp worst-case nominal Z lower bounds when all pairwise |rho_ij| <= r:")
    for rho_max in (0, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0):
        x = unrestricted_correlation_cap(L, R, delta, rho_max)
        print(f"  r={rho_max:>5.3f}: SE_delta,max={x['SE_max']:.4f}, Z_min={x['Z_min_nominal']:.3f}")
    sample = nonnegative_within_group_cap(L, R, delta, 0.40, 0.15)
    print("\nWithin-group positive-equicorrelation stress case (rho_L=0.40, rho_R=0.15; cross covariance=0):")
    print(f"  SE_L={sample['SE_ladder']:.4f}, SE_R={sample['SE_comparison']:.4f}, Z_nominal={sample['Z_nominal']:.3f}")
    print("\nWrote", output)
    print("WARNING: these are covariance-set sensitivity bounds for fixed IVW weights, not calibrated p-values or a final physical inference.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
