#!/usr/bin/env python3
"""Phase 7: provenance and covariance sensitivity audit for the 88-row H0 compilation.

Standard-library-only; designed for Python 3.10+ and Termux/Android.
The upstream source is pinned to an immutable Git commit. The script does not
re-fit primary cosmological likelihoods and does not interpret nominal Z values
as calibrated p-values when covariance assumptions are uncertain.

Usage:
  python hubble_tension_phase7_provenance_audit.py --self-test
  python hubble_tension_phase7_provenance_audit.py

A full run needs network access to fetch the pinned upstream table.
"""
from __future__ import annotations

import csv
import math
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

SOURCE_COMMIT = "5a0ddb5418a58796b39d6c1783df99ab362a9839"
SOURCE_URL = (
    "https://raw.githubusercontent.com/ipantos/Dissecting-the-Hubble-Tension/"
    + SOURCE_COMMIT
    + "/scripts/Mean_Gaussian%28v2%29.py"
)
OUT_CSV = "hubble_tension_phase7_scenarios.csv"
OUT_LEVERAGE = "hubble_tension_phase7_source_leverage.csv"


@dataclass(frozen=True)
class Estimate:
    row_id: str
    author: str
    year: str
    h0: float
    sigma: float
    method: str = ""


def fetch_source() -> str:
    request = urllib.request.Request(
        SOURCE_URL, headers={"User-Agent": "H0ProvenanceAudit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def extract_table(source: str, name: str) -> str:
    pattern = rf"latex_data_{re.escape(name)}\s*=\s*r?\"\"\"([\s\S]*?)\"\"\""
    match = re.search(pattern, source)
    if not match:
        raise ValueError(f"Could not find latex_data_{name!r} in pinned source")
    return match.group(1)


def parse_sigma(text: str) -> float:
    """Use the upstream script's RMS symmetrization for asymmetric 1-sigma errors."""
    match = re.search(r"\^\{\+?([\d.]+)\}_\{-?([\d.]+)\}", text)
    if match:
        up, down = float(match.group(1)), float(match.group(2))
        return math.sqrt((up * up + down * down) / 2.0)
    match = re.search(r"[\d.]+", text)
    if not match:
        raise ValueError(f"No uncertainty found in {text!r}")
    return float(match.group(0))


def parse_table(raw: str) -> list[Estimate]:
    output: list[Estimate] = []
    for line in raw.strip().splitlines():
        if "&" not in line:
            continue
        parts = [part.strip() for part in line.split("&")]
        if len(parts) < 5:
            continue
        try:
            estimate = Estimate(
                row_id=parts[0],
                author=parts[1],
                year=parts[2],
                h0=float(parts[3]),
                sigma=parse_sigma(parts[4]),
                method=" & ".join(parts[5:]).strip(),
            )
        except (ValueError, IndexError):
            continue
        if math.isfinite(estimate.h0) and math.isfinite(estimate.sigma) and estimate.sigma > 0:
            output.append(estimate)
    return output


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Gauss-Jordan solve with partial pivoting; adequate for these 30/58 matrices."""
    n = len(matrix)
    if n == 0 or len(rhs) != n or any(len(row) != n for row in matrix):
        raise ValueError("Matrix dimensions do not match")
    aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1e-14:
            raise ValueError("Covariance matrix is singular or ill-conditioned")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        for j in range(col, n + 1):
            aug[col][j] /= scale
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor:
                for j in range(col, n + 1):
                    aug[row][j] -= factor * aug[col][j]
    return [row[n] for row in aug]


def ivw(data: list[Estimate]) -> dict[str, float]:
    raw_weights = [1.0 / (item.sigma * item.sigma) for item in data]
    total = sum(raw_weights)
    weights = [weight / total for weight in raw_weights]
    mean = sum(weight * item.h0 for weight, item in zip(weights, data))
    return {"n": len(data), "mean": mean, "se": math.sqrt(1.0 / total), "negative_weights": 0.0}


def covariance(data: list[Estimate], rho: float) -> list[list[float]]:
    n = len(data)
    if not -1.0 / max(n - 1, 1) < rho < 1.0:
        raise ValueError("rho does not define a positive-definite equicorrelation matrix")
    return [
        [data[i].sigma * data[j].sigma * (1.0 if i == j else rho) for j in range(n)]
        for i in range(n)
    ]


def gls(data: list[Estimate], rho: float) -> tuple[dict[str, float], list[float]]:
    c = covariance(data, rho)
    c_inv_one = solve_linear(c, [1.0] * len(data))
    normalization = sum(c_inv_one)
    if normalization <= 0:
        raise ArithmeticError("Non-positive GLS normalization")
    weights = [x / normalization for x in c_inv_one]
    mean = sum(weight * item.h0 for weight, item in zip(weights, data))
    return {
        "n": len(data),
        "mean": mean,
        "se": math.sqrt(1.0 / normalization),
        "negative_weights": float(sum(weight < -1e-10 for weight in weights)),
    }, weights


def is_preliminary_pogosian(item: Estimate) -> bool:
    # The exact 69.37 +/- 0.65 row is labelled as 2025/private communication
    # in the paper, but the repository's compact table labels the row 2024 and
    # points to arXiv:2405.20306 (the earlier DESI DR1 analysis).
    return item.author.lower().startswith("pogosian") and abs(item.h0 - 69.37) < 0.01 and abs(item.sigma - 0.65) < 0.02


def is_garcia_desi_dr2(item: Estimate) -> bool:
    return item.author.lower().startswith("garcia") and abs(item.h0 - 70.03) < 0.03


def is_du_mapage_desi_dr2(item: Estimate) -> bool:
    return item.author.lower().startswith("du") and abs(item.h0 - 71.59) < 0.03


def nominal_z(left: dict[str, float], right: dict[str, float]) -> tuple[float, float, float]:
    delta = left["mean"] - right["mean"]
    se_delta = math.sqrt(left["se"] ** 2 + right["se"] ** 2)  # zero cross-group covariance
    return delta, se_delta, abs(delta) / se_delta


def make_row(label: str, estimator: str, left: dict[str, float], right: dict[str, float], rho_l: float, rho_r: float) -> dict[str, object]:
    delta, se_delta, z = nominal_z(left, right)
    return {
        "scenario": label,
        "estimator": estimator,
        "n_distance_ladder": int(left["n"]),
        "n_comparison": int(right["n"]),
        "rho_ladder": rho_l,
        "rho_comparison": rho_r,
        "H0_ladder": left["mean"],
        "SE_ladder": left["se"],
        "H0_comparison": right["mean"],
        "SE_comparison": right["se"],
        "Delta_ladder_minus_comparison": delta,
        "SE_Delta_assuming_zero_cross_group_cov": se_delta,
        "nominal_Z": z,
        "negative_GLS_weights_ladder": int(left.get("negative_weights", 0)),
        "negative_GLS_weights_comparison": int(right.get("negative_weights", 0)),
        "interpretation": "Summary-level sensitivity; not a primary-likelihood refit or calibrated p-value under unknown covariance.",
    }


def run(source: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    ladder = parse_table(extract_table(source, "distance_ladder"))
    other = parse_table(extract_table(source, "all"))
    if (len(ladder), len(other)) != (30, 58):
        raise RuntimeError(f"Pinned source table count changed: ladder={len(ladder)}, other={len(other)}")

    preliminary = [item for item in other if is_preliminary_pogosian(item)]
    if len(preliminary) != 1:
        raise RuntimeError(f"Expected one preliminary Pogosian row; found {len(preliminary)}")
    desi_family_predicates = (is_preliminary_pogosian, is_garcia_desi_dr2, is_du_mapage_desi_dr2)
    family = [item for item in other if any(pred(item) for pred in desi_family_predicates)]
    if len(family) != 3:
        raise RuntimeError(f"Expected 3 DESI DR2 family rows for this audit; found {len(family)}")

    rows: list[dict[str, object]] = []
    # Independent inverse-variance means (the published summary-table baseline).
    l_ivw, r_ivw = ivw(ladder), ivw(other)
    rows.append(make_row("all_88_rows", "inverse-variance weighted, independence", l_ivw, r_ivw, 0.0, 0.0))
    no_pre = [x for x in other if not is_preliminary_pogosian(x)]
    rows.append(make_row("exclude_preliminary_Pogosian_69.37", "inverse-variance weighted, independence", l_ivw, ivw(no_pre), 0.0, 0.0))
    no_family = [x for x in other if not any(pred(x) for pred in desi_family_predicates)]
    rows.append(make_row("leave_out_three_DESI_DR2_family_rows", "inverse-variance weighted, independence", l_ivw, ivw(no_family), 0.0, 0.0))

    # Replicate the paper's illustrative two-block equicorrelation model.
    l_gls, l_weights = gls(ladder, 0.40)
    r_gls, r_weights = gls(other, 0.15)
    rows.append(make_row("all_88_rows", "GLS with assumed equicorrelation", l_gls, r_gls, 0.40, 0.15))
    r_no_pre, _ = gls(no_pre, 0.15)
    rows.append(make_row("exclude_preliminary_Pogosian_69.37", "GLS with assumed equicorrelation", l_gls, r_no_pre, 0.40, 0.15))
    r_no_family, _ = gls(no_family, 0.15)
    rows.append(make_row("leave_out_three_DESI_DR2_family_rows", "GLS with assumed equicorrelation", l_gls, r_no_family, 0.40, 0.15))

    raw = [1.0 / (item.sigma * item.sigma) for item in other]
    total = sum(raw)
    by_key = []
    for item in other:
        if is_preliminary_pogosian(item):
            label, pred, status = "Pogosian 69.37 +/- 0.65", is_preliminary_pogosian, "private communication / seminar value; paper says publication forthcoming"
        elif is_garcia_desi_dr2(item):
            label, pred, status = "Garcia Escudero et al. 70.03 +/- 0.97", is_garcia_desi_dr2, "public paper; uses DESI DR2 plus independent additional likelihoods"
        elif is_du_mapage_desi_dr2(item):
            label, pred, status = "Du et al. MAPAge 71.59 +/- 0.94", is_du_mapage_desi_dr2, "public paper; uses DESI DR2 plus DESY5/SGL/CC/GRB inputs"
        else:
            continue
        by_key.append({
            "estimate": label,
            "row_id_in_58_table": item.row_id,
            "repository_year_field": item.year,
            "H0": item.h0,
            "sigma": item.sigma,
            "method_and_citation_in_code": item.method,
            "inverse_variance_weight_pct_of_58": 100.0 * (1.0 / item.sigma**2) / total,
            "status_and_overlap_note": status,
        })
    if len(by_key) != 3:
        raise RuntimeError("Could not construct the three-row provenance/leverage table")
    return rows, by_key


def write_csv(path: str, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def self_test() -> None:
    sample = [Estimate("1", "A", "2020", 70.0, 1.0), Estimate("2", "B", "2021", 72.0, 2.0)]
    a = ivw(sample)
    b, weights = gls(sample, 0.0)
    assert abs(a["mean"] - b["mean"]) < 1e-10
    assert abs(a["se"] - b["se"]) < 1e-10
    assert abs(sum(weights) - 1.0) < 1e-10
    # An equicorrelation matrix with positive rho is SPD and the GLS variance is positive.
    c, w = gls(sample, 0.2)
    assert c["se"] > 0 and abs(sum(w) - 1.0) < 1e-10
    assert abs(nominal_z({"mean": 72.0, "se": 1.0}, {"mean": 70.0, "se": 1.0})[2] - math.sqrt(2.0)) < 1e-10
    print("SELF-TEST PASS: IVW/GLS consistency at rho=0, normalized weights, positive covariance case, and contrast arithmetic.")


def main() -> int:
    if "--self-test" in sys.argv:
        self_test()
        return 0
    try:
        source = fetch_source()
        scenarios, leverage = run(source)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        print("The self-test is available offline: python hubble_tension_phase7_provenance_audit.py --self-test", file=sys.stderr)
        return 2
    write_csv(OUT_CSV, scenarios)
    write_csv(OUT_LEVERAGE, leverage)
    print("Pinned upstream source:", SOURCE_URL)
    print("Rows checked: distance ladder=30, comparison=58")
    print("Wrote", OUT_CSV, "and", OUT_LEVERAGE)
    for row in scenarios:
        print(f"{row['scenario']} | {row['estimator']} | Delta={row['Delta_ladder_minus_comparison']:.4f} +/- {row['SE_Delta_assuming_zero_cross_group_cov']:.4f} | Z={row['nominal_Z']:.3f}")
    print("Note: Z values are nominal summary-level Gaussian equivalents under the specified covariance assumptions.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
