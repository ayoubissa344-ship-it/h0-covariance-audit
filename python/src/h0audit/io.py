"""Data input: CSV reader and parser for the pinned upstream LaTeX-style table."""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Dataset:
    ids: list[str]
    yL: np.ndarray
    sL: np.ndarray
    yR: np.ndarray
    sR: np.ndarray


def read_dataset(path: str | Path) -> Dataset:
    ids, g = [], {"L": ([], []), "R": ([], [])}
    with open(path, newline="", encoding="utf-8") as fh:
        rows = (ln for ln in fh if ln.strip() and not ln.startswith("#"))
        for rec in csv.DictReader(rows):
            grp = rec["group"].strip()
            if grp not in g:
                raise ValueError(f"group must be L or R, got {grp!r}")
            y, s = float(rec["h0"]), float(rec["sigma"])
            if not (math.isfinite(y) and math.isfinite(s) and s > 0):
                raise ValueError(f"invalid row {rec}")
            g[grp][0].append(y)
            g[grp][1].append(s)
            ids.append(rec.get("id", str(len(ids) + 1)))
    if min(len(g["L"][0]), len(g["R"][0])) < 2:
        raise ValueError("need >= 2 rows per group")
    return Dataset(ids, *(np.array(v) for v in (g["L"][0], g["L"][1], g["R"][0], g["R"][1])))


_ASYM = re.compile(r"\^\{\+?([\d.]+)\}_\{-?([\d.]+)\}")


def parse_sigma(text: str) -> float:
    """RMS symmetrisation of asymmetric errors, as in the pinned upstream script."""
    m = _ASYM.search(text)
    if m:
        up, dn = float(m.group(1)), float(m.group(2))
        return math.sqrt((up * up + dn * dn) / 2.0)
    m = re.search(r"[\d.]+", text)
    if not m:
        raise ValueError(f"no uncertainty in {text!r}")
    return float(m.group(0))


def parse_latex_table(raw: str):
    """Return list of (row_id, author, year, h0, sigma, method) from an '&'-separated table."""
    out = []
    for line in raw.strip().splitlines():
        parts = [p.strip() for p in line.split("&")]
        if len(parts) < 5:
            continue
        try:
            h0, sg = float(parts[3]), parse_sigma(parts[4])
        except (ValueError, IndexError):
            continue
        if math.isfinite(h0) and math.isfinite(sg) and sg > 0:
            out.append((parts[0], parts[1], parts[2], h0, sg, " & ".join(parts[5:]).strip()))
    return out
