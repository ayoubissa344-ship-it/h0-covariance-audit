#!/usr/bin/env python3
"""Download the pinned upstream table and write data/h0_compilation.csv (needs network).

Verifies (a) the Git blob SHA of the downloaded file and (b) the 30/58 row counts before writing.
Usage: python scripts/fetch_compilation.py [--out data/h0_compilation.csv]
"""
import argparse, hashlib, re, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "src"))
from h0audit.io import parse_latex_table  # noqa: E402

COMMIT = "5a0ddb5418a58796b39d6c1783df99ab362a9839"
BLOB_SHA = "9c335e27532127d4f10dc9a557de35d2aadb51cc"
URL = ("https://raw.githubusercontent.com/ipantos/Dissecting-the-Hubble-Tension/"
       f"{COMMIT}/scripts/Mean_Gaussian%28v2%29.py")


def table(src: str, name: str) -> str:
    m = re.search(rf'latex_data_{re.escape(name)}\s*=\s*r?"""([\s\S]*?)"""', src)
    if not m:
        raise SystemExit(f"table latex_data_{name} not found in pinned source")
    return m.group(1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(ROOT / "data" / "h0_compilation.csv")); a = ap.parse_args()
    req = urllib.request.Request(URL, headers={"User-Agent": "h0audit/1.0"})
    raw = urllib.request.urlopen(req, timeout=30).read()
    sha = hashlib.sha1(b"blob %d\0" % len(raw) + raw).hexdigest()
    if sha != BLOB_SHA:
        raise SystemExit(f"blob SHA mismatch: got {sha}, expected {BLOB_SHA}")
    src = raw.decode("utf-8")
    L, R = parse_latex_table(table(src, "distance_ladder")), parse_latex_table(table(src, "all"))
    if (len(L), len(R)) != (30, 58):
        raise SystemExit(f"unexpected row counts: ladder={len(L)}, comparison={len(R)} (expected 30/58)")
    with open(a.out, "w", encoding="utf-8", newline="") as f:
        f.write(f"# Pinned source: commit {COMMIT}, blob {BLOB_SHA}; sigma symmetrised as upstream (RMS of asymmetric errors).\n")
        f.write("id,group,h0,sigma,label\n")
        for g, rows in (("L", L), ("R", R)):
            for i, (rid, author, year, h0, sg, _m) in enumerate(rows, 1):
                label = f"{author} ({year})".replace(",", ";")
                f.write(f"{g}{i:02d},{g},{h0},{sg:.6f},{label}\n")
    print("wrote", a.out, "- rows: ladder=30, comparison=58; blob SHA verified")


if __name__ == "__main__":
    main()
