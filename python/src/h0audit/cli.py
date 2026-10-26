"""`python -m h0audit.cli summary DATA.csv` -- pure-Python mirror of the C++ CLI (summary, cap)."""
from __future__ import annotations

import argparse
import json

from . import core
from .io import read_dataset


def summary(d, rl, rr):
    Li, Ri = core.ivw(d.yL, d.sL), core.ivw(d.yR, d.sR)
    Lg, Rg = core.gls_equicorr(d.yL, d.sL, rl), core.gls_equicorr(d.yR, d.sR, rr)
    ci, cg = core.contrast(Li, Ri), core.contrast(Lg, Rg)
    lim = core.cross_rho_limit(d.yL.size, rl, d.yR.size, rr)
    scan = []
    for k in range(-10, 11):
        rlr = lim * k / 10
        cov = core.cross_cov(Lg.weights, d.sL, Rg.weights, d.sR, rlr)
        scan.append({"rho_lr": rlr, "kappa": cov / (Lg.se * Rg.se), "z": core.contrast(Lg, Rg, cov).z})
    return {
        "nL": int(d.yL.size), "nR": int(d.yR.size),
        "ivw_L": Li.mean, "ivw_L_se": Li.se, "ivw_R": Ri.mean, "ivw_R_se": Ri.se,
        "ivw_delta": ci.delta, "ivw_se_delta": ci.se, "ivw_z": ci.z,
        "gls_L": Lg.mean, "gls_L_se": Lg.se, "gls_R": Rg.mean, "gls_R_se": Rg.se,
        "gls_delta": cg.delta, "gls_se_delta": cg.se, "gls_z": cg.z,
        "neg_weights_L": Lg.n_negative, "neg_weights_R": Rg.n_negative,
        "neg_threshold_L": core.negative_weight_threshold(d.sL, rl),
        "neg_threshold_R": core.negative_weight_threshold(d.sR, rr),
        "cross_rho_limit": lim, "rho_lr_scan": scan,
        "kappa_envelope": {"z_at_kappa_minus1": core.contrast_kappa(Lg, Rg, -1).z,
                           "z_at_kappa_plus1": core.contrast_kappa(Lg, Rg, 1).z},
    }


def main(argv=None):
    p = argparse.ArgumentParser(prog="h0audit")
    p.add_argument("cmd", choices=["summary"])
    p.add_argument("data")
    p.add_argument("--rho-l", type=float, default=0.40)
    p.add_argument("--rho-r", type=float, default=0.15)
    a = p.parse_args(argv)
    print(json.dumps(summary(read_dataset(a.data), a.rho_l, a.rho_r), indent=2))


if __name__ == "__main__":
    main()
