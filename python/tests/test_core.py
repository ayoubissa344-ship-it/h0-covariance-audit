"""Run with:  python -m unittest discover -s python/tests   (pytest also works)."""
import json
import subprocess
import unittest
from pathlib import Path

import numpy as np

from h0audit import core
from h0audit.cli import summary
from h0audit.io import parse_latex_table, parse_sigma, read_dataset

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "data" / "synthetic_fixture.csv"
BIN = ROOT / "build" / "h0audit"


def equicorr_C(s, rho):
    C = np.outer(s, s) * rho
    np.fill_diagonal(C, s**2)
    return C


class TestCore(unittest.TestCase):
    def test_closed_form_matches_dense(self):
        for rho in (0.0, 0.15, 0.4, 0.9, -0.01):
            rng = np.random.default_rng(1)
            s, y = rng.uniform(0.4, 3, 30), rng.normal(70, 2, 30)
            a, b = core.gls_equicorr(y, s, rho), core.gls_dense(y, equicorr_C(s, rho))
            self.assertAlmostEqual(a.mean, b.mean, places=9)
            self.assertAlmostEqual(a.se, b.se, places=9)
            self.assertEqual(a.n_negative, b.n_negative)

    def test_rho0_is_ivw(self):
        rng = np.random.default_rng(2)
        s, y = rng.uniform(0.4, 3, 58), rng.normal(70, 2, 58)
        a, b = core.ivw(y, s), core.gls_equicorr(y, s, 0.0)
        self.assertAlmostEqual(a.mean, b.mean, places=10)
        self.assertAlmostEqual(a.se, b.se, places=10)

    def test_negative_weight_threshold_is_exact(self):
        rng = np.random.default_rng(3)
        s, y = rng.uniform(0.4, 3, 58), rng.normal(70, 2, 58)
        g, th = core.gls_equicorr(y, s, 0.15), core.negative_weight_threshold(s, 0.15)
        self.assertTrue(np.array_equal(g.weights < 0, s > th))

    def test_design_effect_equal_sigma(self):
        n, rho, s0 = 40, 0.2, 1.3
        g = core.gls_equicorr(np.full(n, 70.0), np.full(n, s0), rho)
        self.assertAlmostEqual(g.se**2, s0**2 * (1 + (n - 1) * rho) / n, places=12)
        self.assertAlmostEqual(core.information_ratio(np.full(n, s0), rho), 1 / (1 + (n - 1) * rho), places=12)

    def test_cross_limit_is_the_psd_boundary(self):
        nL, nR, rl, rr = 12, 20, 0.4, 0.15
        lim = core.cross_rho_limit(nL, rl, nR, rr)
        mineig = lambda r: np.linalg.eigvalsh(core.block_equicorr_matrix(nL, rl, nR, rr, r)).min()
        self.assertAlmostEqual(mineig(lim), 0.0, places=11)
        self.assertGreater(mineig(0.99 * lim), 0)
        self.assertLess(mineig(1.01 * lim), 0)

    def test_capped_bound_attained_and_never_exceeded(self):
        rng = np.random.default_rng(4)
        q, r = rng.normal(size=6), 0.3
        vmax = core.capped_variance_max(q, r)
        sg = np.sign(q)
        R = (1 - r) * np.eye(6) + r * np.outer(sg, sg)
        self.assertAlmostEqual(q @ R @ q, vmax, places=12)
        for _ in range(3000):
            A = rng.uniform(-r, r, (6, 6)); A = (A + A.T) / 2; np.fill_diagonal(A, 1)
            if np.linalg.eigvalsh(A).min() > 0:
                self.assertLessEqual(q @ A @ q, vmax + 1e-12)

    def test_kappa_monotone_in_z(self):
        L, R = core.Combined(72.3, 0.5713, np.ones(1), 0), core.Combined(69.52, 0.4204, np.ones(1), 0)
        zs = [core.contrast_kappa(L, R, k).z for k in np.linspace(-1, 0.99, 50)]
        self.assertTrue(all(b > a for a, b in zip(zs, zs[1:])))

    def test_packaged_group_level_numbers(self):
        # Reproduces the data-free kappa envelope quoted in the manuscript from the packaged summaries.
        L, R = core.Combined(72.3020299956, 0.5712933539, np.ones(1), 0), core.Combined(69.5205256510, 0.4203840062, np.ones(1), 0)
        self.assertAlmostEqual(core.contrast_kappa(L, R, 0).z, 3.9215, places=3)
        self.assertAlmostEqual(core.contrast_kappa(L, R, -1).z, 2.8048, places=3)
        self.assertAlmostEqual(core.contrast_kappa(L, R, 1).z, 18.4316, places=3)
        self.assertAlmostEqual(core.cross_rho_limit(30, 0.4, 58, 0.15), 0.26297, places=5)

    def test_parse(self):
        self.assertEqual(parse_sigma("0.5"), 0.5)
        self.assertAlmostEqual(parse_sigma("^{+0.6}_{-0.2}"), ((0.36 + 0.04) / 2) ** 0.5)
        rows = parse_latex_table("1 & A & 2020 & 70.1 & 1.2 & m\nbad line\n2 & B & 2021 & 72.0 & 0.8 & n")
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][3], 72.0)

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            core.ivw([70.0, 71.0], [1.0, 0.0])
        with self.assertRaises(ValueError):
            core.gls_equicorr([70.0, 71.0, 72.0], [1, 1, 1], -0.6)

    @unittest.skipUnless(BIN.exists(), "C++ binary not built (run `make cpp`)")
    def test_cpp_matches_python(self):
        out = subprocess.run([str(BIN), "summary", str(FIX)], capture_output=True, text=True, check=True)
        c, p = json.loads(out.stdout), summary(read_dataset(FIX), 0.40, 0.15)
        for k in ("ivw_z", "gls_z", "gls_L", "gls_R", "gls_se_delta", "cross_rho_limit", "neg_threshold_L", "neg_threshold_R"):
            self.assertAlmostEqual(c[k] / p[k], 1.0, places=9, msg=k)
        self.assertEqual((c["neg_weights_L"], c["neg_weights_R"]), (p["neg_weights_L"], p["neg_weights_R"]))
        for a, b in zip(c["rho_lr_scan"], p["rho_lr_scan"]):
            self.assertAlmostEqual(a["z"] / b["z"], 1.0, places=9)


class TestLKJ(unittest.TestCase):
    def test_valid_and_marginal_variance(self):
        from h0audit.lkj import lkj_corr
        rng, d, eta = np.random.default_rng(5), 8, 1.0
        draws = [lkj_corr(d, eta, rng) for _ in range(4000)]
        self.assertTrue(all(np.linalg.eigvalsh(R).min() > -1e-12 for R in draws))
        v = np.mean([R[1, 5] ** 2 for R in draws])
        self.assertAlmostEqual(v, 1 / (2 * eta - 1 + d), delta=0.01)


    def test_lkj_quadratic_form_moments(self):
        # Proposition 7: E[q'Rq]=sum q^2 ; Var = 2[(sum q^2)^2 - sum q^4]/(d+2eta-1)
        from h0audit.lkj import lkj_corr
        rng = np.random.default_rng(11)
        d, eta = 12, 2.0
        q = rng.normal(size=d)
        v = np.array([q @ lkj_corr(d, eta, rng) @ q for _ in range(40000)])
        q2, q4 = np.sum(q**2), np.sum(q**4)
        self.assertAlmostEqual(v.mean() / q2, 1.0, delta=0.01)
        self.assertAlmostEqual(v.var() / (2 * (q2**2 - q4) / (d + 2 * eta - 1)), 1.0, delta=0.05)


class TestReconstruction(unittest.TestCase):
    def test_reconstructed_table_matches_packaged_results(self):
        d = read_dataset(ROOT / "data" / "h0_compilation_reconstructed.csv")
        self.assertEqual((d.yL.size, d.yR.size), (30, 58))
        Li, Lg = core.ivw(d.yL, d.sL), core.gls_equicorr(d.yL, d.sL, 0.40)
        self.assertAlmostEqual(Li.mean, 72.7324645383, places=8)   # ladder: exact
        self.assertAlmostEqual(Lg.mean, 72.3020299956, places=8)
        self.assertEqual(Lg.n_negative, 16)
        Ri, Rg = core.ivw(d.yR, d.sR), core.gls_equicorr(d.yR, d.sR, 0.15)
        self.assertAlmostEqual(Ri.mean, 69.3722261931, delta=0.006)  # comparison: provisional
        self.assertAlmostEqual(Rg.se, 0.4203840062, delta=0.001)
        self.assertEqual(Rg.n_negative, 38)


if __name__ == "__main__":
    unittest.main()
