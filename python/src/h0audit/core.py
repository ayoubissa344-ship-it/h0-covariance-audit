"""Reference (NumPy) implementation. The C++ library mirrors these functions; the test-suite
cross-validates the two. Notation follows docs/theory.md."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Combined:
    mean: float
    se: float
    weights: np.ndarray
    n_negative: int


@dataclass(frozen=True)
class Contrast:
    delta: float
    var: float
    se: float
    z: float


def _check(y, s):
    y, s = np.asarray(y, float), np.asarray(s, float)
    if y.shape != s.shape or y.ndim != 1 or y.size == 0:
        raise ValueError("y and s must be equal-length 1-D arrays")
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(s)) and np.all(s > 0)):
        raise ValueError("non-finite estimate or non-positive sigma")
    return y, s


def ivw(y, s) -> Combined:
    y, s = _check(y, s)
    a = 1.0 / s**2
    w = a / a.sum()
    return Combined(float(w @ y), float(a.sum() ** -0.5), w, 0)


def gls_equicorr(y, s, rho: float) -> Combined:
    """O(n) closed form via Sherman-Morrison (Theorem 1 of docs/theory.md)."""
    y, s = _check(y, s)
    n = y.size
    if not (-1.0 / (n - 1) < rho < 1.0 if n > 1 else True):
        raise ValueError("rho outside (-1/(n-1), 1)")
    a = 1.0 / s
    c = rho / (1.0 + (n - 1) * rho)
    u = (a**2 - c * a * a.sum()) / (1.0 - rho)
    N = u.sum()
    if N <= 0:
        raise ArithmeticError("non-positive GLS normalisation")
    w = u / N
    return Combined(float(w @ y), float(N**-0.5), w, int(np.sum(w < -1e-12)))


def gls_dense(y, C) -> Combined:
    y = np.asarray(y, float)
    C = np.asarray(C, float)
    L = np.linalg.cholesky(C)  # raises LinAlgError if not PD
    u = np.linalg.solve(C, np.ones_like(y)) if L is None else np.linalg.solve(C, np.ones_like(y))
    N = u.sum()
    if N <= 0:
        raise ArithmeticError("non-positive GLS normalisation")
    w = u / N
    return Combined(float(w @ y), float(N**-0.5), w, int(np.sum(w < -1e-12)))


def negative_weight_threshold(s, rho: float) -> float:
    """Rows with sigma above this value get a negative equicorrelation-GLS weight (rho > 0)."""
    s = np.asarray(s, float)
    if rho <= 0:
        return float("inf")
    return float((1.0 + (s.size - 1) * rho) / (rho * np.sum(1.0 / s)))


def information_ratio(s, rho: float) -> float:
    """Var_independent / Var_GLS. Equals 1/(1+(n-1)rho) for equal sigmas; with strongly
    heteroscedastic rows it can exceed 1 (correlated errors can be partly cancelled)."""
    s = np.asarray(s, float)
    n = s.size
    a = 1.0 / s
    c = rho / (1.0 + (n - 1) * rho)
    N = (np.sum(a**2) - c * a.sum() ** 2) / (1.0 - rho)
    return float(N / np.sum(a**2))


def cross_rho_limit(nL: int, rhoL: float, nR: int, rhoR: float) -> float:
    """Exact bound |rho_LR| <= sqrt((1+(nL-1)rhoL)(1+(nR-1)rhoR)/(nL nR))."""
    a, b = 1 + (nL - 1) * rhoL, 1 + (nR - 1) * rhoR
    if not (a > 0 and b > 0 and rhoL < 1 and rhoR < 1):
        raise ValueError("blocks are not positive semidefinite")
    return float(np.sqrt(a * b / (nL * nR)))


def cross_cov(wL, sL, wR, sR, rho_lr: float) -> float:
    return float(rho_lr * (np.asarray(wL) @ np.asarray(sL)) * (np.asarray(wR) @ np.asarray(sR)))


def contrast(L: Combined, R: Combined, cov_lr: float = 0.0) -> Contrast:
    var = L.se**2 + R.se**2 - 2.0 * cov_lr
    if var <= 0:
        raise ValueError("non-positive contrast variance (covariance infeasible)")
    d = L.mean - R.mean
    return Contrast(d, var, var**0.5, abs(d) / var**0.5)


def contrast_kappa(L: Combined, R: Combined, kappa: float) -> Contrast:
    return contrast(L, R, kappa * L.se * R.se)


def capped_variance_max(q, r: float) -> float:
    """Sharp max of q^T R q over correlation matrices with |R_ij| <= r (i != j)."""
    q = np.asarray(q, float)
    if not 0.0 <= r <= 1.0:
        raise ValueError("r must lie in [0,1]")
    return float((1 - r) * np.sum(q**2) + r * np.sum(np.abs(q)) ** 2)


def block_equicorr_matrix(nL, rhoL, nR, rhoR, rho_lr) -> np.ndarray:
    d = nL + nR
    R = np.empty((d, d))
    R[:nL, :nL] = rhoL
    R[nL:, nL:] = rhoR
    R[:nL, nL:] = rho_lr
    R[nL:, :nL] = rho_lr
    np.fill_diagonal(R, 1.0)
    return R
