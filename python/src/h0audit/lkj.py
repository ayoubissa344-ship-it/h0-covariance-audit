"""Pure-Python/NumPy LKJ sampler (C-vine) -- reference and speed baseline for the C++ version."""
from __future__ import annotations

import numpy as np


def lkj_corr(d: int, eta: float, rng: np.random.Generator) -> np.ndarray:
    P = np.zeros((d, d))
    S = np.eye(d)
    beta = eta + 0.5 * (d - 1)
    for k in range(d - 1):
        beta -= 0.5
        for i in range(k + 1, d):
            p = 2.0 * rng.beta(beta, beta) - 1.0
            P[k, i] = p
            for l in range(k - 1, -1, -1):
                p = p * np.sqrt((1 - P[l, i] ** 2) * (1 - P[l, k] ** 2)) + P[l, i] * P[l, k]
            S[k, i] = S[i, k] = p
    return S


def mc_variances(q, eta: float, draws: int, seed: int = 1) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q = np.asarray(q, float)
    return np.array([q @ lkj_corr(q.size, eta, rng) @ q for _ in range(draws)])
