// h0audit/lkj.hpp -- LKJ(eta) random correlation matrices (C-vine construction of
// Lewandowski, Kurowicka & Joe 2009) and a deterministic, thread-count-independent
// Monte Carlo driver for the contrast variance q^T R q.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

#include "linalg.hpp"

namespace h0 {

inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

template <class Rng>
inline double rand_beta(Rng& g, double a, double b) {
    std::gamma_distribution<double> ga(a, 1.0), gb(b, 1.0);
    const double x = ga(g), y = gb(g);
    return x / (x + y);
}

/// Draw R ~ LKJ(eta) of dimension d >= 2 (eta > 0). eta = 1 is uniform over the elliptope;
/// larger eta concentrates mass near the identity. Marginal law of R_ij:
/// 2*Beta(eta-1+d/2, eta-1+d/2) - 1.
template <class Rng>
inline Mat lkj_corr(std::size_t d, double eta, Rng& g) {
    Mat P(d), S(d);
    for (std::size_t i = 0; i < d; ++i) S(i, i) = 1.0;
    double beta = eta + 0.5 * static_cast<double>(d - 1);
    for (std::size_t k = 0; k + 1 < d; ++k) {
        beta -= 0.5;
        for (std::size_t i = k + 1; i < d; ++i) {
            double p = 2.0 * rand_beta(g, beta, beta) - 1.0;
            P(k, i) = p;
            for (std::size_t l = k; l-- > 0;)
                p = p * std::sqrt((1.0 - P(l, i) * P(l, i)) * (1.0 - P(l, k) * P(l, k))) + P(l, i) * P(l, k);
            S(k, i) = S(i, k) = p;
        }
    }
    return S;
}

inline double quad_form(const Mat& R, const Vec& q, std::size_t off = 0) {
    const std::size_t d = R.n;
    double v = 0.0;
    for (std::size_t i = 0; i < d; ++i) {
        v += q[off + i] * q[off + i];
        for (std::size_t j = i + 1; j < d; ++j) v += 2.0 * q[off + i] * q[off + j] * R(i, j);
    }
    return v;
}

enum class Structure { Joint, BlockDiagonal };

/// Monte Carlo draws of Var(Delta) = q^T R q, q = (a_L s_L, -a_R s_R).
/// Draws are grouped in fixed-size chunks, each seeded from (seed, chunk index), so the
/// output is bit-identical for any thread count.
inline Vec mc_variances(const Vec& q, std::size_t nL, double eta, Structure st, std::size_t draws,
                        std::uint64_t seed, unsigned threads) {
    const std::size_t d = q.size(), nR = d - nL;
    constexpr std::size_t kChunk = 256;
    const std::size_t nchunks = (draws + kChunk - 1) / kChunk;
    Vec out(draws);
    if (threads == 0) threads = std::max(1u, std::thread::hardware_concurrency());
    auto work = [&](unsigned tid) {
        for (std::size_t c = tid; c < nchunks; c += threads) {
            std::mt19937_64 g(splitmix64(seed ^ splitmix64(c + 1)));
            const std::size_t lo = c * kChunk, hi = std::min(draws, lo + kChunk);
            for (std::size_t m = lo; m < hi; ++m) {
                if (st == Structure::Joint) {
                    out[m] = quad_form(lkj_corr(d, eta, g), q);
                } else {
                    double v = 0.0;
                    if (nL >= 2) v += quad_form(lkj_corr(nL, eta, g), q, 0); else v += q[0] * q[0];
                    if (nR >= 2) v += quad_form(lkj_corr(nR, eta, g), q, nL); else v += q[nL] * q[nL];
                    out[m] = v;
                }
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    return out;
}

inline double quantile_sorted(const Vec& v, double p) {
    if (v.empty()) return std::nan("");
    const double x = p * static_cast<double>(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(x));
    const double f = x - static_cast<double>(i);
    return i + 1 < v.size() ? v[i] * (1.0 - f) + v[i + 1] * f : v[i];
}

}  // namespace h0
