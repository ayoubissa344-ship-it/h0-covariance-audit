// h0audit/gls.hpp -- inverse-variance and generalised least squares combination of
// summary-level estimates, with closed forms for the equicorrelation model.
//
// Notation (see docs/theory.md):  y_i estimates, s_i standard errors,
//   C = S R S, R_ii = 1, R_ij = rho (i != j).
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "linalg.hpp"

namespace h0 {

struct Combined {
    double mean = 0.0;
    double se = 0.0;
    Vec weights;              // normalised, sum to 1 (may be negative under GLS)
    std::size_t n_negative = 0;
};

inline void check_inputs(const Vec& y, const Vec& s) {
    if (y.size() != s.size() || y.empty()) throw std::invalid_argument("y/s size mismatch or empty");
    for (std::size_t i = 0; i < y.size(); ++i) {
        if (!std::isfinite(y[i]) || !std::isfinite(s[i]) || !(s[i] > 0.0))
            throw std::invalid_argument("non-finite estimate or non-positive sigma");
    }
}

/// Diagonal inverse-variance weighted mean (independence assumption).
inline Combined ivw(const Vec& y, const Vec& s) {
    check_inputs(y, s);
    Combined r;
    double tot = 0.0;
    r.weights.resize(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) { r.weights[i] = 1.0 / (s[i] * s[i]); tot += r.weights[i]; }
    for (double& w : r.weights) w /= tot;
    for (std::size_t i = 0; i < y.size(); ++i) r.mean += r.weights[i] * y[i];
    r.se = std::sqrt(1.0 / tot);
    return r;
}

/// Admissible interval for a common correlation among n >= 2 variables: (-1/(n-1), 1).
inline double rho_lower_limit(std::size_t n) { return n < 2 ? -1.0 : -1.0 / static_cast<double>(n - 1); }

/// Equicorrelation GLS in O(n) via the Sherman-Morrison identity:
///   (C^{-1} 1)_i = (1/(1-rho)) * ( a_i^2 - c a_i A ),  a_i = 1/s_i, A = sum a_j,
///   c = rho / (1 + (n-1) rho).
inline Combined gls_equicorr(const Vec& y, const Vec& s, double rho) {
    check_inputs(y, s);
    const std::size_t n = y.size();
    if (!(rho > rho_lower_limit(n) && rho < 1.0)) throw std::invalid_argument("rho outside (-1/(n-1), 1)");
    const double c = rho / (1.0 + static_cast<double>(n - 1) * rho);
    double A = 0.0;
    for (double si : s) A += 1.0 / si;
    Vec u(n);
    double N = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double a = 1.0 / s[i];
        u[i] = (a * a - c * a * A) / (1.0 - rho);
        N += u[i];
    }
    if (!(N > 0.0)) throw std::runtime_error("non-positive GLS normalisation");
    Combined r;
    r.weights.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        r.weights[i] = u[i] / N;
        r.mean += r.weights[i] * y[i];
        if (r.weights[i] < -1e-12) ++r.n_negative;
    }
    r.se = std::sqrt(1.0 / N);
    return r;
}

/// Row uncertainty above which an equicorrelation-GLS weight is negative (rho > 0):
///   s* = (1 + (n-1) rho) / (rho * sum_j 1/s_j).
inline double negative_weight_threshold(const Vec& s, double rho) {
    if (!(rho > 0.0)) return std::numeric_limits<double>::infinity();
    double A = 0.0;
    for (double si : s) A += 1.0 / si;
    return (1.0 + static_cast<double>(s.size() - 1) * rho) / (rho * A);
}

/// Dense GLS for an arbitrary SPD covariance C (Cholesky, O(n^3)).
inline Combined gls_dense(const Vec& y, const Mat& C) {
    if (C.n != y.size() || y.empty()) throw std::invalid_argument("gls_dense: size mismatch");
    Mat L;
    if (!cholesky(C, L)) throw std::runtime_error("covariance matrix is not positive definite");
    Vec ones(y.size(), 1.0);
    Vec u = chol_solve(L, ones);
    double N = 0.0;
    for (double v : u) N += v;
    if (!(N > 0.0)) throw std::runtime_error("non-positive GLS normalisation");
    Combined r;
    r.weights.resize(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        r.weights[i] = u[i] / N;
        r.mean += r.weights[i] * y[i];
        if (r.weights[i] < -1e-12) ++r.n_negative;
    }
    r.se = std::sqrt(1.0 / N);
    return r;
}

/// Exact PSD limit on a constant cross-group correlation rho_LR when the two blocks are
/// equicorrelated with rho_L, rho_R (sizes nL, nR):
///   rho_LR^2 <= (1+(nL-1)rho_L)(1+(nR-1)rho_R) / (nL nR).
inline double cross_rho_limit(std::size_t nL, double rhoL, std::size_t nR, double rhoR) {
    const double a = 1.0 + static_cast<double>(nL - 1) * rhoL;
    const double b = 1.0 + static_cast<double>(nR - 1) * rhoR;
    if (!(a > 0.0 && b > 0.0 && rhoL < 1.0 && rhoR < 1.0)) throw std::invalid_argument("blocks not PSD");
    return std::sqrt(a * b / (static_cast<double>(nL) * static_cast<double>(nR)));
}

/// Cov(H_L, H_R) for linear estimators with weight vectors wL, wR when the cross-block
/// covariance is rho_LR * s_L s_R^T.
inline double cross_cov(const Vec& wL, const Vec& sL, const Vec& wR, const Vec& sR, double rhoLR) {
    return rhoLR * dot(wL, sL) * dot(wR, sR);
}

struct Contrast {
    double delta = 0.0, var = 0.0, se = 0.0, z = 0.0;
};

/// Contrast Delta = H_L - H_R with Var = VL + VR - 2 Cov.
inline Contrast contrast(const Combined& L, const Combined& R, double cov_LR = 0.0) {
    Contrast c;
    c.delta = L.mean - R.mean;
    c.var = L.se * L.se + R.se * R.se - 2.0 * cov_LR;
    if (!(c.var > 0.0)) throw std::runtime_error("non-positive contrast variance (covariance infeasible)");
    c.se = std::sqrt(c.var);
    c.z = std::fabs(c.delta) / c.se;
    return c;
}

/// Contrast under a correlation kappa = Corr(H_L, H_R) in [-1, 1] between the pooled means.
inline Contrast contrast_kappa(const Combined& L, const Combined& R, double kappa) {
    return contrast(L, R, kappa * L.se * R.se);
}

/// Sharp worst-case (largest) variance of q^T R q over correlation matrices with
/// |R_ij| <= r for all i != j:  (1-r) sum q_i^2 + r (sum |q_i|)^2  (attained by a rank-one sign pattern).
inline double capped_variance_max(const Vec& q, double r) {
    if (!(r >= 0.0 && r <= 1.0)) throw std::invalid_argument("r must lie in [0,1]");
    double q2 = 0.0, q1 = 0.0;
    for (double x : q) { q2 += x * x; q1 += std::fabs(x); }
    return (1.0 - r) * q2 + r * q1 * q1;
}

}  // namespace h0
