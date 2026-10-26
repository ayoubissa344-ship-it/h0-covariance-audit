// h0audit/linalg.hpp -- minimal dense linear algebra (header-only, no dependencies).
#pragma once
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace h0 {

using Vec = std::vector<double>;

/// Row-major square matrix.
struct Mat {
    std::size_t n = 0;
    Vec a;
    Mat() = default;
    explicit Mat(std::size_t n_) : n(n_), a(n_ * n_, 0.0) {}
    double& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
    double operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }
};

inline double dot(const Vec& x, const Vec& y) {
    if (x.size() != y.size()) throw std::invalid_argument("dot: size mismatch");
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += x[i] * y[i];
    return s;
}

/// Cholesky factorisation A = L L^T (lower). Returns false if A is not numerically
/// positive definite (pivot <= tol). Only the lower triangle of A is read.
inline bool cholesky(const Mat& A, Mat& L, double tol = 1e-13) {
    const std::size_t n = A.n;
    L = Mat(n);
    for (std::size_t j = 0; j < n; ++j) {
        double d = A(j, j);
        for (std::size_t k = 0; k < j; ++k) d -= L(j, k) * L(j, k);
        if (!(d > tol)) return false;
        const double ljj = std::sqrt(d);
        L(j, j) = ljj;
        for (std::size_t i = j + 1; i < n; ++i) {
            double s = A(i, j);
            for (std::size_t k = 0; k < j; ++k) s -= L(i, k) * L(j, k);
            L(i, j) = s / ljj;
        }
    }
    return true;
}

/// Solve (L L^T) x = b given the Cholesky factor L.
inline Vec chol_solve(const Mat& L, const Vec& b) {
    const std::size_t n = L.n;
    Vec y(n), x(n);
    for (std::size_t i = 0; i < n; ++i) {
        double s = b[i];
        for (std::size_t k = 0; k < i; ++k) s -= L(i, k) * y[k];
        y[i] = s / L(i, i);
    }
    for (std::size_t ii = n; ii-- > 0;) {
        double s = y[ii];
        for (std::size_t k = ii + 1; k < n; ++k) s -= L(k, ii) * x[k];
        x[ii] = s / L(ii, ii);
    }
    return x;
}

}  // namespace h0
