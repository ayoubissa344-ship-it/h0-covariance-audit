// Self-contained unit tests (no framework). Exit code 0 iff all pass.
#include <cmath>
#include <iostream>
#include <random>

#include "h0audit/gls.hpp"
#include "h0audit/lkj.hpp"

using namespace h0;
static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { ++failures; std::cerr << "FAIL: " << msg << "  [" #cond "]\n"; } } while (0)
static bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(b)); }

static void make_data(std::size_t n, std::uint64_t seed, Vec& y, Vec& s) {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> us(0.4, 3.0);
    std::normal_distribution<double> nd(70.0, 2.0);
    y.resize(n); s.resize(n);
    for (std::size_t i = 0; i < n; ++i) { s[i] = us(g); y[i] = nd(g); }
}

int main() {
    // 1. Closed-form equicorrelation GLS == dense Cholesky GLS.
    for (double rho : {0.0, 0.15, 0.40, 0.9, -0.01}) {
        Vec y, s; make_data(30, 7, y, s);
        Mat C(30);
        for (std::size_t i = 0; i < 30; ++i) for (std::size_t j = 0; j < 30; ++j) C(i, j) = s[i] * s[j] * (i == j ? 1.0 : rho);
        const Combined a = gls_equicorr(y, s, rho), b = gls_dense(y, C);
        CHECK(close(a.mean, b.mean, 1e-9), "mean mismatch rho=" << rho);
        CHECK(close(a.se, b.se, 1e-9), "se mismatch rho=" << rho);
        CHECK(a.n_negative == b.n_negative, "negative-weight count mismatch");
    }
    // 2. rho = 0 reduces to IVW.
    { Vec y, s; make_data(58, 3, y, s); const Combined a = ivw(y, s), b = gls_equicorr(y, s, 0.0);
      CHECK(close(a.mean, b.mean, 1e-12) && close(a.se, b.se, 1e-12), "rho=0 != IVW"); }
    // 3. Negative-weight threshold is exact.
    { Vec y, s; make_data(58, 11, y, s); const double rho = 0.15; const Combined g = gls_equicorr(y, s, rho);
      const double th = negative_weight_threshold(s, rho);
      for (std::size_t i = 0; i < s.size(); ++i) {
          if (std::fabs(s[i] - th) < 1e-9) continue;
          CHECK((g.weights[i] < 0.0) == (s[i] > th), "threshold rule violated at i=" << i);
      } }
    // 4. Equal-sigma design effect: Var = s^2 (1+(n-1)rho)/n.
    { Vec y(40, 70.0), s(40, 1.3); const double rho = 0.2; const Combined g = gls_equicorr(y, s, rho);
      CHECK(close(g.se * g.se, 1.3 * 1.3 * (1.0 + 39 * rho) / 40.0, 1e-12), "design effect"); }
    // 5. Cross-group PSD limit is tight (Cholesky succeeds just inside, fails just outside).
    { const std::size_t nL = 12, nR = 20; const double rl = 0.4, rr = 0.15;
      const double lim = cross_rho_limit(nL, rl, nR, rr);
      auto build = [&](double rlr) { Mat R(nL + nR);
        for (std::size_t i = 0; i < nL + nR; ++i) for (std::size_t j = 0; j < nL + nR; ++j) {
            const bool li = i < nL, lj = j < nL;
            R(i, j) = i == j ? 1.0 : (li && lj ? rl : (!li && !lj ? rr : rlr)); }
        return R; };
      Mat L; CHECK(cholesky(build(0.999 * lim), L, 1e-14), "just inside limit should be PD");
      CHECK(!cholesky(build(1.001 * lim), L, 1e-14), "just outside limit should fail");
      CHECK(cholesky(build(-0.999 * lim), L, 1e-14), "negative side inside");
      CHECK(!cholesky(build(-1.001 * lim), L, 1e-14), "negative side outside"); }
    // 6. Capped-variance bound: attained by rank-one sign matrix, never exceeded by random PSD with |R_ij|<=r.
    { Vec q = {0.3, -0.5, 0.2, -0.1, 0.4}; const double r = 0.25; const double vmax = capped_variance_max(q, r);
      Mat R(5); for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { const double si = q[i] >= 0 ? 1 : -1, sj = q[j] >= 0 ? 1 : -1; R(i, j) = i == j ? 1.0 : r * si * sj; }
      CHECK(close(quad_form(R, q), vmax, 1e-12), "bound not attained");
      std::mt19937_64 g(5); std::uniform_real_distribution<double> u(-r, r); double worst = 0;
      for (int t = 0; t < 20000; ++t) { Mat A(5); for (int i = 0; i < 5; ++i) { A(i, i) = 1; for (int j = 0; j < i; ++j) A(i, j) = A(j, i) = u(g); }
        Mat L; if (!cholesky(A, L)) continue; worst = std::max(worst, quad_form(A, q)); }
      CHECK(worst <= vmax + 1e-12, "random feasible matrix exceeded bound"); }
    // 7. LKJ draws are valid correlation matrices; marginal variance of R_ij is 1/(2 eta - 1 + d).
    { std::mt19937_64 g(99); const std::size_t d = 10; const double eta = 1.0; double s2 = 0, s1 = 0; const int M = 20000; int pd = 0;
      for (int t = 0; t < M; ++t) { const Mat R = lkj_corr(d, eta, g); Mat L; pd += cholesky(R, L, 1e-14); s1 += R(2, 7); s2 += R(2, 7) * R(2, 7); }
      CHECK(pd == M, "LKJ draw not PD");
      CHECK(std::fabs(s1 / M) < 0.01, "LKJ mean != 0");
      CHECK(std::fabs(s2 / M - 1.0 / (2 * eta - 1 + d)) < 0.004, "LKJ marginal variance off: " << s2 / M); }
    // 8. Monte Carlo is bit-identical for any thread count; E[q^T R q] = sum q^2 (E R = I).
    { Vec q; std::mt19937_64 g(1); std::normal_distribution<double> nd(0, 0.1); for (int i = 0; i < 30; ++i) q.push_back(nd(g));
      const Vec a = mc_variances(q, 12, 2.0, Structure::Joint, 3000, 42, 1), b = mc_variances(q, 12, 2.0, Structure::Joint, 3000, 42, 4);
      CHECK(a == b, "thread-count dependence");
      double m = 0, q2 = 0; for (double v : a) m += v; m /= a.size(); for (double x : q) q2 += x * x;
      CHECK(std::fabs(m - q2) < 0.1 * q2, "MC mean far from sum q^2"); }
    if (failures == 0) std::cout << "ALL C++ TESTS PASSED\n";
    return failures == 0 ? 0 : 1;
}
