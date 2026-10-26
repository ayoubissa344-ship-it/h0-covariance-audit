// h0audit -- command-line driver. All output is JSON on stdout.
//   h0audit summary DATA.csv [--rho-l 0.40] [--rho-r 0.15]
//   h0audit kappa   DATA.csv [--rho-l 0.40] [--rho-r 0.15]
//   h0audit cap     DATA.csv --r 0.10
//   h0audit lkj     DATA.csv [--eta 2] [--draws 100000] [--seed 1] [--threads 0] [--structure joint|block]
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "h0audit/gls.hpp"
#include "h0audit/io.hpp"
#include "h0audit/lkj.hpp"

using namespace h0;

static std::map<std::string, std::string> parse_flags(int argc, char** argv, int start) {
    std::map<std::string, std::string> m;
    for (int i = start; i + 1 < argc; i += 2) {
        std::string k = argv[i];
        if (k.rfind("--", 0) != 0) throw std::invalid_argument("expected --flag, got " + k);
        m[k.substr(2)] = argv[i + 1];
    }
    return m;
}
static double getd(const std::map<std::string, std::string>& m, const char* k, double def) {
    auto it = m.find(k);
    return it == m.end() ? def : std::stod(it->second);
}
static void jnum(std::ostream& o, const char* k, double v, bool last = false) {
    o << "\"" << k << "\": " << std::setprecision(12) << v << (last ? "" : ", ");
}
static Vec fixed_ivw_q(const Dataset& d, double& delta) {
    const Combined L = ivw(d.yL, d.sL), R = ivw(d.yR, d.sR);
    delta = L.mean - R.mean;
    Vec q;
    for (std::size_t i = 0; i < d.yL.size(); ++i) q.push_back(L.weights[i] * d.sL[i]);
    for (std::size_t j = 0; j < d.yR.size(); ++j) q.push_back(-R.weights[j] * d.sR[j]);
    return q;
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) { std::cerr << "usage: h0audit {summary|kappa|cap|lkj} DATA.csv [--flag value ...]\n"; return 2; }
        const std::string cmd = argv[1];
        const Dataset d = read_dataset(argv[2]);
        const auto fl = parse_flags(argc, argv, 3);
        std::ostringstream o;
        o << "{";
        if (cmd == "summary" || cmd == "kappa") {
            const double rl = getd(fl, "rho-l", 0.40), rr = getd(fl, "rho-r", 0.15);
            const Combined Li = ivw(d.yL, d.sL), Ri = ivw(d.yR, d.sR);
            const Combined Lg = gls_equicorr(d.yL, d.sL, rl), Rg = gls_equicorr(d.yR, d.sR, rr);
            const Contrast ci = contrast(Li, Ri), cg = contrast(Lg, Rg);
            const double lim = cross_rho_limit(d.yL.size(), rl, d.yR.size(), rr);
            o << "\"nL\": " << d.yL.size() << ", \"nR\": " << d.yR.size() << ", ";
            jnum(o, "ivw_L", Li.mean); jnum(o, "ivw_L_se", Li.se); jnum(o, "ivw_R", Ri.mean); jnum(o, "ivw_R_se", Ri.se);
            jnum(o, "ivw_delta", ci.delta); jnum(o, "ivw_se_delta", ci.se); jnum(o, "ivw_z", ci.z);
            jnum(o, "gls_L", Lg.mean); jnum(o, "gls_L_se", Lg.se); jnum(o, "gls_R", Rg.mean); jnum(o, "gls_R_se", Rg.se);
            jnum(o, "gls_delta", cg.delta); jnum(o, "gls_se_delta", cg.se); jnum(o, "gls_z", cg.z);
            o << "\"neg_weights_L\": " << Lg.n_negative << ", \"neg_weights_R\": " << Rg.n_negative << ", ";
            jnum(o, "neg_threshold_L", negative_weight_threshold(d.sL, rl));
            jnum(o, "neg_threshold_R", negative_weight_threshold(d.sR, rr));
            jnum(o, "cross_rho_limit", lim);
            // Feasible cross-group correlation grid (exact PSD limit) -> effect on Z.
            o << "\"rho_lr_scan\": [";
            const int K = 10;
            for (int k = -K; k <= K; ++k) {
                const double rlr = lim * k / K;
                const double cov = cross_cov(Lg.weights, d.sL, Rg.weights, d.sR, rlr);
                const Contrast c = contrast(Lg, Rg, cov);
                o << (k == -K ? "" : ", ") << "{";
                jnum(o, "rho_lr", rlr); jnum(o, "kappa", cov / (Lg.se * Rg.se)); jnum(o, "z", c.z, true);
                o << "}";
            }
            o << "], ";
            o << "\"kappa_envelope\": {";
            jnum(o, "z_at_kappa_minus1", contrast_kappa(Lg, Rg, -1.0).z);
            jnum(o, "z_at_kappa_plus1", contrast_kappa(Lg, Rg, 1.0).z, true);
            o << "}";
        } else if (cmd == "cap") {
            const double r = getd(fl, "r", 0.1);
            double delta;
            const Vec q = fixed_ivw_q(d, delta);
            const double v = capped_variance_max(q, r);
            jnum(o, "r", r); jnum(o, "delta", delta); jnum(o, "var_max", v); jnum(o, "se_max", std::sqrt(v));
            jnum(o, "z_min", std::fabs(delta) / std::sqrt(v), true);
        } else if (cmd == "lkj") {
            const double eta = getd(fl, "eta", 2.0);
            const std::size_t draws = static_cast<std::size_t>(getd(fl, "draws", 100000));
            const std::uint64_t seed = static_cast<std::uint64_t>(getd(fl, "seed", 1));
            const unsigned thr = static_cast<unsigned>(getd(fl, "threads", 0));
            const std::string st = fl.count("structure") ? fl.at("structure") : "joint";
            double delta;
            const Vec q = fixed_ivw_q(d, delta);
            Vec v = mc_variances(q, d.yL.size(), eta, st == "block" ? Structure::BlockDiagonal : Structure::Joint, draws, seed, thr);
            Vec z(v.size());
            std::size_t g2 = 0, g3 = 0, g5 = 0;
            double mv = 0.0;
            for (std::size_t i = 0; i < v.size(); ++i) {
                z[i] = std::fabs(delta) / std::sqrt(std::max(v[i], 1e-300));
                g2 += z[i] > 2; g3 += z[i] > 3; g5 += z[i] > 5; mv += v[i];
            }
            std::sort(z.begin(), z.end());
            const double n = static_cast<double>(z.size());
            o << "\"structure\": \"" << st << "\", \"draws\": " << draws << ", \"seed\": " << seed << ", ";
            jnum(o, "eta", eta); jnum(o, "delta", delta); jnum(o, "mean_var", mv / n);
            jnum(o, "z_q025", quantile_sorted(z, 0.025)); jnum(o, "z_q16", quantile_sorted(z, 0.16));
            jnum(o, "z_median", quantile_sorted(z, 0.5)); jnum(o, "z_q84", quantile_sorted(z, 0.84));
            jnum(o, "z_q975", quantile_sorted(z, 0.975));
            jnum(o, "p_z_gt_2", g2 / n); jnum(o, "p_z_gt_3", g3 / n); jnum(o, "p_z_gt_5", g5 / n, true);
        } else {
            std::cerr << "unknown command " << cmd << "\n";
            return 2;
        }
        o << "}\n";
        std::cout << o.str();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
