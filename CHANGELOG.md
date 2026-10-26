# Changelog
## 1.0.0 — 2026-10-25
- Header-only C++17 library (`cpp/include/h0audit`): IVW, closed-form equicorrelated GLS, dense Cholesky GLS, negative-weight threshold, exact cross-group PSD limit, κ-response, sharp entrywise-cap bound, LKJ(η) prior-predictive Monte Carlo (deterministic, multithreaded).
- NumPy mirror (`python/src/h0audit`) cross-validated against the C++ CLI.
- Provisional offline reconstruction of the 88-row table (`data/h0_compilation_reconstructed.csv`, per-row status) with validation against packaged results; LKJ moment theorem; results in `results/reconstructed/`.
- Manuscript (LaTeX + PDF) rewritten around the proved results; legacy tables preserved under `results/legacy_v1/`.
