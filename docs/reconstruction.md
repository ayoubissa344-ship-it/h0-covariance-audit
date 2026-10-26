# Provisional reconstruction of the 88-row table

**Why.** The pinned upstream commit could not be fetched where this repository was prepared (no network in the build sandbox; the fetch tool refused the GitHub URL).
**What.** `scripts/reconstruct_from_paper.py` rebuilds all 88 rows (30 ladder `L`, 58 comparison `R`) from Appendix A, Tables X–XI of arXiv:2601.00650v2, applying the paper's Eq. (4) RMS rule to asymmetric errors, and writes `data/h0_compilation_reconstructed.csv` with a per-row `status`:

| status | rows | meaning |
|---|---|---|
| exact | 55 | value and symmetric error as printed |
| symmetrised | 27 | RMS of two printed errors |
| approx | 6 | one error component unreadable in the extracted text and estimated: Abbott 2021 (lower), Gao 2023 (lower), Fung 2023 (lower), Bahr-Kalus 2025 ΛCDM (errors), Paic 2025 (lower), Jaiswal 2024 (lower) |

**Validation** (`python scripts/validate_reconstruction.py` → `results/reconstructed/validation.csv`): ladder IVW/GLS reproduced to 1e-9 and 16 negative weights;
comparison IVW mean differs by 0.0034, SE by 0.0003, GLS mean by 0.0009, 38 negative weights reproduced; cap bounds agree with `results/legacy_v1/hubble_tension_robust_bounds.csv` to <0.01.
Two entries were checked against their abstracts: Dominguez et al. 2023 (61.9 +2.9 −2.4) and Gonzalez et al. 2024 (72.7 +6.3 −5.6, PLB 857).

**Use.** Treat every number computed from this file as provisional. Replace it by the authoritative table (`make data`) and rerun `results/reconstructed/*` commands (see `docs/reproducibility.md`) before citing.
