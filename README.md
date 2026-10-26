# h0-covariance-audit

Covariance- and provenance-sensitivity tools for **summary-level** Hubble-constant ($H_0$) compilations,
with a companion manuscript. The reference case is the 88-entry sound-horizon-free compilation of
Pantos & Perivolaropoulos (Phys. Rev. D 114, 063507; arXiv:2601.00650).

> **Scope.** This is a summary-level audit, not a joint-likelihood reanalysis. It does **not** show that the
> Hubble tension has disappeared, that its true significance equals any particular value, or that any
> measurement should be excluded. All $Z$ values are *nominal, conditional* statistics.

## What is implemented

| Result (see [`docs/theory.md`](docs/theory.md)) | C++ (`cpp/include/h0audit`) | Python (`python/src/h0audit`) |
|---|---|---|
| IVW and equicorrelated GLS, $O(n)$ closed form (Sherman–Morrison) | `gls.hpp` | `core.py` |
| Exact threshold for negative GLS weights | `negative_weight_threshold` | same |
| Exact PSD limit on a cross-group correlation $\lvert\rho_{LR}\rvert$ | `cross_rho_limit` | same |
| Response of $Z$ to $\kappa=\mathrm{Corr}(\hat H_L,\hat H_R)$ | `contrast_kappa` | same |
| Sharp worst-case variance under an entrywise correlation cap | `capped_variance_max` | same |
| LKJ($\eta$) prior-predictive Monte Carlo of $Z$ (multithreaded, deterministic) | `lkj.hpp` | `lkj.py` (slow reference) |

## Quick start

```bash
make cpp          # builds build/h0audit and build/h0audit_tests (g++ >= 9, C++17; no dependencies)
make test         # C++ unit tests + Python unit tests + C++/Python cross-check
make data         # needs network: downloads the pinned upstream table, verifies its Git blob SHA and 30/58 counts
build/h0audit summary data/h0_compilation.csv --rho-l 0.40 --rho-r 0.15
build/h0audit cap     data/h0_compilation.csv --r 0.10
build/h0audit lkj     data/h0_compilation.csv --eta 2 --draws 200000 --seed 1 --structure joint
```

Python only: `pip install -e . && python -m h0audit.cli summary data/h0_compilation.csv`.
`data/h0_compilation_reconstructed.csv` is a **provisional reconstruction** of the 88 rows from the source paper's appendix (ladder exact; comparison group within 0.004; 6 rows with an estimated error component; see `docs/reconstruction.md`). `data/synthetic_fixture.csv` is **synthetic** (clearly labelled) and exists only so the test-suite and
benchmarks run offline. Never quote numbers computed from it as physical results.

## Data format
`id,group,h0,sigma[,label]` with `group ∈ {L, R}` (L = distance ladder, R = comparison); lines starting with `#` are ignored.

## Repository layout
```
cpp/            header-only library, CLI, unit tests
python/         NumPy mirror, LKJ reference, unit tests
scripts/        fetch_compilation.py, make_figures.py, scripts/legacy/ (original audit scripts)
results/legacy_v1/  tables produced by the original audit scripts (kept verbatim)
provenance/     source-audit and dependency-edge tables
paper/          manuscript.tex / .pdf, references.bib, figures, cover-letter template
docs/           theory, methodology, reproducibility, limitations, API, pre-submission checklist
```

## Where each number in the manuscript comes from
See [`docs/reproducibility.md`](docs/reproducibility.md). In short: Tables 1–3 and Fig. 1 come from `results/legacy_v1/*.csv`;
Fig. 2, the $\kappa$ envelope, $\rho_{LR}$ limit and $n_{\rm eff}$ values are closed-form and covered by tests;
Section 5 and Tables 4–5 use the provisional reconstruction (`make reconstruct`); regenerate them from the authoritative table (`make data`) before citing.

## Citation and licence
See `CITATION.cff`. MIT licence. Please also cite the source compilation (Pantos & Perivolaropoulos 2026).
