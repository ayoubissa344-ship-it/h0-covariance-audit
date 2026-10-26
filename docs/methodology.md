# Methodology
1. **Freeze the table.** `make data` downloads the pinned upstream file, checks its Git blob SHA (`9c335e27…`), checks 30/58 rows, and writes `data/h0_compilation.csv`.
2. **Symmetrise** asymmetric errors by RMS of the upper/lower errors (upstream convention; discards non-Gaussian information).
3. **Reproduce** IVW and equicorrelated GLS (`h0audit summary`), compare to `results/legacy_v1/hubble_tension_phase7_scenarios.csv`.
4. **Vary assumptions explicitly:** estimator (`results/legacy_v1/hubble_tension_estimator_sensitivity.csv`), cross-group dependence (`rho_lr_scan`, κ envelope), entrywise caps (`cap`), random correlation structures (`lkj`).
5. **Report conditionally.** Every $Z$ is stated together with the estimator, within-group $\rho$, cross-group covariance and cap/prior it assumes.
6. **Never delete rows to move the answer.** Exclusion scenarios are influence diagnostics only; $\rho$ is not re-tuned after exclusion.
7. **Provenance** (`provenance/*.csv`): each flagged row has a source layer, version, public pointer and an audit observation; overlap is recorded as a graph edge and only quantified if joint covariance products support it.
