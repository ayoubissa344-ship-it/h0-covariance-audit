# API summary
**C++** (`#include "h0audit/gls.hpp"`, namespace `h0`): `ivw`, `gls_equicorr`, `gls_dense`, `negative_weight_threshold`, `cross_rho_limit`, `cross_cov`, `contrast`, `contrast_kappa`, `capped_variance_max`;
(`lkj.hpp`) `lkj_corr`, `mc_variances`, `quantile_sorted`; (`io.hpp`) `read_dataset`.
**Python** (`import h0audit`): the same names in snake_case plus `information_ratio`, `block_equicorr_matrix`, `read_dataset`, `parse_latex_table`.
**CLI:** `h0audit {summary|kappa|cap|lkj} DATA.csv [--flag value]`; JSON on stdout. `lkj` flags: `--eta --draws --seed --threads --structure joint|block`.
All functions validate inputs (finite values, $s_i>0$, admissible $\rho$) and throw/raise on violation.
