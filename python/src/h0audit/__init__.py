"""h0audit: covariance-sensitivity tools for summary-level H0 compilations."""
from .core import (
    Combined, Contrast, ivw, gls_equicorr, gls_dense, negative_weight_threshold,
    cross_rho_limit, cross_cov, contrast, contrast_kappa, capped_variance_max,
    block_equicorr_matrix, information_ratio,
)
from .io import Dataset, read_dataset, parse_sigma, parse_latex_table

__version__ = "1.0.0"
__all__ = [n for n in dir() if not n.startswith("_")]
