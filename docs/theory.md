# Theory

Notation: group of $n$ estimates $y_i$ with errors $s_i$; $a_i=1/s_i$, $A=\sum a_i$; $C=SRS$, $S=\mathrm{diag}(s)$;
$R_\rho=(1-\rho)I+\rho\mathbf 1\mathbf 1^\top$ with $\rho\in(-1/(n-1),1)$; $c=\rho/(1+(n-1)\rho)$.

## T1. Closed-form equicorrelated GLS
$(C^{-1}\mathbf 1)_i=a_i(a_i-cA)/(1-\rho)$, $N=\mathbf 1^\top C^{-1}\mathbf 1=(\sum a_i^2-cA^2)/(1-\rho)>0$, $w_i=(C^{-1}\mathbf 1)_i/N$, $\mathrm{Var}=1/N$.
*Proof.* Sherman–Morrison: $R_\rho^{-1}=(I-c\mathbf 1\mathbf 1^\top)/(1-\rho)$; $C^{-1}=S^{-1}R_\rho^{-1}S^{-1}$.
$N>0$: for $\rho>0$, $A^2\le n\sum a_i^2$ and $cn<1$. ∎  Test: `closed_form_matches_dense`.

## T2. Sign of the weights
For $\rho>0$: $w_i<0 \iff s_i>s^\*=(1+(n-1)\rho)/(\rho\sum_j s_j^{-1})$. No negative weights if all $s_i$ are equal.

## T3. Design effect
Equal errors: $\mathrm{Var}=s^2(1+(n-1)\rho)/n$; $n_{\rm eff}=n/(1+(n-1)\rho)\to1/\rho$. (Idealisation; it is a consequence of the *assumed* $\rho$.)
For unequal errors the ratio $\mathrm{Var}_{\rm indep}/\mathrm{Var}_{\rm GLS}$ (`information_ratio`) can exceed 1 — correlated errors can be partly cancelled.

## T4. Exact PSD limit on a cross-group correlation
Equicorrelated blocks, constant cross correlation $\rho_{LR}$: PSD iff $\rho_g\in[-1/(n_g-1),1]$ and
$\rho_{LR}^2\le[1+(n_L-1)\rho_L][1+(n_R-1)\rho_R]/(n_Ln_R)$.
*Proof.* Block-orthogonal-to-$\mathbf 1$ eigenvalues are $1-\rho_g$; the invariant plane gives a $2\times2$ matrix whose determinant must be $\ge0$. ∎
For $(30,0.40,58,0.15)$: $0.2630$. Pooled-mean covariance: $\mathrm{Cov}(\hat H_L,\hat H_R)=\rho_{LR}(w_L\!\cdot s_L)(w_R\!\cdot s_R)$.

## T5. κ-response
$Z(\kappa)=|\Delta|/\sqrt{SE_L^2+SE_R^2-2\kappa SE_LSE_R}$ is strictly increasing in $\kappa\in[-1,1)$.
Envelope $[\,|\Delta|/(SE_L+SE_R),\ |\Delta|/|SE_L-SE_R|\,]$. For the packaged GLS means: $[2.805,\ 18.43]$, $Z(0)=3.9215$.
Consequence: setting cross-group covariance to zero is *not* conservative in general.

## T6. Sharp entrywise-cap bound
$\max\{q^\top Rq:\ R\succeq0,\ R_{ii}=1,\ |R_{ij}|\le r\}=(1-r)\sum q_i^2+r(\sum|q_i|)^2$, attained at $R=(1-r)I+r\,\sigma\sigma^\top$, $\sigma=\mathrm{sign}(q)$.
*Proof.* Upper bound termwise; $R$ is PSD for $r\in[0,1]$. ∎

## T7. Identifiability
$y_i=H_\*+\alpha_{g(i)}+\varepsilon_i$ has design rank 2 for 3 parameters; $c^\top\theta$ is estimable iff $c^\top\in\mathrm{row}(X)$.
$\alpha_L-\alpha_R$ is estimable; $H_\*$ and the offsets are not without an anchor. With zero cross-covariance the BLUE of the contrast is
$\hat H_L^{\rm GLS}-\hat H_R^{\rm GLS}$; with non-zero cross-covariance the joint Aitken estimator is BLUE and the difference of group fits is not.

## T8. Factor representation of $\rho$
$R_\rho$ is the correlation of $\varepsilon_i=s_i(\sqrt\rho\,\eta+\sqrt{1-\rho}\,\xi_i)$: $V=Z\Omega Z^\top+D$ with $Z=s$, $\Omega=\rho$, $D=(1-\rho)S^2$.
Documented dependency models replace $Z$ by a 0/1 design of shared datasets/calibrations and constrain $\Omega\succeq0$ from likelihood products.

## T9. Moments of the contrast variance under LKJ
$E[q^\top Rq]=\sum q_i^2$, $\mathrm{Var}(q^\top Rq)=2[(\sum q_i^2)^2-\sum q_i^4]/(d+2\eta-1)$.
*Proof.* LKJ is invariant under $R\mapsto DRD$, $D=\mathrm{diag}(\pm1)$, so distinct off-diagonal pairs are uncorrelated; each has variance $1/(d+2\eta-1)$. ∎
Hence random (zero-mean) dependence averages out: relative sd of the variance ≤ $\sqrt{2/(d+2\eta-1)}$ (≈15 % for $d=88,\eta=2$). Test: `test_lkj_quadratic_form_moments`.

## LKJ prior-predictive tool
$R\sim\mathrm{LKJ}(\eta)$ by the C-vine construction (Lewandowski–Kurowicka–Joe 2009). Marginal $R_{ij}\sim2\,\mathrm{Beta}(\eta-1+d/2,\eta-1+d/2)-1$, variance $1/(2\eta-1+d)$ (tested).
$E[R]=I$, so $E[q^\top Rq]=\sum q_i^2$. The output is the distribution of $Z=|\Delta|/\sqrt{q^\top Rq}$ for **fixed IVW weights**; it explores sensitivity to unknown correlations and is not a physical model.
