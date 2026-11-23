# Causal Atlas Scientific Methods & Numerical Formulations

This document details the clean-room mathematical, statistical, and geometric algorithms implemented in `ca-core` and `ca-app/src/svg_export.rs`.

---

## 1. Over-Representation & Total-Hypothesis Benjamini-Hochberg FDR

### One-Sided Fisher's Exact Test (`ca_core::stats::fishers_exact_right_tail`)
Given gene universe size $N$, target set size $M$, active dataset DEG count $n$, and intersection $k = |T \cap D|$:
$$p = \sum_{j = k}^{\min(M, n)} \frac{\binom{M}{j}\binom{N - M}{n - j}}{\binom{N}{n}}$$
Log-factorials $\ln(m!)$ are tabulated exactly for $m \le 100{,}000$ and evaluated via Stirling's series above $100{,}000$, using log-sum-exp rescaling (`max_log`) to prevent floating-point underflow or overflow.

### Total-Hypothesis Benjamini-Hochberg Step-Up FDR (`benjamini_hochberg_with_total`)
In genome-wide pathway (`canonical.rs`), regulator (`ura.rs`), and phenotype (`dea.rs`) sweeps, candidate sets with insufficient overlap (e.g., `overlap == 0` in Canonical Pathways or `overlap_count < 2` in URA/DEA) represent tested null hypotheses whose enrichment cannot pass significance (`p = 1.0` for zero overlap).

Rather than applying Benjamini-Hochberg solely to the filtered subset of $K$ reported hits, **Causal Atlas** adjusts across the **full count of tested hypotheses $m \ge K$** (all non-empty pathways with `pw_size > 0` in Canonical Pathways; all molecular regulators with `regulon_size >= 2` in URA; all phenotype terms with `gene_map.len() >= 2` in DEA):
$$q_{(r)} = \min_{j \ge r} \min\!\left(1,\; \frac{m}{j}\, p_{(j)}\right), \qquad r = 1, \dots, K$$
Because the $m - K$ omitted null hypotheses satisfy $p = 1.0$, their BH adjusted values equal $1.0$ and act solely through the global multiple-testing burden $m$.

---

## 2. Signed Activation Z-Score, Unsigned Exclusion & Target Concordance

### Krämer et al. (2014) Weighted Activation Z-Score (`kramer_activation_z`)
For a regulator or pathway with signed interactions $(w_i, s^{\text{edge}}_i, s^{\text{data}}_i)$ where $w_i > 0$ is finite, $s^{\text{edge}}_i \in \{-1, +1\}$, and $s^{\text{data}}_i = \operatorname{sgn}(\log_2 \mathrm{FC}_i) \in \{-1, +1\}$:
$$x_i = s^{\text{edge}}_i \, s^{\text{data}}_i \in \{-1, +1\}, \qquad Z = \frac{\sum_{i \in S} w_i x_i}{\sqrt{\sum_{i \in S} w_i^2}}$$

- **Unsigned Member Exclusion**: Interactions with `EdgeSign::Ambiguous` (`s_edge == 0`), `expected_dir == 0`, or zero fold change (`s_data == 0`) are strictly excluded from $S$, so unsigned binding links neither dilute the denominator $\sqrt{\sum w_i^2}$ nor distort $\mu_1$ (`mean_pw_sign` / regulon mean sign).
- **Skewness Bias-Corrected Z-Score**: Defining $\mu_1$ as the mean expected sign across all *signed* regulon/pathway members in the knowledge graph and $\mu_2$ as the mean expression polarity across active DEGs, let $\mu = \operatorname{clamp}(\mu_1 \mu_2, -0.99, 0.99)$. When $|\mu| > 0.25$ (`skew_flagged = true`), the bias-corrected Z-score is:
  $$Z_{\text{corr}} = \frac{\sum_{i \in S} w_i (x_i - \mu)}{\sqrt{\max(1 - \mu^2, 0.05) \sum_{i \in S} w_i^2}}$$
- **Two-Decimal State Consistency**: `PredictedState::from_z` rounds `z` to two decimal places (`round_z_2dp`), classifying `Activated` ($\ge +2.00$), `Inhibited` ($\le -2.00$), or `Neutral` (`Unresolved`), matching displayed `{:+.2}` labels identically.

### URA Target Concordance Evaluated Against `z.signum()` (`ura.rs`)
For every measured target in an upstream regulator's `target_evidences`, concordance (`consistent`) is evaluated relative to the regulator's inferred activation direction $\operatorname{sgn}(Z)$ rather than assuming every regulator is activated:
$$\text{consistent}_i \iff (Z \ne 0) \;\wedge\; (s^{\text{edge}}_i \ne 0) \;\wedge\; (s^{\text{data}}_i \ne 0) \;\wedge\; \left(s^{\text{edge}}_i \, s^{\text{data}}_i \, \operatorname{sgn}(Z) > 0\right)$$
Consequently, when an upstream regulator is predicted **Inhibited** ($Z \le -2.0$), downregulated targets ($s^{\text{data}}_i = -1$) linked via activating edges ($s^{\text{edge}}_i = +1$)—as well as upregulated targets linked via repressing edges ($s^{\text{edge}}_i = -1$)—correctly score `consistent = true`. Downstream **Regulator Effects** (`dea.rs`) directly consumes `t.consistent` for all activated and inhibited regulators.

---

## 3. Complementary Error Function Continued-Fraction Tail (`stats.rs`)

`ca_core::stats::erfc_approx(x)` provides two-regime evaluation for `z_to_two_sided_p(Z) = erfc(|Z| / sqrt(2))`:

1. **Core Regime ($|x| \le 5$)**: Abramowitz & Stegun (7.1.26) 5-term rational polynomial in $t = (1 + 0.3275911 |x|)^{-1}$ multiplied by $e^{-x^2}$.
2. **Extreme Tail Regime ($|x| > 5$)**: Standard fixed-degree polynomial approximations suffer increasing relative error in extreme p-value tails ($|Z| > 7$). For $|x| > 5$, `erfc(|x|)` is evaluated via the incomplete-gamma continued fraction $\operatorname{erfc}(z) = \Gamma(\tfrac{1}{2}, z^2) / \sqrt{\pi}$:
   $$\operatorname{erfc}(z) = \frac{z\, e^{-z^2}}{\sqrt{\pi}} \cdot \cfrac{1}{z^2 + \frac{1}{2} - \cfrac{1 \cdot \frac{1}{2}}{z^2 + \frac{5}{2} - \cfrac{2 \cdot \frac{3}{2}}{z^2 + \frac{9}{2} - \cdots}}}$$
   computed using modified Lentz's method (`TINY = 1e-300`, `EPS = 1e-15`, up to 128 steps). For $z > 5$ ($z^2 > 25$), convergence completes in a small number of iterations to near IEEE-754 `f64` machine precision down through `erfc(20) ≈ 5.40e-176`.

---

## 4. Rank-Based Signed CMap / L1000 Connectivity (`perturb.rs`)

`run_perturbation_match` guards against `NaN` / `±Inf` fold-changes by filtering to finite `(symbol, log2_fc)` measurements, deduplicating by strongest magnitude, and sorting descending by `log2_fc` into ranks $r \in \{1, \dots, N\}$.

For a signature gene set $T$ of size $t = |T|$ with sorted unique 1-based ranks $V(1) < V(2) < \dots < V(t)$, the unweighted two-sample Kolmogorov-Smirnov / Subramanian running-sum enrichment score $ES(T) \in [-1, +1]$ is evaluated in $O(t \log t)$ time:
$$ES^+ = \max_{j=1..t}\!\left(\frac{j}{t} - \frac{V(j) - j}{N - t}\right), \qquad ES^- = \min_{j=1..t}\!\left(\frac{j - 1}{t} - \frac{V(j) - j}{N - t}\right)$$
choosing $ES^+$ if $ES^+ > |ES^-|$, $ES^-$ if $|ES^-| > ES^+$, and $0$ on an exact tie (or flat query profile).

Given $ES_{\text{up}}$ and $ES_{\text{down}}$:
- When both sets match and have **opposite signs** ($ES_{\text{up}} \cdot ES_{\text{down}} < 0$), the connectivity score is $100 \cdot (ES_{\text{up}} - ES_{\text{down}}) / 2 \in [-100, +100]$.
- When both sets move in the same direction ($ES_{\text{up}} \cdot ES_{\text{down}} \ge 0$), the score cancels to $0.0$ per Lamb et al. (2006).
- Single-sided signatures fall back cleanly to $100 \cdot ES_{\text{up}}$ or $-100 \cdot ES_{\text{down}}$.

Rank-based scoring guarantees strict invariance to extreme fold-change outliers and `NaN` rows.

---

## 5. Multi-Hop CNA & Partial Fisher-Yates Permutation Sampling (`cna.rs`)

Multi-Hop Causal Network Analysis finds shortest directed paths of depth $K \in \{1, 2, 3\}$ through signaling edges ($E_a$) terminating in a transcriptional edge ($E_t$), collapses consistent virtual edges $(r, v)$, prunes nested depth hypotheses (`overlap_cnt > prev_overlap_cnt` and `p_val < prev_best_p * 0.95`), and computes degree-corrected network-bias p-values stratified by $\lfloor \log_2(\text{in\_degree}_t(g)) \rfloor$.

Within each in-degree bin pool of size $P$ requiring $d < P$ draws without replacement, **Causal Atlas** performs **partial Fisher-Yates** sampling in-place over a reusable per-bin vector using a deterministic `SplitMix64` PRNG:
```text
for i in 0..d:
    j = i + rng.gen_index(P - i)
    work_pool.swap(i, j)
```
Bins with $d \ge P$ are pre-accumulated outside the permutation loop. This eliminates `BTreeSet` heap allocations and coupon-collector slowdown while guaranteeing exact uniform $\binom{P}{d}^{-1}$ sampling across all `cna_permutations` iterations.

---

## 6. Molecule Activity Predictor (MAP) & SVG Geometry Clipping

- **MAP Propagation (`canonical.rs::simulate_map`)**: Propagates observed expression anchors `(fc / 3.0).tanh()` and pinned user/drug clamps `[-1.0, +1.0]` along signed pathway edges over `map_iterations` steps with `map_damping`.
- **SVG Node-Box Edge Clipping (`ca-app/src/svg_export.rs`)**: Each node renders as a `74x30` rounded `<rect>` centered at `(cx, cy)` (`[cx ± 37, cy ± 15]`). Prior to emitting `<path d="M sx sy Q mx my ex ey">`, both endpoints are clipped via ray-box intersection along the Quadratic Bezier endpoint tangents `(mx - cx, my - cy)` against `[cx ± (37 + 1.2), cy ± (15 + 1.2)]` so `#arrow-act` and `#arrow-inh` markers render cleanly outside node borders, accompanied by a `<linearGradient>` MAP activity colorbar legend.
