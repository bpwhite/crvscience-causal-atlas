#![forbid(unsafe_code)]

use std::sync::OnceLock;

/// Precomputed `ln(k!)` table up to 100,000 with Stirling expansion fallback above 100,000.
fn ln_factorial(n: usize) -> f64 {
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    let tbl = TABLE.get_or_init(|| {
        let max_n = 100_000usize;
        let mut v = Vec::with_capacity(max_n + 1);
        v.push(0.0);
        let mut acc = 0.0;
        for i in 1..=max_n {
            acc += (i as f64).ln();
            v.push(acc);
        }
        v
    });
    if n < tbl.len() {
        tbl[n]
    } else {
        let x = n as f64;
        x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln() + 1.0 / (12.0 * x)
    }
}

fn ln_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

/// Exact one-sided right-tail Fisher's Exact Test (hypergeometric enrichment p-value):
/// `P(X >= overlap)` given `universe` (`N`), `target_set_size` (`M`), and `dataset_size` (`n`).
pub fn fishers_exact_right_tail(
    overlap: usize,
    target_set_size: usize,
    dataset_size: usize,
    universe_size: usize,
) -> f64 {
    if overlap == 0 {
        return 1.0;
    }
    let n_univ = universe_size.max(target_set_size.max(dataset_size));
    if n_univ == 0 || target_set_size == 0 || dataset_size == 0 {
        return 1.0;
    }
    let max_k = target_set_size.min(dataset_size);
    if overlap > max_k {
        return 0.0;
    }
    let min_k = (target_set_size + dataset_size).saturating_sub(n_univ);
    let start_k = overlap.max(min_k);

    let denom = ln_choose(n_univ, dataset_size);
    let mut log_terms = Vec::with_capacity(max_k - start_k + 1);
    let mut max_log = f64::NEG_INFINITY;

    for k in start_k..=max_k {
        let lp = ln_choose(target_set_size, k)
            + ln_choose(n_univ - target_set_size, dataset_size - k)
            - denom;
        if lp > max_log {
            max_log = lp;
        }
        log_terms.push(lp);
    }

    if !max_log.is_finite() {
        return 1.0;
    }

    let sum_scaled: f64 = log_terms.into_iter().map(|lp| (lp - max_log).exp()).sum();
    (max_log + sum_scaled.ln()).exp().clamp(0.0, 1.0)
}

/// Benjamini-Hochberg step-up False Discovery Rate adjustment (`q-values`).
pub fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let m = p_values.len();
    if m == 0 {
        return Vec::new();
    }
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .copied()
        .enumerate()
        .map(|(i, p)| (i, if p.is_finite() { p.clamp(0.0, 1.0) } else { 1.0 }))
        .collect();

    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut q_sorted = vec![1.0; m];
    let mut running_min = 1.0f64;
    let m_f = m as f64;

    for rank_rev in (0..m).rev() {
        let rank = (rank_rev + 1) as f64;
        let p = indexed[rank_rev].1;
        let raw_q = (p * m_f / rank).min(1.0);
        running_min = running_min.min(raw_q);
        q_sorted[rank_rev] = running_min;
    }

    let mut out = vec![1.0; m];
    for (rank_idx, &(orig_idx, _)) in indexed.iter().enumerate() {
        out[orig_idx] = q_sorted[rank_idx];
    }
    out
}

/// High-accuracy complementary error function `erfc(x)` (Abramowitz & Stegun 7.1.26).
pub fn erfc_approx(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * z);
    let poly = t
        * (0.254829592
            + t * (-0.284496736
                + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let r = poly * (-z * z).exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

/// Two-sided standard normal p-value from a Z-score: `2 * (1 - Phi(|z|))`.
pub fn z_to_two_sided_p(z: f64) -> f64 {
    if !z.is_finite() {
        return 1.0;
    }
    erfc_approx(z.abs() / std::f64::consts::SQRT_2).clamp(0.0, 1.0)
}

/// Krämer et al. (2014) weighted signed Activation Z-score and skewness-adjusted Z-score.
///
/// - `observations`: slice of `(weight_i, regulator_edge_sign_i, dataset_gene_sign_i)` where signs are `+1` or `-1`.
/// - `regulator_mean_sign` (`mu_1`): mean sign across ALL signed targets of regulator `r` in the knowledge graph.
/// - `dataset_mean_sign` (`mu_2`): mean sign across ALL active genes in dataset `D`.
///
/// Returns `(z_score, bias_corrected_z, is_skew_flagged)`.
pub fn kramer_activation_z(
    observations: &[(f64, i8, i8)],
    regulator_mean_sign: f64,
    dataset_mean_sign: f64,
) -> (f64, f64, bool) {
    if observations.is_empty() {
        return (0.0, 0.0, false);
    }
    let mut num = 0.0f64;
    let mut sum_w2 = 0.0f64;
    let mut num_corr = 0.0f64;
    let mu = (regulator_mean_sign * dataset_mean_sign).clamp(-0.99, 0.99);
    let skew_flagged = mu.abs() > 0.25;

    for &(w, s_edge, s_data) in observations {
        if s_edge == 0 || s_data == 0 || w <= 0.0 {
            continue;
        }
        let x_i = f64::from(s_edge) * f64::from(s_data);
        num += w * x_i;
        num_corr += w * (x_i - mu);
        sum_w2 += w * w;
    }

    if sum_w2 <= 1e-12 {
        return (0.0, 0.0, skew_flagged);
    }

    let z = num / sum_w2.sqrt();
    let denom_corr = ((1.0 - mu * mu).max(0.05) * sum_w2).sqrt();
    let z_corr = num_corr / denom_corr;
    (z, z_corr, skew_flagged)
}

/// Deterministic fast SplitMix64 PRNG so degree-preserving Monte Carlo permutations are reproducible
/// without requiring external C/OS entropy crates.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn gen_index(&mut self, upper_exclusive: usize) -> usize {
        if upper_exclusive <= 1 {
            0
        } else {
            (self.next_u64() % (upper_exclusive as u64)) as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fishers_exact_matches_known_2x2_contingency() {
        // Universe = 20, Pathway M = 6, Dataset n = 5, Overlap k = 4
        // Exact right-tail = [C(6,4)*C(14,1) + C(6,5)*C(14,0)] / C(20,5) = (15*14 + 6*1) / 15504 = 216 / 15504
        let p = fishers_exact_right_tail(4, 6, 5, 20);
        let expected = 216.0 / 15504.0;
        assert!((p - expected).abs() < 1e-10, "got {p}, expected {expected}");
    }

    #[test]
    fn bh_fdr_preserves_monotonicity() {
        let raw = vec![0.001, 0.04, 0.01, 0.20];
        let q = benjamini_hochberg(&raw);
        assert!(q[0] < q[2]);
        assert!(q[2] < q[1]);
        assert!(q[1] < q[3]);
    }

    #[test]
    fn kramer_z_detects_activation_and_inhibition() {
        let obs_act = vec![(1.0, 1, 1), (1.0, -1, -1), (1.0, 1, 1), (1.0, 1, 1)];
        let (z_act, _, _) = kramer_activation_z(&obs_act, 0.0, 0.0);
        assert!((z_act - 2.0).abs() < 1e-9);

        let obs_inh = vec![(1.0, 1, -1), (1.0, -1, 1), (1.0, 1, -1), (1.0, 1, -1)];
        let (z_inh, _, _) = kramer_activation_z(&obs_inh, 0.0, 0.0);
        assert!((z_inh + 2.0).abs() < 1e-9);
    }
}
