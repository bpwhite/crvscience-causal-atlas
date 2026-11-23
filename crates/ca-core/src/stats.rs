// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

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
    benjamini_hochberg_with_total(p_values, p_values.len())
}

/// Benjamini-Hochberg step-up False Discovery Rate adjustment (`q-values`) across
/// `total_hypotheses` (`m >= p_values.len()`).
///
/// When non-overlapping hypotheses with `p = 1.0` are omitted prior to constructing
/// `p_values`, passing the full tested hypothesis count `m` preserves exact multiple-testing
/// correction without materializing trailing `1.0` entries.
pub fn benjamini_hochberg_with_total(p_values: &[f64], total_hypotheses: usize) -> Vec<f64> {
    let k = p_values.len();
    if k == 0 {
        return Vec::new();
    }
    let m = total_hypotheses.max(k);
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .copied()
        .enumerate()
        .map(|(i, p)| (i, if p.is_finite() { p.clamp(0.0, 1.0) } else { 1.0 }))
        .collect();

    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut q_sorted = vec![1.0; k];
    let mut running_min = 1.0f64;
    let m_f = m as f64;

    for rank_rev in (0..k).rev() {
        let rank = (rank_rev + 1) as f64;
        let p = indexed[rank_rev].1;
        let raw_q = (p * m_f / rank).min(1.0);
        running_min = running_min.min(raw_q);
        q_sorted[rank_rev] = running_min;
    }

    let mut out = vec![1.0; k];
    for (rank_idx, &(orig_idx, _)) in indexed.iter().enumerate() {
        out[orig_idx] = q_sorted[rank_idx];
    }
    out
}

/// Evaluates `erfc(z)` for `z > 0` via the incomplete-gamma continued fraction
/// `Gamma(1/2, z^2) / sqrt(pi)` using modified Lentz's method.
fn erfc_continued_fraction(z: f64) -> f64 {
    const TINY: f64 = 1e-300;
    const EPS: f64 = 1e-15;
    let x = z * z;
    let mut b = x + 0.5;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b.max(TINY);
    let mut h = d;

    for k in 1..=128 {
        let kf = k as f64;
        let a = -kf * (kf - 0.5);
        b += 2.0;
        d = a * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + a / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() <= EPS {
            break;
        }
    }

    let inv_sqrt_pi = std::f64::consts::FRAC_2_SQRT_PI * 0.5;
    (h * z * inv_sqrt_pi * (-x).exp()).clamp(0.0, 1.0)
}

/// High-accuracy complementary error function `erfc(x)`:
/// - Uses Abramowitz & Stegun 7.1.26 rational polynomial for `|x| <= 5`.
/// - Uses incomplete-gamma continued fraction tail (`erfc_continued_fraction`) for `|x| > 5`
///   to retain relative precision in extreme p-value tails.
pub fn erfc_approx(x: f64) -> f64 {
    if !x.is_finite() {
        return if x.is_nan() {
            f64::NAN
        } else if x > 0.0 {
            0.0
        } else {
            2.0
        };
    }
    let z = x.abs();
    let r = if z > 5.0 {
        erfc_continued_fraction(z)
    } else {
        let t = 1.0 / (1.0 + 0.3275911 * z);
        let poly = t
            * (0.254829592
                + t * (-0.284496736
                    + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
        poly * (-z * z).exp()
    };
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
    let raw_mu = if regulator_mean_sign.is_finite() && dataset_mean_sign.is_finite() {
        regulator_mean_sign * dataset_mean_sign
    } else {
        0.0
    };
    let mu = raw_mu.clamp(-0.99, 0.99);
    let skew_flagged = mu.abs() > 0.25;

    for &(w, s_edge, s_data) in observations {
        if s_edge == 0 || s_data == 0 || !w.is_finite() || w <= 0.0 {
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

    #[test]
    fn kramer_z_skips_non_finite_weights_and_mean_signs() {
        let obs = vec![
            (f64::NAN, 1, 1),
            (f64::INFINITY, 1, 1),
            (f64::NEG_INFINITY, 1, 1),
            (1.0, 1, 1),
            (1.0, -1, -1),
            (1.0, 1, 1),
            (1.0, 1, 1),
        ];
        let (z, z_corr, _) = kramer_activation_z(&obs, 0.0, 0.0);
        assert!(z.is_finite() && (z - 2.0).abs() < 1e-9, "got z={z}");
        assert!(z_corr.is_finite());

        let (z_only_bad, z_corr_bad, _) =
            kramer_activation_z(&[(f64::NAN, 1, 1), (f64::INFINITY, 1, -1)], f64::NAN, 0.0);
        assert_eq!(z_only_bad, 0.0);
        assert_eq!(z_corr_bad, 0.0);
    }

    #[test]
    fn bh_fdr_with_total_accounts_for_filtered_null_hypotheses() {
        let subset = [0.001, 0.01];
        let mut padded = vec![0.001, 0.01];
        padded.extend(std::iter::repeat_n(1.0, 8));

        let q_total = benjamini_hochberg_with_total(&subset, 10);
        let q_padded = benjamini_hochberg(&padded);
        assert!((q_total[0] - q_padded[0]).abs() < 1e-12);
        assert!((q_total[1] - q_padded[1]).abs() < 1e-12);
        assert!((q_total[0] - 0.01).abs() < 1e-12);
        assert!((q_total[1] - 0.05).abs() < 1e-12);
    }

    #[test]
    fn erfc_approx_extreme_tail_continued_fraction_precision() {
        // Crossover behavior around |x| = 5.0 (true derivative |d/dx ln erfc(x)| ~ 2x = 10 at x=5)
        let at_5 = erfc_approx(5.0);
        let left = erfc_approx(4.9999);
        let right = erfc_approx(5.0001);
        assert!(left > at_5 && at_5 > right);
        assert!((left - right).abs() / left < 0.05);

        // Known reference values in extreme tail (|x| > 5):
        // erfc(5.0001) = 1.535893491062367e-12
        // erfc(6.0)    = 2.1519736712498913e-17
        // erfc(10.0)   = 2.088487583762545e-45
        // erfc(20.0)   = 5.395865611607901e-176
        assert!(
            (right - 1.535893491062367e-12).abs() / 1.535893491062367e-12 < 1e-12,
            "right={right:e}"
        );
        let e6 = erfc_approx(6.0);
        assert!((e6 - 2.1519736712498913e-17).abs() / 2.1519736712498913e-17 < 1e-12, "e6={e6:e}");
        let e10 = erfc_approx(10.0);
        assert!((e10 - 2.088487583762545e-45).abs() / 2.088487583762545e-45 < 1e-12, "e10={e10:e}");
        let e20 = erfc_approx(20.0);
        assert!(e20 > 0.0 && (e20 - 5.3958656116079005e-176).abs() / 5.3958656116079005e-176 < 1e-11, "e20={e20:e}");

        // Negative extreme tail approaches 2.0 cleanly and z_to_two_sided_p stays positive & accurate
        assert!((erfc_approx(-6.0) - (2.0 - e6)).abs() < 1e-15);
        let p10 = z_to_two_sided_p(10.0);
        assert!((p10 - 1.523970604832105e-23).abs() / 1.523970604832105e-23 < 1e-12, "p10={p10:e}");
    }
}
