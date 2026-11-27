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

use crate::dataset::ActiveDataset;
use crate::graph::{KnowledgeGraph, PerturbationKind};
use crate::CoreAnalysisReport;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerturbationHit {
    pub signature_id: String,
    pub perturbagen: String,
    pub kind: PerturbationKind,
    pub cell_context: String,
    pub mechanism_of_action: String,
    pub source_db: String,
    pub source_url: String,
    /// Connectivity score in `[-100.0, +100.0]`:
    /// Negative (`<= -40`) indicates a therapeutic **Reverser**; positive (`>= +40`) indicates a **Mimicker**.
    pub connectivity_score: f64,
    pub mode: String,
    pub matched_up: Vec<String>,
    pub matched_down: Vec<String>,
    #[serde(default)]
    pub target_clamps: Vec<(String, i8)>,
}

/// Two-sample unweighted Kolmogorov-Smirnov / Subramanian GSEA running-sum enrichment score
/// over `n` ranked genes (`1..=n`), returning a value in `[-1.0, +1.0]`.
fn ks_rank_enrichment(sorted_ranks: &[usize], n: usize) -> f64 {
    let t = sorted_ranks.len();
    if t == 0 || n == 0 || t >= n {
        return 0.0;
    }
    let t_f = t as f64;
    let miss_denom = (n - t) as f64;
    let mut max_pos = 0.0f64;
    let mut min_neg = 0.0f64;

    for (idx, &v) in sorted_ranks.iter().enumerate() {
        let j = (idx + 1) as f64;
        let misses_before = (v.saturating_sub(idx + 1)) as f64;
        let run_before = ((j - 1.0) / t_f) - (misses_before / miss_denom);
        let run_after = (j / t_f) - (misses_before / miss_denom);
        if run_after > max_pos {
            max_pos = run_after;
        }
        if run_before < min_neg {
            min_neg = run_before;
        }
    }

    if (max_pos - min_neg.abs()).abs() <= 1e-12 {
        0.0
    } else if max_pos > min_neg.abs() {
        max_pos
    } else {
        min_neg
    }
}

/// Evaluates rank-based signed CMap / L1000 connectivity (two-sample KS enrichment over the
/// finite descending `log2FC` profile) between the query profile and open-source perturbation
/// signatures (SigCom LINCS / L1000 consensus up & down sets).
pub fn run_perturbation_match(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
) -> Vec<PerturbationHit> {
    let mut hits = Vec::new();
    if active.all_ranked.is_empty() {
        return hits;
    }

    // Guard NaN/Inf fold-changes and deduplicate symbols keeping strongest magnitude.
    let mut finite_map: BTreeMap<String, f64> = BTreeMap::new();
    for (g, fc) in &active.all_ranked {
        let sym = g.trim().to_uppercase();
        if sym.is_empty() || !fc.is_finite() {
            continue;
        }
        finite_map
            .entry(sym)
            .and_modify(|prev| {
                if fc.abs() > prev.abs() {
                    *prev = *fc;
                }
            })
            .or_insert(*fc);
    }

    if finite_map.is_empty() {
        return hits;
    }

    let mut ranked: Vec<(String, f64)> = finite_map.into_iter().collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let n_ranked = ranked.len();
    let flat_profile = match (ranked.first(), ranked.last()) {
        (Some(first), Some(last)) => (first.1 - last.1).abs() <= 1e-12,
        _ => true,
    };

    let gene_rank: BTreeMap<String, usize> = ranked
        .into_iter()
        .enumerate()
        .map(|(idx, (g, _))| (g, idx + 1))
        .collect();

    for sig in &kg.signatures {
        let mut up_ranks = Vec::new();
        let mut down_ranks = Vec::new();
        let mut matched_up = Vec::new();
        let mut matched_down = Vec::new();

        for g in &sig.up_genes {
            if let Some(&r) = gene_rank.get(g) {
                up_ranks.push(r);
                if active
                    .degs
                    .get(g)
                    .is_some_and(|m| m.log2_fc.is_finite())
                {
                    matched_up.push(g.clone());
                }
            }
        }
        for g in &sig.down_genes {
            if let Some(&r) = gene_rank.get(g) {
                down_ranks.push(r);
                if active
                    .degs
                    .get(g)
                    .is_some_and(|m| m.log2_fc.is_finite())
                {
                    matched_down.push(g.clone());
                }
            }
        }

        let total_overlap = matched_up.len() + matched_down.len();
        if total_overlap == 0 {
            continue;
        }

        up_ranks.sort_unstable();
        up_ranks.dedup();
        down_ranks.sort_unstable();
        down_ranks.dedup();

        let raw_score = if flat_profile {
            0.0
        } else {
            let es_up = ks_rank_enrichment(&up_ranks, n_ranked);
            let es_down = ks_rank_enrichment(&down_ranks, n_ranked);
            if !up_ranks.is_empty() && !down_ranks.is_empty() {
                if es_up * es_down < 0.0 {
                    (es_up - es_down) * 0.5
                } else {
                    0.0
                }
            } else if !up_ranks.is_empty() && sig.down_genes.is_empty() {
                es_up
            } else if !down_ranks.is_empty() && sig.up_genes.is_empty() {
                -es_down
            } else {
                (es_up - es_down) * 0.5
            }
        };

        // Scale signed KS connectivity into [-100, +100].
        let raw_conn = (raw_score * 100.0).clamp(-100.0, 100.0);
        let mode = if raw_conn <= -25.0 {
            "Reverser (Candidate Antagonist)".to_string()
        } else if raw_conn >= 25.0 {
            "Mimicker (Concordant Agonist)".to_string()
        } else {
            "Weak / Orthogonal".to_string()
        };

        hits.push(PerturbationHit {
            signature_id: sig.id.clone(),
            perturbagen: sig.perturbagen.clone(),
            kind: sig.kind,
            cell_context: sig.cell_context.clone(),
            mechanism_of_action: sig.mechanism_of_action.clone(),
            source_db: sig.source_db.clone(),
            source_url: sig.source_url.clone(),
            connectivity_score: raw_conn,
            mode,
            matched_up,
            matched_down,
            target_clamps: sig.target_clamps.clone(),
        });
    }

    hits.sort_by(|a, b| {
        b.connectivity_score
            .abs()
            .total_cmp(&a.connectivity_score.abs())
            .then_with(|| {
                (b.matched_up.len() + b.matched_down.len())
                    .cmp(&(a.matched_up.len() + a.matched_down.len()))
            })
            .then_with(|| a.signature_id.cmp(&b.signature_id))
    });

    hits
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ComparisonMetric {
    CanonicalPathwaysZ,
    UpstreamRegulatorsZ,
    DownstreamEffectsZ,
}

impl ComparisonMetric {
    pub fn label(self) -> &'static str {
        match self {
            Self::CanonicalPathwaysZ => "Canonical Pathways (Activation Z)",
            Self::UpstreamRegulatorsZ => "Upstream Regulators (Activation Z)",
            Self::DownstreamEffectsZ => "Downstream Effects (Activation Z)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonMatrix {
    pub metric: ComparisonMetric,
    pub cohort_names: Vec<String>,
    pub entity_labels: Vec<String>,
    /// Row-major `[entity_idx][cohort_idx]` values (`Activation Z-score`).
    pub values: Vec<Vec<f64>>,
    /// Row-major `[entity_idx][cohort_idx]` `-log10(FDR)` significance weights.
    pub significance: Vec<Vec<f64>>,
}

/// Aligns multiple `CoreAnalysisReport` cohorts side-by-side and orders rows by total absolute Z magnitude
/// and nearest-neighbor Euclidean similarity for visual heatmap clustering.
pub fn run_comparison_analysis(
    reports: &[CoreAnalysisReport],
    metric: ComparisonMetric,
    max_rows: usize,
) -> ComparisonMatrix {
    let cohort_names: Vec<String> = reports.iter().map(|r| r.dataset_name.clone()).collect();
    let n_cohorts = reports.len();

    // Map `entity_label -> Vec<(z, neg_log10_q)>`
    let mut table: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    let mut all_keys = BTreeSet::new();

    for (c_idx, rep) in reports.iter().enumerate() {
        match metric {
            ComparisonMetric::CanonicalPathwaysZ => {
                for hit in &rep.canonical_pathways {
                    let key = hit.pathway_name.clone();
                    all_keys.insert(key.clone());
                    let row = table
                        .entry(key)
                        .or_insert_with(|| vec![(0.0, 0.0); n_cohorts]);
                    row[c_idx] = (hit.activation_z.unwrap_or(0.0), hit.neg_log10_fdr);
                }
            }
            ComparisonMetric::UpstreamRegulatorsZ => {
                for hit in &rep.upstream_regulators {
                    let key = format!("{} ({})", hit.regulator_id, hit.kind.label());
                    all_keys.insert(key.clone());
                    let row = table
                        .entry(key)
                        .or_insert_with(|| vec![(0.0, 0.0); n_cohorts]);
                    row[c_idx] = (hit.activation_z, -hit.fdr.max(1e-300).log10());
                }
            }
            ComparisonMetric::DownstreamEffectsZ => {
                for hit in &rep.downstream_effects {
                    let key = hit.term_name.clone();
                    all_keys.insert(key.clone());
                    let row = table
                        .entry(key)
                        .or_insert_with(|| vec![(0.0, 0.0); n_cohorts]);
                    row[c_idx] = (hit.activation_z, hit.neg_log10_p);
                }
            }
        }
    }

    let mut ranked_rows: Vec<(String, Vec<(f64, f64)>, f64)> = table
        .into_iter()
        .map(|(k, v)| {
            let energy: f64 = v.iter().map(|(z, s)| z.abs() + 0.1 * s).sum();
            (k, v, energy)
        })
        .collect();

    ranked_rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    ranked_rows.truncate(max_rows.max(5));

    // Greedy nearest-neighbor leaf ordering so related patterns cluster together in the heatmap.
    let mut ordered = Vec::with_capacity(ranked_rows.len());
    if !ranked_rows.is_empty() {
        ordered.push(ranked_rows.remove(0));
        while !ranked_rows.is_empty() {
            let Some(last_vec) = ordered.last().map(|r| &r.1) else {
                break;
            };
            let mut best_i = 0usize;
            let mut best_dist = f64::INFINITY;
            for (i, cand) in ranked_rows.iter().enumerate() {
                let dist: f64 = last_vec
                    .iter()
                    .zip(&cand.1)
                    .map(|((z1, _), (z2, _))| (z1 - z2).powi(2))
                    .sum();
                if dist < best_dist {
                    best_dist = dist;
                    best_i = i;
                }
            }
            ordered.push(ranked_rows.remove(best_i));
        }
    }

    let mut entity_labels = Vec::with_capacity(ordered.len());
    let mut values = Vec::with_capacity(ordered.len());
    let mut significance = Vec::with_capacity(ordered.len());

    for (label, cells, _) in ordered {
        entity_labels.push(label);
        values.push(cells.iter().map(|(z, _)| *z).collect());
        significance.push(cells.iter().map(|(_, s)| *s).collect());
    }

    ComparisonMatrix {
        metric,
        cohort_names,
        entity_labels,
        values,
        significance,
    }
}

/// Per-regulator rank & activation Z-score diff between two perturbation/cohort reports.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegulatorRankDiffEntry {
    pub regulator_id: String,
    pub kind: String,
    pub rank_a: Option<usize>,
    pub rank_b: Option<usize>,
    pub rank_shift: Option<i64>,
    pub z_a: Option<f64>,
    pub z_b: Option<f64>,
    pub delta_z: f64,
    pub fdr_a: Option<f64>,
    pub fdr_b: Option<f64>,
    pub direction_flip: bool,
    pub status: String,
}

/// Multi-perturbation regulator rank-comparison diff digest comparing condition A vs. B.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerturbationRankDiffDigest {
    pub schema_version: String,
    pub cohort_a_name: String,
    pub cohort_b_name: String,
    pub regulators_a: usize,
    pub regulators_b: usize,
    pub shared_count: usize,
    pub direction_flip_count: usize,
    pub spearman_rank_rho: f64,
    pub entries: Vec<RegulatorRankDiffEntry>,
}

/// Computes a ranked differential regulator comparison digest between two `CoreAnalysisReport` runs.
pub fn compute_regulator_rank_diff(
    report_a: &CoreAnalysisReport,
    report_b: &CoreAnalysisReport,
) -> PerturbationRankDiffDigest {
    let mut rank_a_list: Vec<&crate::ura::UpstreamRegulatorHit> =
        report_a.upstream_regulators.iter().collect();
    rank_a_list.sort_by(|x, y| {
        y.activation_z
            .abs()
            .total_cmp(&x.activation_z.abs())
            .then_with(|| x.regulator_id.cmp(&y.regulator_id))
    });
    let mut rank_b_list: Vec<&crate::ura::UpstreamRegulatorHit> =
        report_b.upstream_regulators.iter().collect();
    rank_b_list.sort_by(|x, y| {
        y.activation_z
            .abs()
            .total_cmp(&x.activation_z.abs())
            .then_with(|| x.regulator_id.cmp(&y.regulator_id))
    });

    let mut map_a: BTreeMap<String, (usize, &crate::ura::UpstreamRegulatorHit)> = BTreeMap::new();
    for (idx, hit) in rank_a_list.iter().enumerate() {
        map_a.entry(hit.regulator_id.clone()).or_insert((idx + 1, *hit));
    }
    let mut map_b: BTreeMap<String, (usize, &crate::ura::UpstreamRegulatorHit)> = BTreeMap::new();
    for (idx, hit) in rank_b_list.iter().enumerate() {
        map_b.entry(hit.regulator_id.clone()).or_insert((idx + 1, *hit));
    }

    let mut all_ids: BTreeSet<String> = map_a.keys().cloned().collect();
    all_ids.extend(map_b.keys().cloned());

    let mut entries = Vec::with_capacity(all_ids.len());
    let mut shared_count = 0usize;
    let mut direction_flip_count = 0usize;
    let mut shared_d2_sum = 0.0_f64;

    for reg_id in all_ids {
        let a_opt = map_a.get(&reg_id).copied();
        let b_opt = map_b.get(&reg_id).copied();
        let kind = a_opt
            .map(|(_, h)| h.kind.label())
            .or_else(|| b_opt.map(|(_, h)| h.kind.label()))
            .unwrap_or("Regulator")
            .to_string();
        let rank_a = a_opt.map(|(r, _)| r);
        let rank_b = b_opt.map(|(r, _)| r);
        let rank_shift = match (rank_a, rank_b) {
            (Some(ra), Some(rb)) => Some((ra as i64) - (rb as i64)),
            _ => None,
        };
        let z_a = a_opt.map(|(_, h)| h.activation_z);
        let z_b = b_opt.map(|(_, h)| h.activation_z);
        let fdr_a = a_opt.map(|(_, h)| h.fdr);
        let fdr_b = b_opt.map(|(_, h)| h.fdr);
        let delta_z = z_b.unwrap_or(0.0) - z_a.unwrap_or(0.0);

        let direction_flip = match (z_a, z_b) {
            (Some(za), Some(zb)) => za.abs() >= 2.0 && zb.abs() >= 2.0 && za * zb < 0.0,
            _ => false,
        };

        let status = match (a_opt, b_opt) {
            (Some(_), Some(_)) => {
                shared_count += 1;
                if let (Some(ra), Some(rb)) = (rank_a, rank_b) {
                    let d = (ra as f64) - (rb as f64);
                    shared_d2_sum += d * d;
                }
                if direction_flip {
                    direction_flip_count += 1;
                    String::from("Reversed")
                } else {
                    String::from("Shared")
                }
            }
            (None, Some(_)) => String::from("Gained_in_B"),
            (Some(_), None) => String::from("Lost_in_B"),
            (None, None) => String::from("Unobserved"),
        };

        entries.push(RegulatorRankDiffEntry {
            regulator_id: reg_id,
            kind,
            rank_a,
            rank_b,
            rank_shift,
            z_a,
            z_b,
            delta_z,
            fdr_a,
            fdr_b,
            direction_flip,
            status,
        });
    }

    entries.sort_by(|x, y| {
        y.delta_z
            .abs()
            .total_cmp(&x.delta_z.abs())
            .then_with(|| x.regulator_id.cmp(&y.regulator_id))
    });

    let spearman_rank_rho = if shared_count >= 2 {
        let n = shared_count as f64;
        (1.0 - (6.0 * shared_d2_sum) / (n * (n * n - 1.0))).clamp(-1.0, 1.0)
    } else {
        0.0
    };

    PerturbationRankDiffDigest {
        schema_version: String::from("causal-atlas-rank-diff/v1"),
        cohort_a_name: report_a.dataset_name.clone(),
        cohort_b_name: report_b.dataset_name.clone(),
        regulators_a: map_a.len(),
        regulators_b: map_b.len(),
        shared_count,
        direction_flip_count,
        spearman_rank_rho,
        entries,
    }
}

fn json_escape_ca(s: &str) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(s.len() + 4);
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

/// Formats [`PerturbationRankDiffDigest`] as pretty-printed JSON.
pub fn format_regulator_rank_diff_json(digest: &PerturbationRankDiffDigest) -> Result<String, String> {
    use std::fmt::Write as _;
    let fmt_opt_usize = |v: Option<usize>| v.map(|x| x.to_string()).unwrap_or_else(|| String::from("null"));
    let fmt_opt_i64 = |v: Option<i64>| v.map(|x| x.to_string()).unwrap_or_else(|| String::from("null"));
    let fmt_opt_f64 = |v: Option<f64>| {
        v.filter(|x| x.is_finite())
            .map(|x| format!("{x:.6}"))
            .unwrap_or_else(|| String::from("null"))
    };
    let mut out = String::with_capacity(2048);
    let _ = writeln!(out, "{{");
    let _ = writeln!(
        out,
        "  \"schema_version\": \"{}\",",
        json_escape_ca(&digest.schema_version)
    );
    let _ = writeln!(
        out,
        "  \"cohort_a_name\": \"{}\",",
        json_escape_ca(&digest.cohort_a_name)
    );
    let _ = writeln!(
        out,
        "  \"cohort_b_name\": \"{}\",",
        json_escape_ca(&digest.cohort_b_name)
    );
    let _ = writeln!(out, "  \"regulators_a\": {},", digest.regulators_a);
    let _ = writeln!(out, "  \"regulators_b\": {},", digest.regulators_b);
    let _ = writeln!(out, "  \"shared_count\": {},", digest.shared_count);
    let _ = writeln!(
        out,
        "  \"direction_flip_count\": {},",
        digest.direction_flip_count
    );
    let rho = if digest.spearman_rank_rho.is_finite() {
        format!("{:.6}", digest.spearman_rank_rho)
    } else {
        String::from("0.0")
    };
    let _ = writeln!(out, "  \"spearman_rank_rho\": {rho},");
    let _ = writeln!(out, "  \"entries\": [");
    let n = digest.entries.len();
    for (idx, e) in digest.entries.iter().enumerate() {
        let comma = if idx + 1 < n { "," } else { "" };
        let dz = if e.delta_z.is_finite() {
            format!("{:.6}", e.delta_z)
        } else {
            String::from("0.0")
        };
        let _ = writeln!(
            out,
            "    {{\"regulator_id\": \"{}\", \"kind\": \"{}\", \"status\": \"{}\", \"rank_a\": {}, \"rank_b\": {}, \"rank_shift\": {}, \"z_a\": {}, \"z_b\": {}, \"delta_z\": {}, \"fdr_a\": {}, \"fdr_b\": {}, \"direction_flip\": {}}}{}",
            json_escape_ca(&e.regulator_id),
            json_escape_ca(&e.kind),
            json_escape_ca(&e.status),
            fmt_opt_usize(e.rank_a),
            fmt_opt_usize(e.rank_b),
            fmt_opt_i64(e.rank_shift),
            fmt_opt_f64(e.z_a),
            fmt_opt_f64(e.z_b),
            dz,
            fmt_opt_f64(e.fdr_a),
            fmt_opt_f64(e.fdr_b),
            e.direction_flip,
            comma
        );
    }
    let _ = writeln!(out, "  ]");
    let _ = writeln!(out, "}}");
    Ok(out)
}

/// Formats [`PerturbationRankDiffDigest`] as a tidy CSV ledger.
pub fn format_regulator_rank_diff_csv(digest: &PerturbationRankDiffDigest) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(1024);
    out.push_str(
        "regulator_id,kind,status,rank_a,rank_b,rank_shift,z_a,z_b,delta_z,fdr_a,fdr_b,direction_flip\n",
    );
    let fmt_opt_usize = |v: Option<usize>| v.map(|x| x.to_string()).unwrap_or_default();
    let fmt_opt_i64 = |v: Option<i64>| v.map(|x| x.to_string()).unwrap_or_default();
    let fmt_opt_f64 = |v: Option<f64>| v.map(|x| format!("{x:.6}")).unwrap_or_default();
    for e in &digest.entries {
        let _ = writeln!(
            out,
            "{},{},{},{},{},{},{},{},{:.6},{},{},{}",
            e.regulator_id,
            e.kind,
            e.status,
            fmt_opt_usize(e.rank_a),
            fmt_opt_usize(e.rank_b),
            fmt_opt_i64(e.rank_shift),
            fmt_opt_f64(e.z_a),
            fmt_opt_f64(e.z_b),
            e.delta_z,
            fmt_opt_f64(e.fdr_a),
            fmt_opt_f64(e.fdr_b),
            e.direction_flip
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::builtin_cohorts;
    use crate::seed_kb::build_seed_knowledge_graph;
    use crate::CoreAnalysisConfig;

    #[test]
    fn seed_kb_has_expanded_valid_signatures_and_scores_reversers_and_mimickers() {
        let kg = build_seed_knowledge_graph();
        assert!(
            kg.signatures.len() >= 25,
            "Expected >= 25 seed signatures, got {}",
            kg.signatures.len()
        );

        let mut seen_ids = BTreeSet::new();
        let mut seen_kinds = BTreeSet::new();
        for sig in &kg.signatures {
            assert!(
                seen_ids.insert(sig.id.clone()),
                "Duplicate signature ID: {}",
                sig.id
            );
            seen_kinds.insert(format!("{:?}", sig.kind));
            for g in sig
                .up_genes
                .iter()
                .chain(sig.down_genes.iter())
                .chain(sig.target_clamps.iter().map(|(g, _)| g))
            {
                assert!(
                    kg.nodes.contains_key(g),
                    "Signature {} references unknown KB symbol {g}",
                    sig.id
                );
            }
        }
        assert_eq!(
            seen_kinds.len(),
            4,
            "All 4 PerturbationKind variants should be represented"
        );

        let cfg = CoreAnalysisConfig::default();
        for ds in builtin_cohorts() {
            let active = ds.filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
            let hits = run_perturbation_match(&active, &kg);
            let reverser = hits.iter().find(|h| h.connectivity_score <= -25.0);
            let mimicker = hits.iter().find(|h| h.connectivity_score >= 25.0);
            assert!(
                reverser.is_some(),
                "Cohort {} should have at least one strong Reverser",
                ds.name
            );
            assert!(
                mimicker.is_some(),
                "Cohort {} should have at least one strong Mimicker",
                ds.name
            );
            let rev = reverser.unwrap();
            assert!(
                !rev.matched_up.is_empty() || !rev.matched_down.is_empty(),
                "Reverser {} should expose matched genes",
                rev.perturbagen
            );
        }
    }

    #[test]
    fn cmap_guards_nan_fold_changes_and_scores_rank_based_signed_enrichment() {
        use crate::dataset::{ActiveDataset, GeneMeasurement};
        use crate::graph::PerturbationSignature;

        let mut kg = KnowledgeGraph::new();
        for g in ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "NAN_GENE", "INF_GENE"] {
            kg.ensure_gene_node(g, "Test");
        }
        kg.signatures.push(PerturbationSignature {
            id: "REV_1".to_string(),
            perturbagen: "TestReverser".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "Test".to_string(),
            mechanism_of_action: "Reverses up/down".to_string(),
            source_db: "Test".to_string(),
            source_url: String::new(),
            up_genes: vec!["G7".into(), "G8".into()],
            down_genes: vec!["G1".into(), "G2".into()],
            target_clamps: Vec::new(),
        });
        kg.signatures.push(PerturbationSignature {
            id: "SAME_DIR".to_string(),
            perturbagen: "SameDirNoise".to_string(),
            kind: PerturbationKind::SmallMoleculeCompound,
            cell_context: "Test".to_string(),
            mechanism_of_action: "Same-direction should cancel to 0".to_string(),
            source_db: "Test".to_string(),
            source_url: String::new(),
            up_genes: vec!["G1".into()],
            down_genes: vec!["G2".into()],
            target_clamps: Vec::new(),
        });
        kg.rebuild_indices();

        let make_active = |outlier_fc: f64, include_bad: bool| {
            let mut degs = BTreeMap::new();
            let mut ranked = vec![
                ("G1".to_string(), outlier_fc),
                ("G2".to_string(), 3.0),
                ("G3".to_string(), 2.0),
                ("G4".to_string(), 1.0),
                ("G5".to_string(), -0.5),
                ("G6".to_string(), -1.5),
                ("G7".to_string(), -2.5),
                ("G8".to_string(), -3.5),
            ];
            if include_bad {
                ranked.push(("NAN_GENE".to_string(), f64::NAN));
                ranked.push(("INF_GENE".to_string(), f64::INFINITY));
            }
            for (s, fc) in &ranked {
                degs.insert(
                    s.clone(),
                    GeneMeasurement {
                        symbol: s.clone(),
                        log2_fc: *fc,
                        p_value: 1e-4,
                        fdr: 1e-3,
                    },
                );
            }
            ActiveDataset {
                name: "NanGuardTest".to_string(),
                total_measured: ranked.len(),
                mapped_to_graph: ranked.len(),
                degs,
                all_ranked: ranked,
            }
        };

        let hits_base = run_perturbation_match(&make_active(4.0, false), &kg);
        let hits_nan_outlier = run_perturbation_match(&make_active(400.0, true), &kg);

        let rev_base = hits_base.iter().find(|h| h.signature_id == "REV_1").unwrap();
        let rev_nan = hits_nan_outlier.iter().find(|h| h.signature_id == "REV_1").unwrap();
        assert!(rev_base.connectivity_score.is_finite());
        assert!((rev_base.connectivity_score + 100.0).abs() < 1e-9);
        assert!(
            (rev_base.connectivity_score - rev_nan.connectivity_score).abs() < 1e-9,
            "Rank-based CMap score must be invariant to extreme outliers and NaN/Inf rows"
        );

        let same_dir = hits_nan_outlier
            .iter()
            .find(|h| h.signature_id == "SAME_DIR")
            .unwrap();
        assert_eq!(same_dir.connectivity_score, 0.0);
    }

    #[test]
    fn regulator_rank_diff_digest_detects_flips_ranks_and_exports_json_csv() {
        let kg = build_seed_knowledge_graph();
        let cfg = CoreAnalysisConfig::default();
        let cohorts = builtin_cohorts();
        let rep_a = crate::run_core_analysis(&cohorts[0], &kg, &cfg);
        let rep_b = crate::run_core_analysis(&cohorts[1], &kg, &cfg);

        let digest = compute_regulator_rank_diff(&rep_a, &rep_b);
        assert_eq!(digest.schema_version, "causal-atlas-rank-diff/v1");
        assert!(!digest.entries.is_empty());
        assert!(digest.spearman_rank_rho.is_finite());

        let json = format_regulator_rank_diff_json(&digest).expect("serialize rank diff JSON");
        assert!(json.contains("\"schema_version\": \"causal-atlas-rank-diff/v1\""));
        assert!(json.contains("\"spearman_rank_rho\""));

        let csv = format_regulator_rank_diff_csv(&digest);
        assert!(csv.starts_with(
            "regulator_id,kind,status,rank_a,rank_b,rank_shift,z_a,z_b,delta_z,fdr_a,fdr_b,direction_flip\n"
        ));
        assert!(csv.lines().count() >= 2);
    }
}
