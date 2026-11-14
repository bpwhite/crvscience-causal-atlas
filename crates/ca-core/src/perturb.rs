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
}

/// Evaluates directional connectivity between the query expression profile and open-source
/// chemical/genetic perturbation signatures (SigCom LINCS / L1000 consensus up & down sets).
pub fn run_perturbation_match(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
) -> Vec<PerturbationHit> {
    let mut hits = Vec::new();
    if active.all_ranked.is_empty() {
        return hits;
    }

    let gene_fc: BTreeMap<&str, f64> = active
        .all_ranked
        .iter()
        .map(|(g, fc)| (g.as_str(), *fc))
        .collect();

    let max_abs_fc = active
        .all_ranked
        .iter()
        .map(|(_, fc)| fc.abs())
        .fold(1.0f64, f64::max);

    for sig in &kg.signatures {
        let mut up_score = 0.0f64;
        let mut down_score = 0.0f64;
        let mut matched_up = Vec::new();
        let mut matched_down = Vec::new();

        for g in &sig.up_genes {
            if let Some(&fc) = gene_fc.get(g.as_str()) {
                up_score += fc / max_abs_fc;
                if active.degs.contains_key(g) {
                    matched_up.push(g.clone());
                }
            }
        }
        for g in &sig.down_genes {
            if let Some(&fc) = gene_fc.get(g.as_str()) {
                down_score -= fc / max_abs_fc;
                if active.degs.contains_key(g) {
                    matched_down.push(g.clone());
                }
            }
        }

        let total_overlap = matched_up.len() + matched_down.len();
        if total_overlap == 0 {
            continue;
        }

        let denom = (sig.up_genes.len() + sig.down_genes.len()).max(1) as f64;
        // Scale directional concordance into [-100, +100] CMap-like score.
        let raw_conn = ((up_score + down_score) / denom * 135.0).clamp(-100.0, 100.0);
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
        });
    }

    hits.sort_by(|a, b| {
        b.connectivity_score
            .abs()
            .partial_cmp(&a.connectivity_score.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
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
            let last_vec = &ordered.last().unwrap().1;
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
