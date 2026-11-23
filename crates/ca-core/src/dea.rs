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
use crate::graph::{EdgeSign, EntityKind, KnowledgeGraph};
use crate::stats::{benjamini_hochberg_with_total, fishers_exact_right_tail, kramer_activation_z};
use crate::ura::{PredictedState, UpstreamRegulatorHit};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DownstreamEffectHit {
    pub term_id: String,
    pub term_name: String,
    pub category: String,
    pub kind: EntityKind,
    pub source_db: String,
    pub predicted_direction: String,
    pub activation_z: f64,
    pub p_value: f64,
    pub fdr: f64,
    pub neg_log10_p: f64,
    pub overlap_genes: Vec<String>,
}

/// Downstream Effects Analysis (DEA) per Krämer et al. (2014), Section 3.7:
/// Evaluates enrichment and directional increase/decrease of biological processes and diseases (`E_p`).
pub fn run_downstream_effects(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
) -> Vec<DownstreamEffectHit> {
    let universe = &kg.pathway_gene_universe;
    let n_univ = universe.len().max(1);
    let n_dataset = active.degs.len().max(1);
    let mu_2 = active.summary().mean_sign_skew;

    let mut hits = Vec::new();
    let mut p_vals = Vec::new();
    let mut total_tested_terms = 0usize;

    for (proc_id, edge_indices) in &kg.in_edges_p {
        let mut gene_map: BTreeMap<String, (EdgeSign, f64)> = BTreeMap::new();
        for &eidx in edge_indices {
            let e = &kg.edges[eidx];
            gene_map
                .entry(e.source.clone())
                .or_insert((e.sign, e.weight));
        }

        if gene_map.len() < 2 {
            continue;
        }
        total_tested_terms += 1;

        let mut signed_cnt = 0usize;
        let mut signed_sum = 0i32;
        for (s, _) in gene_map.values() {
            if *s != EdgeSign::Ambiguous {
                signed_cnt += 1;
                signed_sum += i32::from(s.as_i8());
            }
        }
        let mu_1 = if signed_cnt > 0 {
            (signed_sum as f64) / (signed_cnt as f64)
        } else {
            0.0
        };

        let mut overlap_genes = Vec::new();
        let mut z_obs = Vec::new();
        for (g, &(sign, weight)) in &gene_map {
            if let Some(meas) = active.degs.get(g) {
                overlap_genes.push(g.clone());
                if sign != EdgeSign::Ambiguous {
                    let s_data = if meas.log2_fc > 0.0 { 1 } else { -1 };
                    z_obs.push((weight, sign.as_i8(), s_data));
                }
            }
        }

        if overlap_genes.len() < 2 {
            continue;
        }

        let p_val = fishers_exact_right_tail(overlap_genes.len(), gene_map.len(), n_dataset, n_univ);
        let (z, _, _) = kramer_activation_z(&z_obs, mu_1, mu_2);

        let (term_name, category, kind, source_db) = kg
            .nodes
            .get(proc_id)
            .map(|n| {
                (
                    n.display_name.clone(),
                    n.family_category.clone(),
                    n.kind,
                    n.provenance_source.clone(),
                )
            })
            .unwrap_or_else(|| {
                (
                    proc_id.clone(),
                    "Biological Process".to_string(),
                    EntityKind::BiologicalProcess,
                    "GO / Open Targets".to_string(),
                )
            });

        let predicted_direction = match PredictedState::from_z(z) {
            PredictedState::Activated => "Increased".to_string(),
            PredictedState::Inhibited => "Decreased".to_string(),
            PredictedState::Neutral => "Affected".to_string(),
        };

        overlap_genes.sort();
        p_vals.push(p_val);
        hits.push(DownstreamEffectHit {
            term_id: proc_id.clone(),
            term_name,
            category,
            kind,
            source_db,
            predicted_direction,
            activation_z: z,
            p_value: p_val,
            fdr: 1.0,
            neg_log10_p: -p_val.max(1e-300).log10(),
            overlap_genes,
        });
    }

    let q_vals = benjamini_hochberg_with_total(&p_vals, total_tested_terms);
    for (h, q) in hits.iter_mut().zip(q_vals) {
        h.fdr = q;
    }

    hits.sort_by(|a, b| {
        a.p_value
            .partial_cmp(&b.p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.activation_z
                    .abs()
                    .partial_cmp(&a.activation_z.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    hits
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TreeMapRect {
    pub term_id: String,
    pub term_name: String,
    pub category: String,
    pub x: f32,
    pub y: f32,
    pub w: f32,
    pub h: f32,
    pub activation_z: f64,
    pub neg_log10_p: f64,
    pub overlap_count: usize,
}

fn slice_dice_layout(
    items: &[(usize, f64)],
    x: f32,
    y: f32,
    w: f32,
    h: f32,
    horizontal: bool,
    out: &mut Vec<(usize, f32, f32, f32, f32)>,
) {
    if items.is_empty() || w <= 0.0 || h <= 0.0 {
        return;
    }
    if items.len() == 1 {
        out.push((items[0].0, x, y, w, h));
        return;
    }
    let total: f64 = items.iter().map(|it| it.1).sum();
    if total <= 0.0 {
        return;
    }

    let mut acc = 0.0f64;
    let half = total * 0.5;
    let mut split_idx = 1;
    for (i, &(_, weight)) in items.iter().enumerate() {
        acc += weight;
        if acc >= half {
            split_idx = (i + 1).min(items.len() - 1);
            break;
        }
    }

    let left_sum: f64 = items[..split_idx].iter().map(|it| it.1).sum();
    let ratio = (left_sum / total).clamp(0.05, 0.95) as f32;

    if horizontal {
        let w1 = w * ratio;
        slice_dice_layout(&items[..split_idx], x, y, w1, h, !horizontal, out);
        slice_dice_layout(&items[split_idx..], x + w1, y, w - w1, h, !horizontal, out);
    } else {
        let h1 = h * ratio;
        slice_dice_layout(&items[..split_idx], x, y, w, h1, !horizontal, out);
        slice_dice_layout(&items[split_idx..], x, y + h1, w, h - h1, !horizontal, out);
    }
}

/// Computes a two-level hierarchical TreeMap (grouped by functional category, sized by `-log10(p)`).
pub fn compute_treemap(
    effects: &[DownstreamEffectHit],
    width: f32,
    height: f32,
) -> Vec<TreeMapRect> {
    if effects.is_empty() {
        return Vec::new();
    }

    let mut by_cat: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, hit) in effects.iter().take(40).enumerate() {
        by_cat.entry(hit.category.clone()).or_default().push(idx);
    }

    let cat_keys: Vec<String> = by_cat.keys().cloned().collect();
    let mut cat_weights: Vec<(usize, f64)> = cat_keys
        .iter()
        .enumerate()
        .map(|(cidx, k)| {
            let sum_w: f64 = by_cat[k]
                .iter()
                .map(|&i| effects[i].neg_log10_p.max(0.5))
                .sum();
            (cidx, sum_w)
        })
        .collect();
    cat_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut cat_boxes = Vec::new();
    slice_dice_layout(&cat_weights, 0.0, 0.0, width, height, width >= height, &mut cat_boxes);

    let mut rects = Vec::new();
    for (cidx, cx, cy, cw, ch) in cat_boxes {
        let cat_name = &cat_keys[cidx];
        let mut leaf_weights: Vec<(usize, f64)> = by_cat[cat_name]
            .iter()
            .map(|&i| (i, effects[i].neg_log10_p.max(0.5)))
            .collect();
        leaf_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut leaf_boxes = Vec::new();
        slice_dice_layout(&leaf_weights, cx, cy, cw, ch, cw >= ch, &mut leaf_boxes);

        for (hit_idx, lx, ly, lw, lh) in leaf_boxes {
            let hit = &effects[hit_idx];
            rects.push(TreeMapRect {
                term_id: hit.term_id.clone(),
                term_name: hit.term_name.clone(),
                category: hit.category.clone(),
                x: lx,
                y: ly,
                w: lw,
                h: lh,
                activation_z: hit.activation_z,
                neg_log10_p: hit.neg_log10_p,
                overlap_count: hit.overlap_genes.len(),
            });
        }
    }

    rects
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegulatorEffectHit {
    pub regulators: Vec<(String, f64)>,
    pub shared_dataset_targets: Vec<(String, f64)>,
    pub downstream_phenotypes: Vec<(String, String, f64)>,
    pub consistency_score: f64,
}

/// Regulator Effects engine: merges top Upstream Regulators with Downstream Effects through
/// shared differentially expressed target molecules and ranks causal 3-tier networks by Consistency Score.
pub fn run_regulator_effects(
    ura_hits: &[UpstreamRegulatorHit],
    dea_hits: &[DownstreamEffectHit],
    active: &ActiveDataset,
    _kg: &KnowledgeGraph,
) -> Vec<RegulatorEffectHit> {
    let sig_ura: Vec<&UpstreamRegulatorHit> = ura_hits
        .iter()
        .filter(|u| u.predicted_state != PredictedState::Neutral && u.p_value <= 0.05)
        .take(15)
        .collect();

    let sig_dea: Vec<&DownstreamEffectHit> = dea_hits
        .iter()
        .filter(|d| d.activation_z.abs() >= 1.5 && d.p_value <= 0.05)
        .take(15)
        .collect();

    let mut results = Vec::new();

    for u in &sig_ura {
        let u_consistent_genes: BTreeSet<&str> = u
            .target_evidences
            .iter()
            .filter(|t| t.consistent)
            .map(|t| t.target_gene.as_str())
            .collect();

        if u_consistent_genes.len() < 2 {
            continue;
        }

        let mut matched_phenos = Vec::new();
        let mut bridged_genes = BTreeSet::new();

        for d in &sig_dea {
            let shared: Vec<&str> = d
                .overlap_genes
                .iter()
                .map(String::as_str)
                .filter(|g| u_consistent_genes.contains(g))
                .collect();

            if shared.len() >= 2 {
                for g in shared {
                    bridged_genes.insert(g.to_string());
                }
                matched_phenos.push((d.term_id.clone(), d.term_name.clone(), d.activation_z));
            }
        }

        if matched_phenos.is_empty() {
            continue;
        }

        matched_phenos.sort_by(|a, b| {
            b.2.abs()
                .partial_cmp(&a.2.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        matched_phenos.truncate(4);

        let target_tuples: Vec<(String, f64)> = bridged_genes
            .into_iter()
            .filter_map(|g| active.degs.get(&g).map(|m| (g, m.log2_fc)))
            .collect();

        let pheno_z_sum: f64 = matched_phenos.iter().map(|p| p.2.abs()).sum();
        let consistency_score =
            (u.activation_z.abs() + pheno_z_sum) * (target_tuples.len() as f64).sqrt();

        results.push(RegulatorEffectHit {
            regulators: vec![(u.regulator_id.clone(), u.activation_z)],
            shared_dataset_targets: target_tuples,
            downstream_phenotypes: matched_phenos,
            consistency_score,
        });
    }

    results.sort_by(|a, b| {
        b.consistency_score
            .partial_cmp(&a.consistency_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}
