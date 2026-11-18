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
use crate::graph::{EdgeClass, EdgeSign, EntityKind, KnowledgeGraph};
use crate::stats::{benjamini_hochberg, fishers_exact_right_tail, kramer_activation_z};
use crate::CoreAnalysisConfig;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PredictedState {
    Activated,
    Inhibited,
    Neutral,
}

/// Rounds an activation Z-score to 2 decimal places matching `{:+.2}` display formatting.
pub fn round_z_2dp(z: f64) -> f64 {
    if !z.is_finite() {
        return 0.0;
    }
    let rounded = format!("{z:.2}").parse::<f64>().unwrap_or(0.0);
    if rounded == 0.0 { 0.0 } else { rounded }
}

/// Formats an activation Z-score to 2 decimal places with explicit sign, normalizing `-0.00` to `+0.00`.
pub fn format_z_2dp(z: f64) -> String {
    let z_disp = round_z_2dp(z);
    format!("{z_disp:+.2}")
}

impl PredictedState {
    pub fn from_z(z: f64) -> Self {
        let z_disp = round_z_2dp(z);
        if z_disp >= 2.0 {
            Self::Activated
        } else if z_disp <= -2.0 {
            Self::Inhibited
        } else {
            Self::Neutral
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Activated => "Activated",
            Self::Inhibited => "Inhibited",
            Self::Neutral => "Unresolved",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetEvidence {
    pub target_gene: String,
    pub edge_sign: EdgeSign,
    pub observed_log2fc: f64,
    pub consistent: bool,
    pub source_db: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UpstreamRegulatorHit {
    pub regulator_id: String,
    pub regulator_name: String,
    pub kind: EntityKind,
    pub observed_log2fc: Option<f64>,
    pub predicted_state: PredictedState,
    pub activation_z: f64,
    pub bias_corrected_z: f64,
    pub skew_flagged: bool,
    pub overlap_count: usize,
    pub regulon_size: usize,
    pub p_value: f64,
    pub fdr: f64,
    pub target_evidences: Vec<TargetEvidence>,
    pub mechanistic_summary: Option<String>,
}

/// Implements Upstream Regulator Analysis (URA) per Krämer et al. (2014), Section 3.3.
pub fn run_upstream_regulators(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
    config: &CoreAnalysisConfig,
) -> Vec<UpstreamRegulatorHit> {
    let universe = &kg.regulated_gene_universe;
    let n_univ = universe.len().max(1);
    let dataset_in_univ: BTreeSet<&String> = active
        .degs
        .keys()
        .filter(|g| universe.contains(*g))
        .collect();
    let n_dataset = dataset_in_univ.len();
    let mu_2 = active.summary().mean_sign_skew;

    let mut hits = Vec::new();
    let mut p_vals = Vec::new();

    for (reg_id, edge_indices) in &kg.out_edges_t {
        if let Some(node) = kg.nodes.get(reg_id) {
            if !node.kind.is_molecular_regulator() {
                continue;
            }
        }

        // Deduplicate targets regulated by `reg_id` along T-edges.
        let mut target_map: BTreeMap<String, (EdgeSign, f64, String)> = BTreeMap::new();
        for &eidx in edge_indices {
            let e = &kg.edges[eidx];
            target_map
                .entry(e.target.clone())
                .and_modify(|cur| {
                    if e.weight > cur.1 {
                        *cur = (e.sign, e.weight, e.source_db.clone());
                    }
                })
                .or_insert_with(|| (e.sign, e.weight, e.source_db.clone()));
        }

        let regulon_size = target_map.len();
        if regulon_size < 2 {
            continue;
        }

        let mut signed_total = 0usize;
        let mut signed_sum = 0i32;
        for (sign, _, _) in target_map.values() {
            if *sign != EdgeSign::Ambiguous {
                signed_total += 1;
                signed_sum += i32::from(sign.as_i8());
            }
        }
        let mu_1 = if signed_total > 0 {
            (signed_sum as f64) / (signed_total as f64)
        } else {
            0.0
        };

        let mut overlap_count = 0usize;
        let mut z_obs = Vec::new();
        let mut target_evidences = Vec::new();

        for (target_sym, &(sign, weight, ref sdb)) in &target_map {
            if let Some(meas) = active.degs.get(target_sym) {
                overlap_count += 1;
                let s_data = if meas.log2_fc > 0.0 { 1 } else { -1 };
                let s_edge = sign.as_i8();
                let consistent = s_edge != 0 && (i32::from(s_edge) * i32::from(s_data) > 0);
                if s_edge != 0 {
                    z_obs.push((weight, s_edge, s_data));
                }
                target_evidences.push(TargetEvidence {
                    target_gene: target_sym.clone(),
                    edge_sign: sign,
                    observed_log2fc: meas.log2_fc,
                    consistent,
                    source_db: sdb.clone(),
                });
            }
        }

        if overlap_count < 2 {
            continue;
        }

        let p_val = fishers_exact_right_tail(overlap_count, regulon_size, n_dataset, n_univ);
        let (z, z_corr, skew_flagged) = kramer_activation_z(&z_obs, mu_1, mu_2);
        let (reg_name, kind) = kg
            .nodes
            .get(reg_id)
            .map(|n| (n.display_name.clone(), n.kind))
            .unwrap_or_else(|| (reg_id.clone(), EntityKind::GeneOther));

        let obs_fc = active.degs.get(reg_id).map(|m| m.log2_fc);
        target_evidences.sort_by(|a, b| {
            b.observed_log2fc
                .abs()
                .partial_cmp(&a.observed_log2fc.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        p_vals.push(p_val);
        hits.push(UpstreamRegulatorHit {
            regulator_id: reg_id.clone(),
            regulator_name: reg_name,
            kind,
            observed_log2fc: obs_fc,
            predicted_state: PredictedState::from_z(z),
            activation_z: z,
            bias_corrected_z: z_corr,
            skew_flagged,
            overlap_count,
            regulon_size,
            p_value: p_val,
            fdr: 1.0,
            target_evidences,
            mechanistic_summary: None,
        });
    }

    let q_vals = benjamini_hochberg(&p_vals);
    for (h, q) in hits.iter_mut().zip(q_vals) {
        h.fdr = q;
    }

    hits.sort_by(|a, b| {
        b.activation_z
            .abs()
            .partial_cmp(&a.activation_z.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.p_value
                    .partial_cmp(&b.p_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    // Annotate top hits with Mechanistic Network summary counts (`targets (regulators)`).
    let snapshot = hits.clone();
    for h in hits.iter_mut().take(25) {
        let mn = build_mechanistic_network(
            &h.regulator_id,
            &snapshot,
            active,
            kg,
            config.mechanistic_breadth,
            config.mechanistic_depth,
        );
        if mn.regulator_ids.len() > 1 {
            h.mechanistic_summary = Some(format!(
                "{} targets ({} regulators)",
                mn.connected_dataset_genes.len(),
                mn.regulator_ids.len()
            ));
        }
    }

    hits
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MechanisticEdge {
    pub from_regulator: String,
    pub to_regulator: String,
    pub sign: EdgeSign,
    pub triangle_overlap: usize,
    pub triangle_p_value: f64,
    pub mechanism: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MechanisticNetwork {
    pub root_regulator: String,
    pub regulator_ids: Vec<String>,
    pub regulator_z_scores: BTreeMap<String, f64>,
    pub edges: Vec<MechanisticEdge>,
    pub connected_dataset_genes: Vec<String>,
}

/// Implements Mechanistic Networks (Krämer et al., 2014, Section 3.4) via Fisher's exact enrichment
/// of causal transitive triangles (`r1 -> r2`, `r1 -> v`, `r2 -> v` for `v in D`).
pub fn build_mechanistic_network(
    root_regulator: &str,
    ura_hits: &[UpstreamRegulatorHit],
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
    breadth_n: usize,
    depth_k: usize,
) -> MechanisticNetwork {
    let n_dataset = active.degs.len().max(1);

    let hit_map: BTreeMap<&str, &UpstreamRegulatorHit> = ura_hits
        .iter()
        .filter(|h| h.p_value <= 0.05 || h.activation_z.abs() >= 1.5)
        .map(|h| (h.regulator_id.as_str(), h))
        .collect();

    let overlap_sets: BTreeMap<&str, BTreeSet<&str>> = hit_map
        .iter()
        .map(|(&rid, &hit)| {
            let s: BTreeSet<&str> = hit
                .target_evidences
                .iter()
                .map(|t| t.target_gene.as_str())
                .collect();
            (rid, s)
        })
        .collect();

    let mut visited_regs = BTreeSet::new();
    visited_regs.insert(root_regulator.to_string());
    let mut frontier = vec![root_regulator.to_string()];
    let mut mech_edges = Vec::new();

    for _hop in 0..depth_k.max(1) {
        let mut next_frontier = Vec::new();
        for r1 in &frontier {
            let Some(o_r1) = overlap_sets.get(r1.as_str()) else {
                continue;
            };
            let Some(out_indices) = kg.out_edges_at.get(r1) else {
                continue;
            };

            let mut candidates: Vec<MechanisticEdge> = Vec::new();
            for &eidx in out_indices {
                let edge = &kg.edges[eidx];
                if edge.edge_class == EdgeClass::Phenotype {
                    continue;
                }
                let r2 = &edge.target;
                if visited_regs.contains(r2) {
                    continue;
                }
                let Some(o_r2) = overlap_sets.get(r2.as_str()) else {
                    continue;
                };

                let shared = o_r1.intersection(o_r2).count();
                if shared == 0 {
                    continue;
                }
                let tri_p = fishers_exact_right_tail(shared, o_r2.len(), o_r1.len(), n_dataset);
                if tri_p <= 0.25 {
                    candidates.push(MechanisticEdge {
                        from_regulator: r1.clone(),
                        to_regulator: r2.clone(),
                        sign: edge.sign,
                        triangle_overlap: shared,
                        triangle_p_value: tri_p,
                        mechanism: edge.mechanism.clone(),
                    });
                }
            }

            candidates.sort_by(|a, b| {
                a.triangle_p_value
                    .partial_cmp(&b.triangle_p_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            for chosen in candidates.into_iter().take(breadth_n.max(1)) {
                if visited_regs.insert(chosen.to_regulator.clone()) {
                    next_frontier.push(chosen.to_regulator.clone());
                }
                mech_edges.push(chosen);
            }
        }
        if next_frontier.is_empty() {
            break;
        }
        frontier = next_frontier;
    }

    let mut regulator_z_scores = BTreeMap::new();
    let mut connected_genes_set = BTreeSet::new();
    for rid in &visited_regs {
        if let Some(h) = hit_map.get(rid.as_str()) {
            regulator_z_scores.insert(rid.clone(), h.activation_z);
            for t in &h.target_evidences {
                connected_genes_set.insert(t.target_gene.clone());
            }
        }
    }

    MechanisticNetwork {
        root_regulator: root_regulator.to_string(),
        regulator_ids: visited_regs.into_iter().collect(),
        regulator_z_scores,
        edges: mech_edges,
        connected_dataset_genes: connected_genes_set.into_iter().collect(),
    }
}

/// Selects up to `max_total` target evidences for Mechanistic Network visualization using
/// polarity-stratified quotas so repressed/inhibited targets are preserved alongside activated
/// targets instead of being truncated by pure `|log2FC|` ordering.
pub fn select_stratified_mechanistic_targets<'a>(
    targets: &'a [TargetEvidence],
    max_total: usize,
) -> Vec<&'a TargetEvidence> {
    if max_total == 0 || targets.is_empty() {
        return Vec::new();
    }

    let mut ordered: Vec<&'a TargetEvidence> = targets.iter().collect();
    ordered.sort_by(|a, b| {
        b.observed_log2fc
            .abs()
            .total_cmp(&a.observed_log2fc.abs())
            .then_with(|| a.target_gene.cmp(&b.target_gene))
    });

    let mut activated: Vec<&'a TargetEvidence> = Vec::new();
    let mut repressed: Vec<&'a TargetEvidence> = Vec::new();
    for t in &ordered {
        let is_repressed = t.edge_sign == EdgeSign::Inhibiting
            || (t.edge_sign == EdgeSign::Ambiguous && t.observed_log2fc < 0.0);
        if is_repressed {
            repressed.push(*t);
        } else {
            activated.push(*t);
        }
    }

    // Fallback: if edge signs are one-sided while dataset log2FC signs span both directions,
    // stratify by observed expression sign so repressed targets remain visible.
    if (activated.is_empty() || repressed.is_empty())
        && ordered.iter().any(|t| t.observed_log2fc < 0.0)
        && ordered.iter().any(|t| t.observed_log2fc >= 0.0)
    {
        activated.clear();
        repressed.clear();
        for t in &ordered {
            if t.observed_log2fc < 0.0 {
                repressed.push(*t);
            } else {
                activated.push(*t);
            }
        }
    }

    let half = max_total / 2;
    let rem = max_total - half;
    let mut take_act = activated.len().min(rem);
    let mut take_rep = repressed.len().min(half);
    let leftover = max_total - (take_act + take_rep);
    if leftover > 0 {
        let extra_act = (activated.len() - take_act).min(leftover);
        take_act += extra_act;
        let rem_after = leftover - extra_act;
        if rem_after > 0 {
            take_rep += (repressed.len() - take_rep).min(rem_after);
        }
    }

    let mut out = Vec::with_capacity(take_act + take_rep);
    out.extend_from_slice(&activated[..take_act]);
    out.extend_from_slice(&repressed[..take_rep]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::builtin_cohorts;
    use crate::seed_kb::build_seed_knowledge_graph;

    #[test]
    fn predicted_state_matches_two_decimal_display_never_contradicts() {
        let cases = [
            (-1.999637, PredictedState::Inhibited),
            (1.998199, PredictedState::Activated),
            (-1.995, PredictedState::Inhibited),
            (1.995, PredictedState::Activated),
            (-1.994, PredictedState::Neutral),
            (1.994, PredictedState::Neutral),
            (-2.0, PredictedState::Inhibited),
            (2.0, PredictedState::Activated),
            (f64::NAN, PredictedState::Neutral),
            (f64::INFINITY, PredictedState::Neutral),
        ];
        for (z, expected) in cases {
            let state = PredictedState::from_z(z);
            assert_eq!(state, expected, "z={z}");
            if z.is_finite() {
                let disp = format_z_2dp(z);
                assert_ne!(disp, "-0.00", "z={z} formatted as negative zero");
                if disp == "+2.00" {
                    assert_eq!(state, PredictedState::Activated);
                } else if disp == "-2.00" {
                    assert_eq!(state, PredictedState::Inhibited);
                } else {
                    assert_eq!(state, PredictedState::Neutral);
                }
            }
        }

        // Negative-zero normalization checks
        for z in [-0.004, -0.0049, -0.0, 0.0, 0.004] {
            let r = round_z_2dp(z);
            assert_eq!(r, 0.0);
            assert!(!r.is_sign_negative(), "round_z_2dp({z}) returned -0.0");
            assert_eq!(format_z_2dp(z), "+0.00");
        }
    }

    #[test]
    fn stratified_mechanistic_targets_preserves_repressed_targets_on_tnf_cohort() {
        let kg = build_seed_knowledge_graph();
        let cohorts = builtin_cohorts();
        let cfg = CoreAnalysisConfig::default();
        let active = cohorts[0].filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
        let hits = run_upstream_regulators(&active, &kg, &cfg);
        let tnf = hits
            .iter()
            .find(|h| h.regulator_id == "TNF")
            .expect("TNF hit present");

        let selected = select_stratified_mechanistic_targets(&tnf.target_evidences, 14);
        assert_eq!(selected.len(), 14);
        let rep_syms: BTreeSet<&str> = selected
            .iter()
            .filter(|t| t.edge_sign == EdgeSign::Inhibiting)
            .map(|t| t.target_gene.as_str())
            .collect();
        for expected in ["KLF2", "NOS3", "THBD", "PPARG"] {
            assert!(
                rep_syms.contains(expected),
                "repressed target {expected} must not be dropped"
            );
        }
        let act_count = selected
            .iter()
            .filter(|t| t.edge_sign == EdgeSign::Activating)
            .count();
        assert_eq!(act_count, 10);
    }

    #[test]
    fn stratified_mechanistic_targets_splits_balanced_and_one_sided_regulons() {
        let mut balanced = Vec::new();
        for i in 0..12 {
            balanced.push(TargetEvidence {
                target_gene: format!("ACT{i:02}"),
                edge_sign: EdgeSign::Activating,
                observed_log2fc: 3.0 - (i as f64) * 0.1,
                consistent: true,
                source_db: "Test".to_string(),
            });
            balanced.push(TargetEvidence {
                target_gene: format!("REP{i:02}"),
                edge_sign: EdgeSign::Inhibiting,
                observed_log2fc: -3.0 + (i as f64) * 0.1,
                consistent: true,
                source_db: "Test".to_string(),
            });
        }
        let picked = select_stratified_mechanistic_targets(&balanced, 14);
        assert_eq!(picked.len(), 14);
        assert_eq!(
            picked
                .iter()
                .filter(|t| t.edge_sign == EdgeSign::Activating)
                .count(),
            7
        );
        assert_eq!(
            picked
                .iter()
                .filter(|t| t.edge_sign == EdgeSign::Inhibiting)
                .count(),
            7
        );

        let one_sided: Vec<TargetEvidence> = (0..16)
            .map(|i| TargetEvidence {
                target_gene: format!("INH{i:02}"),
                edge_sign: EdgeSign::Inhibiting,
                observed_log2fc: -4.0 + (i as f64) * 0.1,
                consistent: true,
                source_db: "Test".to_string(),
            })
            .collect();
        let picked_one = select_stratified_mechanistic_targets(&one_sided, 14);
        assert_eq!(picked_one.len(), 14);
        assert!(select_stratified_mechanistic_targets(&one_sided, 0).is_empty());
    }
}
