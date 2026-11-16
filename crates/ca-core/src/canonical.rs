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
use crate::graph::{KnowledgeGraph, PathwayDef};
use crate::stats::{benjamini_hochberg, fishers_exact_right_tail, kramer_activation_z};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CanonicalPathwayHit {
    pub pathway_id: String,
    pub pathway_name: String,
    pub source_db: String,
    pub category: String,
    pub url: String,
    pub pathway_size: usize,
    pub overlap_count: usize,
    pub overlap_ratio: f64,
    pub p_value: f64,
    pub fdr: f64,
    pub neg_log10_fdr: f64,
    pub activation_z: Option<f64>,
    pub upregulated_genes: Vec<String>,
    pub downregulated_genes: Vec<String>,
}

/// Executes Over-Representation Analysis (one-sided Fisher's exact test + BH-FDR) and
/// signed pathway Activation Z-score across all canonical pathways in `kg`.
pub fn run_canonical_pathways(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
) -> Vec<CanonicalPathwayHit> {
    let universe = &kg.pathway_gene_universe;
    let n_univ = universe.len().max(1);
    let dataset_in_univ: BTreeSet<&String> = active
        .degs
        .keys()
        .filter(|g| universe.contains(*g))
        .collect();
    let n_dataset = dataset_in_univ.len();
    let ds_skew = active.summary().mean_sign_skew;

    let mut raw_hits = Vec::new();
    let mut p_vals = Vec::new();

    for pw in &kg.pathways {
        let unique_members: BTreeSet<&String> = pw
            .members
            .iter()
            .filter(|m| universe.contains(*m))
            .collect();
        let pw_size = unique_members.len();
        if pw_size == 0 {
            continue;
        }

        let mut up_genes = Vec::new();
        let mut down_genes = Vec::new();
        let mut signed_obs = Vec::new();

        for &m in &unique_members {
            if let Some(meas) = active.degs.get(m) {
                if meas.log2_fc > 0.0 {
                    up_genes.push(m.clone());
                } else if meas.log2_fc < 0.0 {
                    down_genes.push(m.clone());
                }
                let s_data = if meas.log2_fc > 0.0 { 1 } else { -1 };
                let s_exp = pw.member_expected_signs.get(m).copied().unwrap_or(1);
                if s_exp != 0 {
                    signed_obs.push((1.0, s_exp, s_data));
                }
            }
        }

        let overlap = up_genes.len() + down_genes.len();
        if overlap == 0 {
            continue;
        }

        let p_val = fishers_exact_right_tail(overlap, pw_size, n_dataset, n_univ);
        let mean_pw_sign = if pw.member_expected_signs.is_empty() {
            1.0
        } else {
            let s: i32 = pw.member_expected_signs.values().map(|&v| i32::from(v)).sum();
            (s as f64) / (pw.member_expected_signs.len() as f64)
        };

        let activation_z = if signed_obs.len() >= 2 {
            let (z, _, _) = kramer_activation_z(&signed_obs, mean_pw_sign, ds_skew);
            Some(z)
        } else {
            None
        };

        up_genes.sort();
        down_genes.sort();
        p_vals.push(p_val);

        raw_hits.push(CanonicalPathwayHit {
            pathway_id: pw.id.clone(),
            pathway_name: pw.name.clone(),
            source_db: pw.source_db.clone(),
            category: pw.category.clone(),
            url: pw.url.clone(),
            pathway_size: pw_size,
            overlap_count: overlap,
            overlap_ratio: (overlap as f64) / (pw_size as f64),
            p_value: p_val,
            fdr: 1.0,
            neg_log10_fdr: 0.0,
            activation_z,
            upregulated_genes: up_genes,
            downregulated_genes: down_genes,
        });
    }

    let q_vals = benjamini_hochberg(&p_vals);
    for (hit, q) in raw_hits.iter_mut().zip(q_vals) {
        hit.fdr = q;
        hit.neg_log10_fdr = -q.max(1e-300).log10();
    }

    raw_hits.sort_by(|a, b| {
        a.p_value
            .partial_cmp(&b.p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.activation_z
                    .unwrap_or(0.0)
                    .abs()
                    .partial_cmp(&a.activation_z.unwrap_or(0.0).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    raw_hits
}

/// Node activity state during Molecule Activity Predictor (MAP) topological propagation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MapNodeState {
    pub node_id: String,
    pub observed_log2fc: Option<f64>,
    pub user_clamp: Option<f64>,
    /// Final signed activity in `[-1.0, +1.0]` (`> +0.15` predicted activated, `< -0.15` predicted inhibited).
    pub predicted_activity: f64,
    pub is_inconsistent: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MapSimulationResult {
    pub pathway_id: String,
    pub states: BTreeMap<String, MapNodeState>,
    /// Edge index -> `(is_consistent, propagated_flux)`
    pub edge_flux: Vec<(bool, f64)>,
}

/// Molecule Activity Predictor (MAP): propagates observed expression states and user in-silico
/// clamps (`+1.0` activate / `-1.0` inhibit) both downstream and upstream along signed pathway topology edges.
pub fn simulate_map(
    pathway: &PathwayDef,
    active: &ActiveDataset,
    user_clamps: &BTreeMap<String, f64>,
    damping: f64,
    iterations: usize,
) -> MapSimulationResult {
    let mut all_nodes: BTreeSet<String> = pathway.members.iter().cloned().collect();
    for e in &pathway.topology_edges {
        all_nodes.insert(e.from.clone());
        all_nodes.insert(e.to.clone());
    }

    let mut states: BTreeMap<String, MapNodeState> = BTreeMap::new();
    for nid in &all_nodes {
        let obs = active.degs.get(nid).map(|m| m.log2_fc);
        let clamp = user_clamps.get(nid).copied();
        let init_val = if let Some(c) = clamp {
            c.clamp(-1.0, 1.0)
        } else if let Some(fc) = obs {
            (fc / 3.0).tanh()
        } else {
            0.0
        };

        states.insert(
            nid.clone(),
            MapNodeState {
                node_id: nid.clone(),
                observed_log2fc: obs,
                user_clamp: clamp,
                predicted_activity: init_val,
                is_inconsistent: false,
            },
        );
    }

    let damp = damping.clamp(0.2, 0.95);
    for _ in 0..iterations.max(1) {
        let prev: BTreeMap<String, f64> = states
            .iter()
            .map(|(k, v)| (k.clone(), v.predicted_activity))
            .collect();

        for (nid, st) in &mut states {
            // User clamps remain strictly pinned; observed DEGs anchor with soft prior unless overridden.
            if let Some(c) = st.user_clamp {
                st.predicted_activity = c.clamp(-1.0, 1.0);
                continue;
            }

            let mut incoming_sum = 0.0f64;
            let mut incoming_cnt = 0.0f64;
            let mut reverse_sum = 0.0f64;
            let mut reverse_cnt = 0.0f64;

            for edge in &pathway.topology_edges {
                let s = edge.sign.as_f64();
                if s == 0.0 {
                    continue;
                }
                if &edge.to == nid {
                    let src_val = prev.get(&edge.from).copied().unwrap_or(0.0);
                    if src_val.abs() > 1e-4 {
                        incoming_sum += src_val * s;
                        incoming_cnt += 1.0;
                    }
                } else if &edge.from == nid && st.observed_log2fc.is_none() {
                    // Weak upstream retro-inference for unmeasured upstream receptors/ligands
                    let dst_val = prev.get(&edge.to).copied().unwrap_or(0.0);
                    if dst_val.abs() > 1e-4 {
                        reverse_sum += dst_val * s * 0.35;
                        reverse_cnt += 0.35;
                    }
                }
            }

            let net_signal = if incoming_cnt + reverse_cnt > 0.0 {
                (incoming_sum + reverse_sum) / (incoming_cnt + reverse_cnt)
            } else {
                0.0
            };

            if let Some(fc) = st.observed_log2fc {
                let obs_anchor = (fc / 3.0).tanh();
                st.predicted_activity = (0.75 * obs_anchor + 0.25 * damp * net_signal).clamp(-1.0, 1.0);
            } else {
                st.predicted_activity = (damp * net_signal).tanh();
            }
        }
    }

    let mut edge_flux = Vec::with_capacity(pathway.topology_edges.len());
    for edge in &pathway.topology_edges {
        let src_a = states
            .get(&edge.from)
            .map(|s| s.predicted_activity)
            .unwrap_or(0.0);
        let dst_a = states
            .get(&edge.to)
            .map(|s| s.predicted_activity)
            .unwrap_or(0.0);
        let s = edge.sign.as_f64();
        let expected_dst = src_a * s;
        let consistent = expected_dst * dst_a >= -0.05;
        if !consistent {
            if let Some(dst_node) = states.get_mut(&edge.to) {
                dst_node.is_inconsistent = true;
            }
        }
        edge_flux.push((consistent, expected_dst));
    }

    MapSimulationResult {
        pathway_id: pathway.id.clone(),
        states,
        edge_flux,
    }
}

/// Molecule Activity Predictor (MAP) with automatic `KnowledgeGraph` causal edge induction
/// so gene-set pathways (e.g. GMT imports) propagate activity across curated interactions.
pub fn simulate_map_in_kg(
    pathway: &PathwayDef,
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
    user_clamps: &BTreeMap<String, f64>,
    damping: f64,
    iterations: usize,
) -> MapSimulationResult {
    let induced = crate::layout::induce_pathway_topology_edges(pathway, kg);
    if induced.len() == pathway.topology_edges.len() {
        return simulate_map(pathway, active, user_clamps, damping, iterations);
    }
    let mut enriched = pathway.clone();
    enriched.topology_edges = induced;
    simulate_map(&enriched, active, user_clamps, damping, iterations)
}
