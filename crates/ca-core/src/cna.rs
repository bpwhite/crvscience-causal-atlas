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
use crate::stats::{fishers_exact_right_tail, kramer_activation_z, SplitMix64};
use crate::ura::PredictedState;
use crate::CoreAnalysisConfig;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalPathHop {
    pub nodes: Vec<String>,
    pub net_sign: i8,
    pub observed_log2fc: f64,
    pub consistent: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalNetworkHit {
    pub root_regulator: String,
    pub root_name: String,
    pub kind: EntityKind,
    pub max_depth_k: usize,
    pub predicted_state: PredictedState,
    pub activation_z: f64,
    pub overlap_p_value: f64,
    pub network_bias_p_value: f64,
    pub target_count: usize,
    pub intermediate_regulators: Vec<String>,
    pub sample_paths: Vec<CausalPathHop>,
}

fn log2_bin(in_deg: usize) -> usize {
    if in_deg <= 1 {
        0
    } else {
        (usize::BITS - (in_deg).leading_zeros() - 1) as usize
    }
}

/// Multi-hop Causal Network Analysis (CNA) per Krämer et al. (2014), Sections 3.5 & 3.6:
/// - Enumerates shortest paths up to length `K` (`1..=3`) across `E_a U E_t` (final link in `E_t`, weight >= `delta`).
/// - Collapses consistent multi-hop paths into signed virtual edges `(r, v)`.
/// - Prunes nested hypotheses via Occam's razor (`p(G_K) < p(G_{k'})` and strictly larger target coverage).
/// - Computes `log2` in-degree binned Monte Carlo network-bias-corrected p-values in parallel via `rayon`.
pub fn run_causal_networks(
    active: &ActiveDataset,
    kg: &KnowledgeGraph,
    config: &CoreAnalysisConfig,
) -> Vec<CausalNetworkHit> {
    let max_k = config.cna_max_depth.clamp(1, 3);
    let delta = config.cna_weight_cutoff;
    let universe: Vec<String> = kg.regulated_gene_universe.iter().cloned().collect();
    let n_univ = universe.len().max(1);
    let dataset_in_univ: BTreeSet<String> = active
        .degs
        .keys()
        .filter(|g| kg.regulated_gene_universe.contains(*g))
        .cloned()
        .collect();
    let n_dataset = dataset_in_univ.len();
    if n_dataset < 2 {
        return Vec::new();
    }
    let mu_2 = active.summary().mean_sign_skew;

    let mut bins_universe: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut bins_dataset_counts: BTreeMap<usize, usize> = BTreeMap::new();
    let mut gene_to_u_idx: BTreeMap<String, usize> = BTreeMap::new();

    for (u_idx, g) in universe.iter().enumerate() {
        gene_to_u_idx.insert(g.clone(), u_idx);
        let deg = kg.in_degree_t.get(g).copied().unwrap_or(1);
        let b = log2_bin(deg);
        bins_universe.entry(b).or_default().push(u_idx);
        if dataset_in_univ.contains(g) {
            *bins_dataset_counts.entry(b).or_insert(0) += 1;
        }
    }

    let candidate_roots: Vec<&String> = kg
        .nodes
        .iter()
        .filter(|(id, node)| {
            node.kind.is_molecular_regulator() && kg.out_edges_at.contains_key(*id)
        })
        .map(|(id, _)| id)
        .collect();

    let perms = config.cna_permutations.clamp(100, 10_000);

    let mut all_hits: Vec<CausalNetworkHit> = candidate_roots
        .par_iter()
        .flat_map(|&root_id| {
            let mut best_gene_depth: BTreeMap<String, usize> = BTreeMap::new();
            let mut gene_signs_at_best: BTreeMap<String, BTreeSet<i8>> = BTreeMap::new();
            let mut gene_witness_path: BTreeMap<String, Vec<String>> = BTreeMap::new();

            let mut best_inter_depth: BTreeMap<String, usize> = BTreeMap::new();
            best_inter_depth.insert(root_id.clone(), 0);

            let mut queue: VecDeque<(String, usize, i8, Vec<String>)> = VecDeque::new();
            queue.push_back((root_id.clone(), 0, 1, vec![root_id.clone()]));

            while let Some((curr, depth, prod_sign, path)) = queue.pop_front() {
                if depth >= max_k {
                    continue;
                }
                let Some(edge_indices) = kg.out_edges_at.get(&curr) else {
                    continue;
                };
                let next_depth = depth + 1;

                for &eidx in edge_indices {
                    let edge = &kg.edges[eidx];
                    if edge.weight < delta || edge.sign == EdgeSign::Ambiguous {
                        continue;
                    }
                    let step_sign = prod_sign * edge.sign.as_i8();
                    let next_node = &edge.target;
                    if path.contains(next_node) {
                        continue;
                    }

                    if edge.edge_class == EdgeClass::Transcription {
                        let prev_best = best_gene_depth.get(next_node).copied();
                        if prev_best.is_none() || prev_best == Some(next_depth) {
                            best_gene_depth.insert(next_node.clone(), next_depth);
                            gene_signs_at_best
                                .entry(next_node.clone())
                                .or_default()
                                .insert(step_sign);
                            gene_witness_path.entry(next_node.clone()).or_insert_with(|| {
                                let mut p = path.clone();
                                p.push(next_node.clone());
                                p
                            });
                        }
                    }

                    if next_depth < max_k && kg.out_edges_at.contains_key(next_node) {
                        let prev_inter = best_inter_depth.get(next_node).copied();
                        if prev_inter.is_none() || prev_inter == Some(next_depth) {
                            best_inter_depth.insert(next_node.clone(), next_depth);
                            let mut next_path = path.clone();
                            next_path.push(next_node.clone());
                            queue.push_back((next_node.clone(), next_depth, step_sign, next_path));
                        }
                    }
                }
            }

            let mut root_hypotheses: Vec<CausalNetworkHit> = Vec::new();
            let mut prev_best_p = 1.0f64;
            let mut prev_overlap_count = 0usize;

            for k in 1..=max_k {
                let mut virtual_targets: Vec<(String, i8)> = Vec::new();
                for (g, &d) in &best_gene_depth {
                    if d <= k {
                        if let Some(signs) = gene_signs_at_best.get(g) {
                            if signs.len() == 1 {
                                let net_s = *signs.iter().next().unwrap_or(&0);
                                if net_s != 0 {
                                    virtual_targets.push((g.clone(), net_s));
                                }
                            }
                        }
                    }
                }

                if virtual_targets.len() < 3 {
                    continue;
                }

                let mut obs = Vec::new();
                let mut overlap_genes = Vec::new();
                let mut sign_sum = 0i32;

                for &(ref g, net_s) in &virtual_targets {
                    sign_sum += i32::from(net_s);
                    if let Some(meas) = active.degs.get(g) {
                        let s_data = if meas.log2_fc > 0.0 { 1 } else { -1 };
                        obs.push((1.0, net_s, s_data));
                        overlap_genes.push((g.clone(), net_s, meas.log2_fc));
                    }
                }

                let overlap_cnt = overlap_genes.len();
                if overlap_cnt < 3 {
                    continue;
                }

                let p_val =
                    fishers_exact_right_tail(overlap_cnt, virtual_targets.len(), n_dataset, n_univ);
                if k > 1 && (overlap_cnt <= prev_overlap_count || p_val >= prev_best_p * 0.95) {
                    continue;
                }

                let mu_1 = (sign_sum as f64) / (virtual_targets.len() as f64);
                let (z, _, _) = kramer_activation_z(&obs, mu_1, mu_2);
                if z.abs() < 1.5 && p_val > 0.05 {
                    continue;
                }

                prev_best_p = prev_best_p.min(p_val);
                prev_overlap_count = overlap_cnt;

                let mut is_target_mask = vec![false; n_univ];
                for (g, _) in &virtual_targets {
                    if let Some(&uidx) = gene_to_u_idx.get(g) {
                        is_target_mask[uidx] = true;
                    }
                }

                let seed = root_id
                    .bytes()
                    .fold(0x517C_C1B7_2722_0A95u64, |acc, b| {
                        acc.wrapping_mul(31).wrapping_add(u64::from(b))
                    })
                    ^ (k as u64);
                let mut rng = SplitMix64::new(seed);
                let mut ge_count = 0usize;

                for _ in 0..perms {
                    let mut sim_overlap = 0usize;
                    for (&bin_id, &draw_cnt) in &bins_dataset_counts {
                        if let Some(pool) = bins_universe.get(&bin_id) {
                            let pool_len = pool.len();
                            if draw_cnt >= pool_len {
                                for &u in pool {
                                    if is_target_mask[u] {
                                        sim_overlap += 1;
                                    }
                                }
                            } else {
                                let mut chosen = BTreeSet::new();
                                while chosen.len() < draw_cnt {
                                    let idx = rng.gen_index(pool_len);
                                    if chosen.insert(idx) && is_target_mask[pool[idx]] {
                                        sim_overlap += 1;
                                    }
                                }
                            }
                        }
                    }
                    if sim_overlap >= overlap_cnt {
                        ge_count += 1;
                    }
                }

                let net_bias_p = ((ge_count + 1) as f64) / ((perms + 1) as f64);

                let mut intermediates = BTreeSet::new();
                let mut sample_paths = Vec::new();
                let root_state_sign = if z >= 0.0 { 1 } else { -1 };

                for (g, net_s, fc) in overlap_genes {
                    if let Some(wpath) = gene_witness_path.get(&g) {
                        for mid in wpath.iter().skip(1).take(wpath.len().saturating_sub(2)) {
                            intermediates.insert(mid.clone());
                        }
                        let s_data = if fc > 0.0 { 1 } else { -1 };
                        let consistent = (net_s * root_state_sign) == s_data;
                        sample_paths.push(CausalPathHop {
                            nodes: wpath.clone(),
                            net_sign: net_s,
                            observed_log2fc: fc,
                            consistent,
                        });
                    }
                }

                sample_paths.sort_by(|a, b| {
                    b.nodes.len().cmp(&a.nodes.len()).then_with(|| {
                        b.observed_log2fc
                            .abs()
                            .partial_cmp(&a.observed_log2fc.abs())
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                });
                sample_paths.truncate(24);

                let (root_name, kind) = kg
                    .nodes
                    .get(root_id)
                    .map(|n| (n.display_name.clone(), n.kind))
                    .unwrap_or_else(|| (root_id.clone(), EntityKind::GeneOther));

                root_hypotheses.push(CausalNetworkHit {
                    root_regulator: root_id.clone(),
                    root_name,
                    kind,
                    max_depth_k: k,
                    predicted_state: PredictedState::from_z(z),
                    activation_z: z,
                    overlap_p_value: p_val,
                    network_bias_p_value: net_bias_p,
                    target_count: overlap_cnt,
                    intermediate_regulators: intermediates.into_iter().collect(),
                    sample_paths,
                });
            }

            root_hypotheses
        })
        .collect();

    all_hits.sort_by(|a, b| {
        b.activation_z
            .abs()
            .partial_cmp(&a.activation_z.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.network_bias_p_value
                    .partial_cmp(&b.network_bias_p_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    all_hits
}
