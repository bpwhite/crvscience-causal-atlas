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

//! # `ca-core` — Clean-Room Causal Pathway & Network Inference Engine
//!
//! Pure-Rust statistical and graph-theoretic algorithms for interpreting high-throughput
//! transcriptomic and proteomic differential expression datasets over directed, signed
//! biological knowledge graphs assembled exclusively from open-source data streams.
//!
//! ## Clean-Room Independence
//! - Zero proprietary SDKs, APIs, or database schemas.
//! - Zero FFI linking to C, C++, R/Bioconductor, Python, or Java runtimes (`#![forbid(unsafe_code)]`).
//! - Implements published, peer-reviewed causal inference formulations:
//!   - Krämer, A., Green, J., Pollard, J. Jr., & Tugendreich, S. (2014). Causal analysis
//!     approaches in Ingenuity Pathway Analysis. *Bioinformatics*, 30(4), 523–530.
//!     <https://pmc.ncbi.nlm.nih.gov/articles/PMC3928520/>
//!   - Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate. *JRSS-B*.
//!   - Shneiderman, B. (1992). Tree visualization with tree-maps: 2-d space-filling approach.

pub mod canonical;
pub mod cna;
pub mod dataset;
pub mod dea;
pub mod graph;
pub mod layout;
pub mod perturb;
pub mod seed_kb;
pub mod stats;
pub mod ura;

use serde::{Deserialize, Serialize};

pub use canonical::{
    run_canonical_pathways, simulate_map, simulate_map_in_kg, CanonicalPathwayHit, MapNodeState,
    MapSimulationResult,
};
pub use cna::{run_causal_networks, CausalNetworkHit, CausalPathHop};
pub use dataset::{
    builtin_cohorts, parse_expression_table, ActiveDataset, DatasetSummary, ExpressionDataset,
    GeneMeasurement,
};
pub use dea::{
    compute_treemap, run_downstream_effects, run_regulator_effects, DownstreamEffectHit,
    RegulatorEffectHit, TreeMapRect,
};
pub use graph::{
    CausalEdge, EdgeClass, EdgeSign, EntityKind, EntityNode, KnowledgeGraph, PathwayDef,
    PathwayReactionEdge, PerturbationKind, PerturbationSignature, SubcellularLocation,
};
pub use layout::{
    compute_pathway_world_layout_cached, compute_pathway_world_layout_uncached,
    compute_tiered_causal_layout, compute_world_canvas_size, induce_pathway_topology_edges,
    recompute_and_cache_pathway_layout, resolve_normalized_pathway_coords, CompartmentBands,
    PathwayWorldLayout, TieredNodeSpec, TieredWorldLayout, MIN_NODE_DX, MIN_NODE_DY,
};
pub use perturb::{
    compute_regulator_rank_diff, format_regulator_rank_diff_csv, format_regulator_rank_diff_json,
    run_comparison_analysis, run_perturbation_match, ComparisonMatrix, ComparisonMetric,
    PerturbationHit, PerturbationRankDiffDigest, RegulatorRankDiffEntry,
};
pub use ura::{
    build_mechanistic_network, format_z_2dp, round_z_2dp, run_upstream_regulators,
    select_stratified_mechanistic_targets, MechanisticEdge, MechanisticNetwork, PredictedState,
    UpstreamRegulatorHit,
};

/// Tunable thresholds and algorithmic parameters for a full Core Analysis run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoreAnalysisConfig {
    /// Minimum absolute log2 fold-change to include a gene in the active differentially expressed set.
    pub abs_log2fc_cutoff: f64,
    /// Maximum Benjamini-Hochberg adjusted p-value (FDR) to include a gene in the active set.
    pub fdr_cutoff: f64,
    /// Minimum edge confidence weight `delta` for multi-hop Causal Network Analysis (default 0.25).
    pub cna_weight_cutoff: f64,
    /// Maximum path length `K` (1..=3) for Causal Network Analysis root regulator search.
    pub cna_max_depth: usize,
    /// Number of degree-preserving Monte Carlo permutations for CNA network-bias p-values.
    pub cna_permutations: usize,
    /// Maximum regulator breadth `N` per hop when constructing Mechanistic Networks.
    pub mechanistic_breadth: usize,
    /// Maximum hop depth `K` when constructing Mechanistic Networks.
    pub mechanistic_depth: usize,
    /// Damping factor for Molecule Activity Predictor (MAP) iterative signed message passing.
    pub map_damping: f64,
    /// Number of propagation iterations for MAP simulation.
    pub map_iterations: usize,
}

impl Default for CoreAnalysisConfig {
    fn default() -> Self {
        Self {
            abs_log2fc_cutoff: 1.0,
            fdr_cutoff: 0.05,
            cna_weight_cutoff: 0.25,
            cna_max_depth: 3,
            cna_permutations: 1_000,
            mechanistic_breadth: 5,
            mechanistic_depth: 3,
            map_damping: 0.82,
            map_iterations: 16,
        }
    }
}

/// Complete output bundle produced by a Core Analysis execution on a single expression dataset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoreAnalysisReport {
    pub dataset_name: String,
    pub config: CoreAnalysisConfig,
    pub summary: DatasetSummary,
    pub canonical_pathways: Vec<CanonicalPathwayHit>,
    pub upstream_regulators: Vec<UpstreamRegulatorHit>,
    pub causal_networks: Vec<CausalNetworkHit>,
    pub downstream_effects: Vec<DownstreamEffectHit>,
    pub treemap_rects: Vec<TreeMapRect>,
    pub regulator_effects: Vec<RegulatorEffectHit>,
    pub perturbation_matches: Vec<PerturbationHit>,
}

/// Executes all causal and enrichment modules over `dataset` and `kg` in a single coordinated run.
pub fn run_core_analysis(
    dataset: &ExpressionDataset,
    kg: &KnowledgeGraph,
    config: &CoreAnalysisConfig,
) -> CoreAnalysisReport {
    let active = dataset.filter_active(config.abs_log2fc_cutoff, config.fdr_cutoff, kg);
    let summary = active.summary();

    let canonical_pathways = run_canonical_pathways(&active, kg);
    let upstream_regulators = run_upstream_regulators(&active, kg, config);
    let causal_networks = run_causal_networks(&active, kg, config);
    let downstream_effects = run_downstream_effects(&active, kg);
    let treemap_rects = compute_treemap(&downstream_effects, 1.0, 1.0);
    let regulator_effects =
        run_regulator_effects(&upstream_regulators, &downstream_effects, &active, kg);
    let perturbation_matches = run_perturbation_match(&active, kg);

    CoreAnalysisReport {
        dataset_name: dataset.name.clone(),
        config: config.clone(),
        summary,
        canonical_pathways,
        upstream_regulators,
        causal_networks,
        downstream_effects,
        treemap_rects,
        regulator_effects,
        perturbation_matches,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_to_end_core_analysis_identifies_expected_drivers() {
        let kg = seed_kb::build_seed_knowledge_graph();
        let cohorts = dataset::builtin_cohorts();
        let tnf_ds = &cohorts[0];
        let cfg = CoreAnalysisConfig {
            cna_permutations: 250,
            ..CoreAnalysisConfig::default()
        };
        let report = run_core_analysis(tnf_ds, &kg, &cfg);

        assert!(report.summary.active_genes >= 15);
        assert!(!report.canonical_pathways.is_empty());
        assert!(!report.upstream_regulators.is_empty());
        assert!(!report.causal_networks.is_empty());
        assert!(!report.downstream_effects.is_empty());
        assert!(!report.perturbation_matches.is_empty());

        // TNF and NFKB1 should rank near the top of activated upstream regulators for the TNF cohort.
        let tnf_hit = report
            .upstream_regulators
            .iter()
            .find(|r| r.regulator_id == "TNF")
            .expect("TNF regulator must be scored");
        assert!(
            tnf_hit.activation_z > 2.0,
            "Expected TNF z-score > 2.0, got {}",
            tnf_hit.activation_z
        );
        assert_eq!(tnf_hit.predicted_state, PredictedState::Activated);
    }
}
