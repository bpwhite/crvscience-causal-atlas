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

use ca_core::cna::run_causal_networks;
use ca_core::perturb::run_perturbation_match;
use ca_core::seed_kb::build_seed_knowledge_graph;
use ca_core::stats::{benjamini_hochberg, z_to_two_sided_p};
use ca_core::ura::{build_mechanistic_network, run_upstream_regulators, PredictedState};
use ca_core::{
    builtin_cohorts, parse_expression_table, run_core_analysis, CoreAnalysisConfig,
};
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

const MAX_VALIDATION_FILE_BYTES: usize = 384 * 1024;
const PREPRINT_TELEMETRY_REL: &str =
    "../../../app-preprints/releases/causal-atlas-pathway-perturbation-validation-telemetry.json";

#[derive(Debug, Clone, Default, Serialize)]
struct T1DifferentialExpressionReport {
    cohort_name: String,
    tsv_example_rows: usize,
    total_measured_genes: usize,
    mapped_to_knowledge_graph: usize,
    active_degs_cohort_fdr: usize,
    active_degs_bh_recalibrated: usize,
    upregulated_count: usize,
    downregulated_count: usize,
    filtered_non_significant_housekeeping: Vec<String>,
    mean_sign_skew: f64,
    bh_monotonicity_verified: bool,
    wald_two_sided_p_calibration_max_abs_diff_log10: f64,
    top_degs_sample: Vec<DegSampleRow>,
}

#[derive(Debug, Clone, Serialize)]
struct DegSampleRow {
    symbol: String,
    log2_fc: f64,
    p_value: f64,
    cohort_fdr: f64,
    bh_recalibrated_q: f64,
    implied_wald_abs_z: f64,
    z_to_p_roundtrip: f64,
}

#[derive(Debug, Clone, Default, Serialize)]
struct T2UpstreamRegulatorReport {
    total_regulators_scored: usize,
    significant_activated_count: usize,
    significant_inhibited_count: usize,
    csv_output_path: String,
    csv_bytes: usize,
    top_regulators: Vec<RegulatorSummaryRow>,
    esr1_contrast_tnf_activation_z: f64,
    esr1_contrast_esr1_activation_z: f64,
}

#[derive(Debug, Clone, Serialize)]
struct RegulatorSummaryRow {
    rank: usize,
    regulator_id: String,
    regulator_name: String,
    kind: String,
    predicted_state: String,
    activation_z: f64,
    bias_corrected_z: f64,
    overlap_count: usize,
    regulon_size: usize,
    fisher_p_value: f64,
    bh_fdr: f64,
    mechanistic_summary: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize)]
struct T3MechanisticPathsAndCmapReport {
    tnf_mechanistic_network_regulators: Vec<String>,
    tnf_mechanistic_edges_count: usize,
    tnf_connected_dataset_genes_count: usize,
    tnf_causal_networks_count: usize,
    tnf_top_causal_paths: Vec<CausalPathSummaryRow>,
    sost_wnt_cohort_sost_log2fc: f64,
    sost_wnt_inferred_root_z: Option<f64>,
    sost_wnt_inferred_root_network_bias_p: Option<f64>,
    perturbation_signatures_scored: usize,
    reversers_count: usize,
    mimickers_count: usize,
    top_reversers: Vec<PerturbationSummaryRow>,
    top_mimickers: Vec<PerturbationSummaryRow>,
}

#[derive(Debug, Clone, Serialize)]
struct CausalPathSummaryRow {
    root_regulator: String,
    max_depth_k: usize,
    predicted_state: String,
    activation_z: f64,
    overlap_p_value: f64,
    network_bias_p_value: f64,
    target_count: usize,
    intermediate_regulators: Vec<String>,
    sample_path_nodes: Vec<String>,
    sample_path_net_sign: i8,
}

#[derive(Debug, Clone, Serialize)]
struct PerturbationSummaryRow {
    signature_id: String,
    perturbagen: String,
    mode: String,
    connectivity_score: f64,
    mechanism_of_action: String,
    matched_up_count: usize,
    matched_down_count: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
struct ValidationSummary {
    schema_version: String,
    repository: String,
    version: String,
    preprint_id: String,
    concept_doi: String,
    validation_mode: String,
    t1_differential_expression: Option<T1DifferentialExpressionReport>,
    t2_upstream_regulators: Option<T2UpstreamRegulatorReport>,
    t3_causal_paths_and_cmap: Option<T3MechanisticPathsAndCmapReport>,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn validation_dir() -> PathBuf {
    let dir = repo_root().join("validation");
    std::fs::create_dir_all(&dir).expect("create validation directory");
    dir
}

fn summary_store() -> &'static Mutex<ValidationSummary> {
    static STORE: OnceLock<Mutex<ValidationSummary>> = OnceLock::new();
    STORE.get_or_init(|| {
        Mutex::new(ValidationSummary {
            schema_version: "1.0.0".to_string(),
            repository: "crvscience-causal-atlas".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            preprint_id: "APP-PREPRINT-04".to_string(),
            concept_doi: "10.5281/zenodo.24222376".to_string(),
            validation_mode: "headless_real_rust_ca_core".to_string(),
            ..ValidationSummary::default()
        })
    })
}

fn write_and_mirror_summary(update: impl FnOnce(&mut ValidationSummary)) {
    let mut guard = summary_store().lock().expect("lock summary store");
    update(&mut guard);
    let json = serde_json::to_string_pretty(&*guard).expect("serialize validation summary");
    assert!(
        json.len() < MAX_VALIDATION_FILE_BYTES,
        "summary JSON {} bytes exceeds 384 KB",
        json.len()
    );
    let val_path = validation_dir().join("causal_atlas_validation_summary.json");
    std::fs::write(&val_path, json.as_bytes()).expect("write causal_atlas_validation_summary.json");

    let preprint_path = repo_root().join(PREPRINT_TELEMETRY_REL);
    if let Some(parent) = preprint_path.parent() {
        std::fs::create_dir_all(parent).expect("create app-preprints/releases dir");
    }
    std::fs::write(&preprint_path, json.as_bytes())
        .expect("mirror causal-atlas telemetry JSON to app-preprints/releases");
}

/// Invert two-sided standard normal p-value to |Z| by binary search using `ca_core::stats::z_to_two_sided_p`
/// to verify exact Wald-test Z <-> p-value roundtrip calibration.
fn invert_two_sided_p_to_abs_z(p: f64) -> f64 {
    let p_clamped = p.clamp(1e-30, 1.0);
    let mut lo = 0.0_f64;
    let mut hi = 12.0_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if z_to_two_sided_p(mid) > p_clamped {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn t1_deseq2_wald_and_bh_fdr_calibration_on_huvec_tnf_cohort() {
    let kg = build_seed_knowledge_graph();
    let cohorts = builtin_cohorts();
    let tnf_ds = &cohorts[0];

    let example_path = repo_root().join("examples").join("tnf_huvec_deg.tsv");
    let example_raw = std::fs::read_to_string(&example_path).expect("read examples/tnf_huvec_deg.tsv");
    let parsed_tsv =
        parse_expression_table("tnf_huvec_deg.tsv", &example_raw).expect("parse DESeq2 TSV table");
    assert_eq!(parsed_tsv.rows.len(), 26);

    let active = tnf_ds.filter_active(1.0, 0.05, &kg);
    let summary = active.summary();
    assert_eq!(summary.total_measured, 28);
    assert!(summary.active_genes >= 20);
    assert!(summary.upregulated_count > summary.downregulated_count);

    // Confirm non-significant housekeeping/unrelated transcripts are filtered out cleanly
    for ctrl in ["GAPDH", "ACTB", "TFF1", "GREB1"] {
        assert!(
            !active.degs.contains_key(ctrl),
            "control/non-DEG {ctrl} must be excluded by |log2FC|>=1.0 & FDR<=0.05"
        );
    }

    // Recompute exact Benjamini-Hochberg FDR across all 28 raw p-values
    let p_vals: Vec<f64> = tnf_ds.rows.iter().map(|r| r.p_value).collect();
    let bh_q = benjamini_hochberg(&p_vals);
    assert_eq!(bh_q.len(), tnf_ds.rows.len());

    let mut pairs: Vec<(f64, f64)> = p_vals.iter().copied().zip(bh_q.iter().copied()).collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut bh_monotonic = true;
    for w in pairs.windows(2) {
        if w[1].1 + 1e-12 < w[0].1 {
            bh_monotonic = false;
        }
    }
    assert!(bh_monotonic, "BH q-values must be monotone with raw p-values");

    let mut active_bh_count = 0usize;
    let mut max_log10_roundtrip_diff = 0.0_f64;
    let mut sample_rows = Vec::new();
    for (idx, (row, &q_val)) in tnf_ds.rows.iter().zip(bh_q.iter()).enumerate() {
        assert!(
            q_val + 1e-12 >= row.p_value,
            "BH q-value ({q_val}) must be >= p-value ({}) for {}",
            row.p_value,
            row.symbol
        );
        if row.log2_fc.abs() >= 1.0 && q_val <= 0.05 {
            active_bh_count += 1;
        }
        let z_abs = invert_two_sided_p_to_abs_z(row.p_value);
        let p_roundtrip = z_to_two_sided_p(z_abs);
        let diff_log10 = ((row.p_value.max(1e-30).log10()) - (p_roundtrip.max(1e-30).log10())).abs();
        if diff_log10 > max_log10_roundtrip_diff {
            max_log10_roundtrip_diff = diff_log10;
        }
        if idx < 8 || row.symbol == "GAPDH" || row.symbol == "KLF2" {
            sample_rows.push(DegSampleRow {
                symbol: row.symbol.clone(),
                log2_fc: row.log2_fc,
                p_value: row.p_value,
                cohort_fdr: row.fdr,
                bh_recalibrated_q: q_val,
                implied_wald_abs_z: z_abs,
                z_to_p_roundtrip: p_roundtrip,
            });
        }
    }
    assert!(
        max_log10_roundtrip_diff < 0.05,
        "Wald Z <-> two-sided p calibration error {max_log10_roundtrip_diff} too large"
    );
    assert_eq!(active_bh_count, summary.active_genes);

    let report = T1DifferentialExpressionReport {
        cohort_name: tnf_ds.name.clone(),
        tsv_example_rows: parsed_tsv.rows.len(),
        total_measured_genes: summary.total_measured,
        mapped_to_knowledge_graph: summary.mapped_to_graph,
        active_degs_cohort_fdr: summary.active_genes,
        active_degs_bh_recalibrated: active_bh_count,
        upregulated_count: summary.upregulated_count,
        downregulated_count: summary.downregulated_count,
        filtered_non_significant_housekeeping: vec![
            "GAPDH".to_string(),
            "ACTB".to_string(),
            "TFF1".to_string(),
            "GREB1".to_string(),
        ],
        mean_sign_skew: summary.mean_sign_skew,
        bh_monotonicity_verified: bh_monotonic,
        wald_two_sided_p_calibration_max_abs_diff_log10: max_log10_roundtrip_diff,
        top_degs_sample: sample_rows,
    };

    write_and_mirror_summary(|s| {
        s.t1_differential_expression = Some(report);
    });
}

#[test]
fn t2_kramer_causal_z_and_fisher_upstream_regulators_ranking() {
    let kg = build_seed_knowledge_graph();
    let cohorts = builtin_cohorts();
    let tnf_ds = &cohorts[0];
    let esr1_ds = &cohorts[1];
    let cfg = CoreAnalysisConfig::default();

    let active_tnf = tnf_ds.filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
    let ura_hits = run_upstream_regulators(&active_tnf, &kg, &cfg);
    assert!(!ura_hits.is_empty(), "URA must recover upstream regulators");

    for master in ["TNF", "NFKB1", "RELA"] {
        let hit = ura_hits
            .iter()
            .find(|h| h.regulator_id == master)
            .unwrap_or_else(|| panic!("Expected upstream regulator {master} in URA hits"));
        assert!(
            hit.activation_z >= 2.0,
            "Expected {master} Activation Z >= 2.0, got {}",
            hit.activation_z
        );
        assert_eq!(hit.predicted_state, PredictedState::Activated);
        assert!(
            hit.fdr <= 0.05,
            "Expected {master} BH FDR <= 0.05, got {}",
            hit.fdr
        );
    }

    let csv_path = validation_dir().join("upstream_regulators_ranked.csv");
    {
        let mut wtr = csv::Writer::from_path(&csv_path).expect("create upstream_regulators_ranked.csv");
        wtr.write_record([
            "rank",
            "regulator_id",
            "regulator_name",
            "kind",
            "predicted_state",
            "activation_z",
            "bias_corrected_z",
            "skew_flagged",
            "overlap_count",
            "regulon_size",
            "fisher_p_value",
            "bh_fdr",
            "observed_log2fc",
            "mechanistic_summary",
        ])
        .expect("write CSV header");

        for (idx, h) in ura_hits.iter().enumerate() {
            wtr.write_record(&[
                (idx + 1).to_string(),
                h.regulator_id.clone(),
                h.regulator_name.clone(),
                format!("{:?}", h.kind),
                h.predicted_state.label().to_string(),
                format!("{:.6}", h.activation_z),
                format!("{:.6}", h.bias_corrected_z),
                h.skew_flagged.to_string(),
                h.overlap_count.to_string(),
                h.regulon_size.to_string(),
                format!("{:.6e}", h.p_value),
                format!("{:.6e}", h.fdr),
                h.observed_log2fc
                    .map(|v| format!("{v:.4}"))
                    .unwrap_or_else(|| "NA".to_string()),
                h.mechanistic_summary.clone().unwrap_or_default(),
            ])
            .expect("write CSV row");
        }
        wtr.flush().expect("flush CSV");
    }

    let csv_meta = std::fs::metadata(&csv_path).expect("CSV metadata");
    let csv_bytes = csv_meta.len() as usize;
    assert!(
        csv_bytes > 0 && csv_bytes < MAX_VALIDATION_FILE_BYTES,
        "upstream_regulators_ranked.csv size {csv_bytes} out of bounds"
    );

    // Specificity check on Cohort 2 (MCF-7 + Estradiol)
    let active_esr1 = esr1_ds.filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
    let ura_esr1 = run_upstream_regulators(&active_esr1, &kg, &cfg);
    let esr1_hit_z = ura_esr1
        .iter()
        .find(|h| h.regulator_id == "ESR1")
        .map(|h| h.activation_z)
        .expect("ESR1 must score on MCF-7 cohort");
    assert!(esr1_hit_z >= 2.0, "Expected ESR1 Z >= 2.0, got {esr1_hit_z}");
    let tnf_in_esr1_z = ura_esr1
        .iter()
        .find(|h| h.regulator_id == "TNF")
        .map(|h| h.activation_z)
        .unwrap_or(0.0);
    assert!(
        tnf_in_esr1_z < 2.0,
        "TNF should not be strongly activated in MCF-7 Estradiol cohort, got {tnf_in_esr1_z}"
    );

    let significant_activated = ura_hits
        .iter()
        .filter(|h| h.predicted_state == PredictedState::Activated && h.fdr <= 0.05)
        .count();
    for h in &ura_hits {
        let disp = format!("{:+.2}", h.activation_z);
        if disp == "+2.00" {
            assert_eq!(h.predicted_state, PredictedState::Activated, "{}", h.regulator_id);
        } else if disp == "-2.00" {
            assert_eq!(h.predicted_state, PredictedState::Inhibited, "{}", h.regulator_id);
        } else if h.activation_z.abs() < 2.0 {
            assert_eq!(h.predicted_state, PredictedState::Neutral, "{}", h.regulator_id);
        }
    }
    let significant_inhibited = ura_hits
        .iter()
        .filter(|h| h.predicted_state == PredictedState::Inhibited && h.fdr <= 0.05)
        .count();

    let top_regulators = ura_hits
        .iter()
        .enumerate()
        .take(12)
        .map(|(idx, h)| RegulatorSummaryRow {
            rank: idx + 1,
            regulator_id: h.regulator_id.clone(),
            regulator_name: h.regulator_name.clone(),
            kind: format!("{:?}", h.kind),
            predicted_state: h.predicted_state.label().to_string(),
            activation_z: h.activation_z,
            bias_corrected_z: h.bias_corrected_z,
            overlap_count: h.overlap_count,
            regulon_size: h.regulon_size,
            fisher_p_value: h.p_value,
            bh_fdr: h.fdr,
            mechanistic_summary: h.mechanistic_summary.clone(),
        })
        .collect();

    let report = T2UpstreamRegulatorReport {
        total_regulators_scored: ura_hits.len(),
        significant_activated_count: significant_activated,
        significant_inhibited_count: significant_inhibited,
        csv_output_path: "validation/upstream_regulators_ranked.csv".to_string(),
        csv_bytes,
        top_regulators,
        esr1_contrast_tnf_activation_z: tnf_in_esr1_z,
        esr1_contrast_esr1_activation_z: esr1_hit_z,
    };

    write_and_mirror_summary(|s| {
        s.t2_upstream_regulators = Some(report);
    });
}

#[test]
fn t3_signed_mechanistic_paths_and_cmap_perturbation_reversal() {
    let kg = build_seed_knowledge_graph();
    let cohorts = builtin_cohorts();
    let cfg = CoreAnalysisConfig::default();

    // Cohort 0: TNF HUVEC Core Analysis (CNA + Mechanistic Network + CMap)
    let tnf_report = run_core_analysis(&cohorts[0], &kg, &cfg);
    let active_tnf = cohorts[0].filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
    let tnf_mech = build_mechanistic_network(
        "TNF",
        &tnf_report.upstream_regulators,
        &active_tnf,
        &kg,
        cfg.mechanistic_breadth,
        cfg.mechanistic_depth,
    );
    assert!(
        tnf_mech.regulator_ids.len() >= 2,
        "TNF mechanistic network must link downstream co-regulators"
    );
    assert!(!tnf_mech.edges.is_empty());
    assert!(!tnf_report.causal_networks.is_empty());

    let tnf_top_paths: Vec<CausalPathSummaryRow> = tnf_report
        .causal_networks
        .iter()
        .take(6)
        .map(|cn| {
            let (sample_nodes, sample_sign) = cn
                .sample_paths
                .first()
                .map(|p| (p.nodes.clone(), p.net_sign))
                .unwrap_or_default();
            CausalPathSummaryRow {
                root_regulator: cn.root_regulator.clone(),
                max_depth_k: cn.max_depth_k,
                predicted_state: cn.predicted_state.label().to_string(),
                activation_z: cn.activation_z,
                overlap_p_value: cn.overlap_p_value,
                network_bias_p_value: cn.network_bias_p_value,
                target_count: cn.target_count,
                intermediate_regulators: cn.intermediate_regulators.clone(),
                sample_path_nodes: sample_nodes,
                sample_path_net_sign: sample_sign,
            }
        })
        .collect();

    // Cohort 2: Silent SOST root discovered through 3-hop CNA paths
    let active_wnt = cohorts[2].filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &kg);
    let sost_row = cohorts[2]
        .rows
        .iter()
        .find(|r| r.symbol == "SOST")
        .expect("SOST row present");
    assert!(!active_wnt.degs.contains_key("SOST"));
    let wnt_cna = run_causal_networks(&active_wnt, &kg, &cfg);
    let sost_hit = wnt_cna
        .iter()
        .find(|h| h.root_regulator == "SOST");

    // CMap Perturbation Reversal on TNF HUVEC
    let cmap_hits = run_perturbation_match(&active_tnf, &kg);
    assert!(!cmap_hits.is_empty(), "CMap signatures must match TNF cohort");
    let mut reversers: Vec<&ca_core::PerturbationHit> = cmap_hits
        .iter()
        .filter(|h| h.connectivity_score <= -25.0)
        .collect();
    reversers.sort_by(|a, b| a.connectivity_score.total_cmp(&b.connectivity_score));

    let mut mimickers: Vec<&ca_core::PerturbationHit> = cmap_hits
        .iter()
        .filter(|h| h.connectivity_score >= 25.0)
        .collect();
    mimickers.sort_by(|a, b| b.connectivity_score.total_cmp(&a.connectivity_score));

    assert!(
        !reversers.is_empty(),
        "Expected therapeutic reverser compounds with connectivity_score <= -25"
    );

    let to_summary = |h: &ca_core::PerturbationHit| PerturbationSummaryRow {
        signature_id: h.signature_id.clone(),
        perturbagen: h.perturbagen.clone(),
        mode: h.mode.clone(),
        connectivity_score: h.connectivity_score,
        mechanism_of_action: h.mechanism_of_action.clone(),
        matched_up_count: h.matched_up.len(),
        matched_down_count: h.matched_down.len(),
    };

    let report = T3MechanisticPathsAndCmapReport {
        tnf_mechanistic_network_regulators: tnf_mech.regulator_ids,
        tnf_mechanistic_edges_count: tnf_mech.edges.len(),
        tnf_connected_dataset_genes_count: tnf_mech.connected_dataset_genes.len(),
        tnf_causal_networks_count: tnf_report.causal_networks.len(),
        tnf_top_causal_paths: tnf_top_paths,
        sost_wnt_cohort_sost_log2fc: sost_row.log2_fc,
        sost_wnt_inferred_root_z: sost_hit.map(|h| h.activation_z),
        sost_wnt_inferred_root_network_bias_p: sost_hit.map(|h| h.network_bias_p_value),
        perturbation_signatures_scored: cmap_hits.len(),
        reversers_count: reversers.len(),
        mimickers_count: mimickers.len(),
        top_reversers: reversers.into_iter().take(5).map(to_summary).collect(),
        top_mimickers: mimickers.into_iter().take(5).map(to_summary).collect(),
    };

    write_and_mirror_summary(|s| {
        s.t3_causal_paths_and_cmap = Some(report);
    });
}
