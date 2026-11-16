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

use crate::graph::KnowledgeGraph;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum DatasetError {
    #[error("Dataset table is empty or missing header row")]
    EmptyInput,
    #[error("Could not identify gene symbol and log2 fold-change columns in header: {0}")]
    MissingRequiredColumns(String),
    #[error("No valid numeric gene measurements parsed from input")]
    NoValidRows,
}

/// Single gene/protein differential measurement from RNA-seq, microarray, or proteomics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneMeasurement {
    pub symbol: String,
    pub log2_fc: f64,
    pub p_value: f64,
    pub fdr: f64,
}

/// Full uploaded or built-in differential expression cohort before threshold filtering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpressionDataset {
    pub name: String,
    pub description: String,
    pub rows: Vec<GeneMeasurement>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetSummary {
    pub total_measured: usize,
    pub mapped_to_graph: usize,
    pub active_genes: usize,
    pub upregulated_count: usize,
    pub downregulated_count: usize,
    pub mean_sign_skew: f64,
}

/// Filtered differentially expressed gene subset `D` mapped against the KnowledgeGraph.
#[derive(Debug, Clone)]
pub struct ActiveDataset {
    pub name: String,
    pub total_measured: usize,
    pub mapped_to_graph: usize,
    /// Map of uppercase HGNC symbol -> GeneMeasurement for active DEGs passing cutoffs.
    pub degs: BTreeMap<String, GeneMeasurement>,
    /// Full ranked list `(symbol, log2_fc)` across all measured genes for GSEA / Perturbation connectivity.
    pub all_ranked: Vec<(String, f64)>,
}

impl ExpressionDataset {
    pub fn filter_active(
        &self,
        abs_log2fc_cutoff: f64,
        fdr_cutoff: f64,
        kg: &KnowledgeGraph,
    ) -> ActiveDataset {
        let mut mapped_to_graph = 0usize;
        let mut degs = BTreeMap::new();
        let mut ranked_map: BTreeMap<String, f64> = BTreeMap::new();

        for row in &self.rows {
            let sym = row.symbol.trim().to_uppercase();
            if sym.is_empty() || !row.log2_fc.is_finite() {
                continue;
            }
            let in_kg = kg.nodes.contains_key(&sym)
                || kg.regulated_gene_universe.contains(&sym)
                || kg.pathway_gene_universe.contains(&sym);
            if in_kg {
                mapped_to_graph += 1;
            }
            // Retain strongest magnitude measurement if duplicate transcripts map to the same gene symbol.
            ranked_map
                .entry(sym.clone())
                .and_modify(|prev| {
                    if row.log2_fc.abs() > prev.abs() {
                        *prev = row.log2_fc;
                    }
                })
                .or_insert(row.log2_fc);

            if row.log2_fc.abs() >= abs_log2fc_cutoff && row.fdr <= fdr_cutoff {
                degs.entry(sym)
                    .and_modify(|prev: &mut GeneMeasurement| {
                        if row.log2_fc.abs() > prev.log2_fc.abs() {
                            *prev = row.clone();
                        }
                    })
                    .or_insert_with(|| GeneMeasurement {
                        symbol: row.symbol.trim().to_uppercase(),
                        log2_fc: row.log2_fc,
                        p_value: row.p_value,
                        fdr: row.fdr,
                    });
            }
        }

        let mut all_ranked: Vec<(String, f64)> = ranked_map.into_iter().collect();
        all_ranked.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.0.cmp(&b.0))
        });

        ActiveDataset {
            name: self.name.clone(),
            total_measured: self.rows.len(),
            mapped_to_graph,
            degs,
            all_ranked,
        }
    }
}

impl ActiveDataset {
    pub fn sign_of(&self, symbol: &str) -> Option<i8> {
        self.degs.get(symbol).map(|m| {
            if m.log2_fc > 0.0 {
                1
            } else if m.log2_fc < 0.0 {
                -1
            } else {
                0
            }
        })
    }

    pub fn summary(&self) -> DatasetSummary {
        let mut up = 0usize;
        let mut down = 0usize;
        for m in self.degs.values() {
            if m.log2_fc > 0.0 {
                up += 1;
            } else if m.log2_fc < 0.0 {
                down += 1;
            }
        }
        let active_genes = self.degs.len();
        let mean_sign_skew = if active_genes > 0 {
            (up as f64 - down as f64) / (active_genes as f64)
        } else {
            0.0
        };
        DatasetSummary {
            total_measured: self.total_measured,
            mapped_to_graph: self.mapped_to_graph,
            active_genes,
            upregulated_count: up,
            downregulated_count: down,
            mean_sign_skew,
        }
    }
}

/// Auto-detects CSV or TSV differential expression exports (DESeq2, edgeR, limma, Scanpy).
pub fn parse_expression_table(
    name: &str,
    raw_text: &str,
) -> Result<ExpressionDataset, DatasetError> {
    let mut lines = raw_text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));

    let header_line = lines.next().ok_or(DatasetError::EmptyInput)?;
    let delimiter = if header_line.contains('\t') {
        b'\t'
    } else if header_line.contains(';') {
        b';'
    } else {
        b','
    };

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .flexible(true)
        .has_headers(true)
        .from_reader(raw_text.as_bytes());

    let headers = rdr
        .headers()
        .map_err(|_| DatasetError::EmptyInput)?
        .iter()
        .map(|h| h.trim().trim_matches('"').to_lowercase())
        .collect::<Vec<_>>();

    let gene_col = headers.iter().position(|h| {
        matches!(
            h.as_str(),
            "gene"
                | "symbol"
                | "gene_symbol"
                | "genesymbol"
                | "hgnc_symbol"
                | "id"
                | "gene_name"
                | "names"
                | ""
        )
    });

    let fc_col = headers.iter().position(|h| {
        matches!(
            h.as_str(),
            "log2foldchange"
                | "log2_fc"
                | "log2fc"
                | "logfc"
                | "logfoldchanges"
                | "avg_log2fc"
                | "lfc"
                | "fc"
        )
    });

    let pval_col = headers.iter().position(|h| {
        matches!(
            h.as_str(),
            "pvalue" | "p_value" | "pval" | "p.value" | "p_val"
        )
    });

    let fdr_col = headers.iter().position(|h| {
        matches!(
            h.as_str(),
            "padj" | "fdr" | "adj.p.val" | "qvalue" | "q_value" | "p_val_adj"
        )
    });

    let (Some(g_idx), Some(fc_idx)) = (gene_col, fc_col) else {
        return Err(DatasetError::MissingRequiredColumns(header_line.to_string()));
    };

    let mut rows = Vec::new();
    for rec in rdr.records().flatten() {
        let Some(sym_raw) = rec.get(g_idx) else {
            continue;
        };
        let symbol = sym_raw.trim().trim_matches('"').to_uppercase();
        if symbol.is_empty() || symbol == "NA" {
            continue;
        }
        let Some(fc_str) = rec.get(fc_idx) else {
            continue;
        };
        let Ok(log2_fc) = fc_str.trim().trim_matches('"').parse::<f64>() else {
            continue;
        };
        if !log2_fc.is_finite() {
            continue;
        }
        let p_value = pval_col
            .and_then(|i| rec.get(i))
            .and_then(|s| s.trim().trim_matches('"').parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v >= 0.0)
            .unwrap_or(0.001);
        let fdr = fdr_col
            .and_then(|i| rec.get(i))
            .and_then(|s| s.trim().trim_matches('"').parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v >= 0.0)
            .unwrap_or(p_value);

        rows.push(GeneMeasurement {
            symbol,
            log2_fc,
            p_value,
            fdr,
        });
    }

    if rows.is_empty() {
        return Err(DatasetError::NoValidRows);
    }

    Ok(ExpressionDataset {
        name: name.to_string(),
        description: format!("Imported differential expression dataset ({} rows)", rows.len()),
        rows,
    })
}

/// Three built-in benchmark expression profiles modeled on canonical perturbation paradigms
/// (TNF endothelial inflammation, ESR1/Estradiol hormone response, and SOST-inhibited Wnt osteoproliferation).
pub fn builtin_cohorts() -> Vec<ExpressionDataset> {
    let tnf_rows = vec![
        ("ICAM1", 4.35, 1.2e-18, 8.0e-17),
        ("VCAM1", 4.82, 3.5e-20, 4.1e-18),
        ("SELE", 5.10, 1.1e-21, 2.0e-19),
        ("CCL2", 3.95, 2.4e-15, 9.5e-14),
        ("CXCL8", 4.60, 4.0e-19, 3.2e-17),
        ("CXCL10", 3.74, 8.2e-14, 2.1e-12),
        ("IL6", 3.88, 1.9e-14, 6.0e-13),
        ("IL1B", 3.12, 6.5e-11, 1.1e-9),
        ("PTGS2", 2.94, 2.3e-10, 3.4e-9),
        ("NFKBIA", 2.75, 1.4e-9, 1.8e-8),
        ("TNFAIP3", 3.28, 3.1e-12, 6.2e-11),
        ("BIRC3", 2.62, 5.5e-9, 6.1e-8),
        ("SOD2", 2.48, 1.9e-8, 1.9e-7),
        ("MMP9", 2.85, 8.0e-10, 1.1e-8),
        ("CSF2", 3.05, 9.4e-11, 1.5e-9),
        ("IRF1", 2.42, 3.1e-8, 2.9e-7),
        ("STAT1", 1.92, 4.2e-6, 2.8e-5),
        ("NOS3", -2.14, 8.5e-7, 6.4e-6),
        ("KLF2", -2.38, 1.1e-7, 9.5e-7),
        ("THBD", -1.85, 9.0e-6, 5.2e-5),
        ("BMPR2", -1.42, 3.1e-4, 1.4e-3),
        ("PPARG", -1.68, 4.8e-5, 2.6e-4),
        ("CCND1", 1.35, 6.0e-4, 2.5e-3),
        ("MYC", 1.28, 1.1e-3, 4.2e-3),
        ("GAPDH", 0.08, 0.74, 0.89),
        ("ACTB", -0.05, 0.81, 0.92),
        ("TFF1", -0.22, 0.42, 0.68),
        ("GREB1", 0.14, 0.58, 0.79),
    ];

    let esr1_rows = vec![
        ("TFF1", 4.45, 2.0e-19, 1.9e-17),
        ("GREB1", 4.18, 8.5e-18, 6.2e-16),
        ("PGR", 3.76, 4.1e-15, 1.8e-13),
        ("CXCL12", 3.22, 7.0e-12, 1.9e-10),
        ("MYC", 2.84, 9.2e-10, 1.4e-8),
        ("CCND1", 2.68, 3.8e-9, 4.9e-8),
        ("AREG", 3.15, 1.3e-11, 3.1e-10),
        ("BCL2", 2.34, 4.6e-8, 4.5e-7),
        ("CTSD", 2.21, 1.2e-7, 1.1e-6),
        ("XBP1", 2.05, 5.1e-7, 4.0e-6),
        ("E2F1", 1.96, 1.1e-6, 8.2e-6),
        ("PCNA", 1.82, 4.4e-6, 2.9e-5),
        ("MKI67", 1.79, 6.2e-6, 3.8e-5),
        ("CDKN1A", -2.28, 7.5e-8, 7.1e-7),
        ("CDKN1B", -1.94, 1.5e-6, 1.0e-5),
        ("FOXO3", -1.76, 9.8e-6, 5.5e-5),
        ("BAX", -1.52, 1.2e-4, 5.8e-4),
        ("ERBB2", -1.48, 1.9e-4, 8.7e-4),
        ("ICAM1", -1.18, 2.1e-3, 7.8e-3),
        ("IL6", -1.32, 8.5e-4, 3.4e-3),
        ("CXCL8", -1.15, 2.8e-3, 9.9e-3),
        ("GAPDH", 0.04, 0.88, 0.95),
    ];

    let wnt_sost_rows = vec![
        ("RUNX2", 2.92, 3.2e-10, 4.5e-9),
        ("SP7", 3.18, 4.8e-11, 8.2e-10),
        ("ALPL", 2.74, 1.8e-9, 2.1e-8),
        ("BGLAP", 3.05, 1.1e-10, 1.6e-9),
        ("COL1A1", 2.46, 2.4e-8, 2.3e-7),
        ("BMP2", 2.12, 9.0e-7, 6.8e-6),
        ("BMP6", 1.95, 3.6e-6, 2.4e-5),
        ("AXIN2", 2.88, 5.1e-10, 6.7e-9),
        ("LEF1", 2.31, 6.2e-8, 5.4e-7),
        ("MYC", 1.88, 6.8e-6, 4.1e-5),
        ("CCND1", 2.04, 1.5e-6, 1.1e-5),
        ("LCN2", 2.58, 7.9e-9, 8.3e-8),
        ("MMP13", 2.41, 3.8e-8, 3.5e-7),
        ("PTGS2", 1.65, 6.4e-5, 3.3e-4),
        ("PPARG", -2.18, 5.4e-7, 4.3e-6),
        ("CEBPA", -1.82, 1.2e-5, 6.5e-5),
        ("CDKN1A", -1.58, 9.5e-5, 4.6e-4),
        // Note: SOST itself is intentionally below significance cutoff (unmeasured / unchanged mRNA)
        // to demonstrate multi-hop Causal Network Analysis root discovery per Krämer et al. (2014).
        ("SOST", -0.35, 0.24, 0.48),
        ("DKK1", -0.41, 0.19, 0.41),
    ];

    let to_ds = |name: &str, desc: &str, items: Vec<(&str, f64, f64, f64)>| ExpressionDataset {
        name: name.to_string(),
        description: desc.to_string(),
        rows: items
            .into_iter()
            .map(|(s, fc, p, q)| GeneMeasurement {
                symbol: s.to_string(),
                log2_fc: fc,
                p_value: p,
                fdr: q,
            })
            .collect(),
    };

    vec![
        to_ds(
            "HUVEC + TNF Inflammatory Cascade",
            "Endothelial acute cytokine response benchmark (NF-kB / STAT1 activation & vascular adhesion program)",
            tnf_rows,
        ),
        to_ds(
            "MCF-7 + Beta-Estradiol Signaling",
            "Estrogen receptor alpha (ESR1) mitogenic transcriptional program & selective anti-estrogen reversal target",
            esr1_rows,
        ),
        to_ds(
            "Spinal Osteoproliferation (Multi-Hop Wnt/SOST)",
            "Osteogenic differentiation profile where root inhibitor SOST is transcriptionally silent but inferred via 3-hop CNA",
            wnt_sost_rows,
        ),
    ]
}
