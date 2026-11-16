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

//! # `ca-streams` — SQLite WAL Cache, Batch Pathway Downloaders & Open-Source Stream Connectors
//!
//! Combines a bundled SQLite 3 (`WAL` mode) persistent knowledge & analysis cache with
//! pure-Rust HTTPS (`ureq` + `rustls`), pure-Rust `.zip` / `.gz` archive extraction, ETag /
//! SHA-256 conditional incremental updates, and parsers for 8 open-source biological formats.

pub mod batch_downloader;
pub mod client;
pub mod parsers;
pub mod sqlite_store;

use ca_core::KnowledgeGraph;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::path::Path;
use thiserror::Error;

pub use batch_downloader::{
    fetch_and_cache_reactome_diagram, run_batch_download_and_update, sync_reactome_gmt_bundle,
    sync_wikipathways_gpml_bundle, BatchProgressEvent, BatchSyncSummary,
};
pub use sqlite_store::{SqliteCacheStats, SqliteKnowledgeStore, StreamCacheMeta};

#[derive(Debug, Error)]
pub enum StreamError {
    #[error("Network stream error: {0}")]
    Network(String),
    #[error("SQLite database error: {0}")]
    Database(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Format parse error: {0}")]
    Parse(String),
    #[error("Serialization error: {0}")]
    Json(#[from] serde_json::Error),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StreamFormat {
    GmtPathways,
    GpmlXml,
    SignorTsv,
    OmniPathTsv,
    StringTsv,
    GoGafTsv,
    OpenTargetsGraphQlJson,
    LincsSignaturesJson,
}

impl StreamFormat {
    pub fn label(self) -> &'static str {
        match self {
            Self::GmtPathways => "GMT Gene Matrix Transposed",
            Self::GpmlXml => "WikiPathways GPML XML",
            Self::SignorTsv => "SIGNOR Signed Causal TSV",
            Self::OmniPathTsv => "OmniPath / CollecTRI TSV",
            Self::StringTsv => "STRING Regulatory TSV",
            Self::GoGafTsv => "Gene Ontology GAF 2.1",
            Self::OpenTargetsGraphQlJson => "Open Targets GraphQL JSON",
            Self::LincsSignaturesJson => "SigCom LINCS Consensus JSON",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StreamDescriptor {
    pub id: String,
    pub name: String,
    pub role_in_analysis: String,
    pub license: String,
    pub homepage_url: String,
    pub endpoint_url: String,
    pub format: StreamFormat,
    pub enabled_by_default: bool,
    pub last_sync_status: String,
    pub last_sha256: Option<String>,
    pub records_ingested: usize,
}

/// Returns the catalog of open-source biological knowledge streams cached in SQLite.
pub fn default_stream_registry() -> Vec<StreamDescriptor> {
    vec![
        StreamDescriptor {
            id: "omnipath_collectri".to_string(),
            name: "OmniPath + CollecTRI Signed TF Regulons".to_string(),
            role_in_analysis: "Upstream Regulator Analysis (E_t transcription factor-target activation/repression)".to_string(),
            license: "CC-BY 4.0 / MIT (Saez-Rodriguez Lab)".to_string(),
            homepage_url: "https://github.com/saezlab/collecTRI".to_string(),
            endpoint_url: "https://omnipathdb.org/interactions?datasets=collectri&fields=sources,references,curation_effort&genesymbols=1&organisms=9606".to_string(),
            format: StreamFormat::OmniPathTsv,
            enabled_by_default: true,
            last_sync_status: "Seeded in SQLite (ready for ETag/SHA-256 sync)".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "omnipath_signaling".to_string(),
            name: "OmniPath Literature-Curated Causal Signaling".to_string(),
            role_in_analysis: "Mechanistic & Multi-Hop Causal Networks (E_a phosphorylation & receptor cascades)".to_string(),
            license: "Open Academic / Constituent CC-BY".to_string(),
            homepage_url: "https://omnipathdb.org/".to_string(),
            endpoint_url: "https://omnipathdb.org/interactions?datasets=omnipath&fields=sources,references,curation_effort&genesymbols=1&organisms=9606".to_string(),
            format: StreamFormat::OmniPathTsv,
            enabled_by_default: true,
            last_sync_status: "Seeded in SQLite (ready for ETag/SHA-256 sync)".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "signor_causal".to_string(),
            name: "SIGNOR 3.0 Causal Network & 114 Signaling Pathway Maps".to_string(),
            role_in_analysis: "Signed E_a / E_t / E_p cause-effect interactions with biochemical mechanisms & PMIDs".to_string(),
            license: "CC-BY 4.0 (University of Rome Tor Vergata)".to_string(),
            homepage_url: "https://signor.uniroma2.it/downloads.php".to_string(),
            endpoint_url: "https://signor.uniroma2.it/getData.php?organism=9606".to_string(),
            format: StreamFormat::SignorTsv,
            enabled_by_default: true,
            last_sync_status: "Seeded in SQLite (ready for ETag/SHA-256 sync)".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "wikipathways_gmt".to_string(),
            name: "WikiPathways Homo sapiens (866+ Curated GPML/GMT Pathways)".to_string(),
            role_in_analysis: "Canonical Pathway ORA & GPML Coordinate Diagrams".to_string(),
            license: "CC0 1.0 Public Domain Dedication".to_string(),
            homepage_url: "https://data.wikipathways.org/current/gpml/".to_string(),
            endpoint_url: "https://data.wikipathways.org/current/gpml/".to_string(),
            format: StreamFormat::GpmlXml,
            enabled_by_default: true,
            last_sync_status: "Use Batch Downloader to cache all 866+ GPML diagrams".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "reactome_content".to_string(),
            name: "Reactome Pathway Knowledgebase (2,647+ Human Pathways)".to_string(),
            role_in_analysis: "Canonical Pathways GMT Zip & On-Demand DiagramJS Layouts".to_string(),
            license: "CC-BY 4.0 / CC0".to_string(),
            homepage_url: "https://reactome.org/ContentService".to_string(),
            endpoint_url: "https://reactome.org/download/current/ReactomePathways.gmt.zip".to_string(),
            format: StreamFormat::GmtPathways,
            enabled_by_default: false,
            last_sync_status: "Use Batch Downloader to cache all 2,647+ Reactome pathways".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "pathway_commons_all".to_string(),
            name: "Pathway Commons v12 Integrated Compendium (4,794 Human Pathways)".to_string(),
            role_in_analysis: "Unified BioPAX / SBGN / GMT compendium across 22 open pathway databases".to_string(),
            license: "Open Academic / Constituent CC-BY & CC0".to_string(),
            homepage_url: "https://www.pathwaycommons.org/".to_string(),
            endpoint_url: "https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.gmt.gz".to_string(),
            format: StreamFormat::GmtPathways,
            enabled_by_default: false,
            last_sync_status: "Ready for HTTPS gzip download into SQLite".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "string_regulatory".to_string(),
            name: "STRING v12 Directed Regulatory & Physical Network".to_string(),
            role_in_analysis: "Neighborhood Context & Protein Complex Expansion".to_string(),
            license: "CC-BY 4.0 (STRING Consortium)".to_string(),
            homepage_url: "https://string-db.org/help/api/".to_string(),
            endpoint_url: "https://string-db.org/api/tsv/network?identifiers=TNF%0dNFKB1%0dESR1%0dCTNNB1%0dRUNX2%0dTP53&species=9606".to_string(),
            format: StreamFormat::StringTsv,
            enabled_by_default: true,
            last_sync_status: "Ready for live HTTPS sync".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "gene_ontology_gaf".to_string(),
            name: "Gene Ontology (GO) & GOA Human Annotations".to_string(),
            role_in_analysis: "Downstream Effects Analysis (Biological Process TreeMap)".to_string(),
            license: "CC-BY 4.0 (Gene Ontology Consortium)".to_string(),
            homepage_url: "https://geneontology.org/docs/download-go-annotations/".to_string(),
            endpoint_url: "https://current.geneontology.org/annotations/goa_human.gaf.gz".to_string(),
            format: StreamFormat::GoGafTsv,
            enabled_by_default: false,
            last_sync_status: "Seeded in SQLite (import local .gaf.gz or sync)".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "open_targets_graphql".to_string(),
            name: "Open Targets Platform v4 Disease & Drug GraphQL".to_string(),
            role_in_analysis: "Disease Phenotype Associations & Clinical Drug Target Modulation".to_string(),
            license: "CC0 1.0 / Apache-2.0 (Open Targets)".to_string(),
            homepage_url: "https://platform-docs.opentargets.org/data-access/graphql-api".to_string(),
            endpoint_url: "https://api.platform.opentargets.org/api/v4/graphql".to_string(),
            format: StreamFormat::OpenTargetsGraphQlJson,
            enabled_by_default: true,
            last_sync_status: "Seeded in SQLite (ready for live GraphQL sync)".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
        StreamDescriptor {
            id: "sigcom_lincs".to_string(),
            name: "NIH LINCS L1000 / SigCom Perturbation Signatures".to_string(),
            role_in_analysis: "Perturbation Match (Therapeutic Reverser & Mimicker Connectivity)".to_string(),
            license: "Open Access / Apache-2.0 (Ma'ayan Lab / NIH LINCS)".to_string(),
            homepage_url: "https://lincsproject.org/".to_string(),
            endpoint_url: "https://maayanlab.cloud/sigcom-lincs".to_string(),
            format: StreamFormat::LincsSignaturesJson,
            enabled_by_default: false,
            last_sync_status: "Consensus signatures cached in SQLite".to_string(),
            last_sha256: None,
            records_ingested: 0,
        },
    ]
}

/// Portable export helper: writes a `KnowledgeGraph` snapshot to a gzip-compressed JSON file (`.ca-kb.json.gz`).
pub fn save_compressed_snapshot(kg: &KnowledgeGraph, path: &Path) -> Result<(), StreamError> {
    let json_bytes = serde_json::to_vec(kg)?;
    let file = std::fs::File::create(path)?;
    let mut gz = GzEncoder::new(file, Compression::fast());
    gz.write_all(&json_bytes)?;
    gz.finish()?;
    Ok(())
}

/// Portable import helper: loads a `KnowledgeGraph` snapshot from `.ca-kb.json.gz` or `.json`.
pub fn load_compressed_snapshot(path: &Path) -> Result<KnowledgeGraph, StreamError> {
    let raw = std::fs::read(path)?;
    let decoded = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut gz = GzDecoder::new(&raw[..]);
        let mut out = Vec::new();
        gz.read_to_end(&mut out)?;
        out
    } else {
        raw
    };
    let mut kg: KnowledgeGraph = serde_json::from_slice(&decoded)?;
    kg.rebuild_indices();
    Ok(kg)
}
