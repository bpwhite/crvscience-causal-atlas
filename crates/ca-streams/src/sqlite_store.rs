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

use crate::StreamError;
use ca_core::{
    CausalEdge, CoreAnalysisConfig, CoreAnalysisReport, EdgeClass, EdgeSign, EntityKind,
    EntityNode, ExpressionDataset, KnowledgeGraph, PathwayDef, PerturbationSignature,
    SubcellularLocation,
};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rusqlite::{params, Connection, OptionalExtension};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

/// Summary telemetry for the local SQLite WAL database & result cache.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SqliteCacheStats {
    pub db_path: String,
    pub file_size_bytes: u64,
    pub wal_mode: bool,
    pub cached_entities: usize,
    pub cached_edges: usize,
    pub cached_pathways: usize,
    pub cached_diagrams_with_coords: usize,
    pub cached_signatures: usize,
    pub synced_streams: usize,
    pub cached_analysis_reports: usize,
    pub kg_revision_hash: String,
}

/// HTTP conditional & TTL metadata stored per open-source data stream.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StreamCacheMeta {
    pub stream_id: String,
    pub etag: Option<String>,
    pub last_modified: Option<String>,
    pub content_sha256: String,
    pub records_count: usize,
    pub synced_at_unix: u64,
    pub ttl_seconds: u64,
}

impl StreamCacheMeta {
    pub fn is_fresh(&self, now_unix: u64) -> bool {
        now_unix.saturating_sub(self.synced_at_unix) < self.ttl_seconds
    }
}

/// Persistent SQLite 3 database (`WAL` journal mode) caching all open-source biological streams,
/// GPML / Reactome diagram coordinates, HTTP ETag/SHA-256 sync state, and Core Analysis reports.
pub struct SqliteKnowledgeStore {
    conn: Connection,
    path: PathBuf,
}

pub fn unix_now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

fn gzip_compress(bytes: &[u8]) -> Result<Vec<u8>, StreamError> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::fast());
    enc.write_all(bytes)?;
    Ok(enc.finish()?)
}

fn gzip_decompress(bytes: &[u8]) -> Result<Vec<u8>, StreamError> {
    let mut dec = GzDecoder::new(bytes);
    let mut out = Vec::new();
    dec.read_to_end(&mut out)?;
    Ok(out)
}

fn hex_sha256(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    let d = h.finalize();
    let mut s = String::with_capacity(64);
    for b in d {
        use std::fmt::Write;
        let _ = write!(&mut s, "{b:02x}");
    }
    s
}

impl SqliteKnowledgeStore {
    /// Opens or creates the SQLite knowledge & analysis cache at `path`, enabling WAL mode.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, StreamError> {
        let p = path.as_ref().to_path_buf();
        if let Some(parent) = p.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let conn = Connection::open(&p).map_err(|e| StreamError::Database(e.to_string()))?;
        let store = Self { conn, path: p };
        store.init_schema()?;
        Ok(store)
    }

    /// Opens an ephemeral in-memory SQLite store (used by unit tests and stateless CLI runs).
    pub fn open_in_memory() -> Result<Self, StreamError> {
        let conn = Connection::open_in_memory().map_err(|e| StreamError::Database(e.to_string()))?;
        let store = Self {
            conn,
            path: PathBuf::from(":memory:"),
        };
        store.init_schema()?;
        Ok(store)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    fn init_schema(&self) -> Result<(), StreamError> {
        self.conn
            .execute_batch(
                "PRAGMA journal_mode = WAL;
                 PRAGMA synchronous = NORMAL;
                 PRAGMA temp_store = MEMORY;
                 PRAGMA foreign_keys = ON;

                 CREATE TABLE IF NOT EXISTS kv_meta (
                     key TEXT PRIMARY KEY,
                     value TEXT NOT NULL
                 );

                 CREATE TABLE IF NOT EXISTS stream_sync_meta (
                     stream_id TEXT PRIMARY KEY,
                     etag TEXT,
                     last_modified TEXT,
                     content_sha256 TEXT NOT NULL,
                     records_count INTEGER NOT NULL,
                     synced_at_unix INTEGER NOT NULL,
                     ttl_seconds INTEGER NOT NULL DEFAULT 604800
                 );

                 CREATE TABLE IF NOT EXISTS entities (
                     id TEXT PRIMARY KEY,
                     display_name TEXT NOT NULL,
                     kind_json TEXT NOT NULL,
                     location_json TEXT NOT NULL,
                     family_category TEXT NOT NULL,
                     external_id TEXT,
                     provenance_source TEXT NOT NULL
                 );

                 CREATE TABLE IF NOT EXISTS causal_edges (
                     source TEXT NOT NULL,
                     target TEXT NOT NULL,
                     edge_class_json TEXT NOT NULL,
                     sign_i8 INTEGER NOT NULL,
                     weight REAL NOT NULL,
                     mechanism TEXT NOT NULL,
                     source_db TEXT NOT NULL,
                     pmids_json TEXT NOT NULL,
                     PRIMARY KEY (source, target, edge_class_json)
                 );
                 CREATE INDEX IF NOT EXISTS idx_edges_source ON causal_edges(source);
                 CREATE INDEX IF NOT EXISTS idx_edges_target ON causal_edges(target);

                 CREATE TABLE IF NOT EXISTS canonical_pathways (
                     id TEXT PRIMARY KEY,
                     name TEXT NOT NULL,
                     source_db TEXT NOT NULL,
                     category TEXT NOT NULL,
                     url TEXT NOT NULL,
                     members_json TEXT NOT NULL,
                     expected_signs_json TEXT NOT NULL,
                     topology_json TEXT NOT NULL,
                     coords_json TEXT NOT NULL,
                     phospho_json TEXT NOT NULL,
                     cell_context_label TEXT NOT NULL,
                     has_coords INTEGER NOT NULL DEFAULT 0,
                     updated_at_unix INTEGER NOT NULL
                 );
                 CREATE INDEX IF NOT EXISTS idx_pathways_source ON canonical_pathways(source_db);

                 CREATE TABLE IF NOT EXISTS perturbation_signatures (
                     id TEXT PRIMARY KEY,
                     signature_json TEXT NOT NULL
                 );

                 CREATE TABLE IF NOT EXISTS analysis_result_cache (
                     cache_key TEXT PRIMARY KEY,
                     dataset_name TEXT NOT NULL,
                     computed_at_unix INTEGER NOT NULL,
                     report_gzip BLOB NOT NULL
                 );",
            )
            .map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(())
    }

    /// Seeds the SQLite database with the starter open-source knowledge graph if currently empty,
    /// and backfills expanded seed perturbation signatures into pre-existing databases.
    pub fn ensure_seeded(&mut self, seed_fn: fn() -> KnowledgeGraph) -> Result<bool, StreamError> {
        let count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM canonical_pathways", [], |r| r.get(0))
            .map_err(|e| StreamError::Database(e.to_string()))?;
        if count == 0 {
            let seed_kg = seed_fn();
            self.upsert_knowledge_graph(&seed_kg)?;
            return Ok(true);
        }

        let sig_count: i64 = self
            .conn
            .query_row("SELECT COUNT(*) FROM perturbation_signatures", [], |r| {
                r.get(0)
            })
            .map_err(|e| StreamError::Database(e.to_string()))?;
        if sig_count < 25 {
            let seed_kg = seed_fn();
            let now = unix_now() as i64;
            let tx = self
                .conn
                .transaction()
                .map_err(|e| StreamError::Database(e.to_string()))?;
            {
                let mut stmt_sig = tx
                    .prepare_cached(
                        "INSERT INTO perturbation_signatures (id, signature_json)
                         VALUES (?1, ?2)
                         ON CONFLICT(id) DO UPDATE SET signature_json = excluded.signature_json",
                    )
                    .map_err(|e| StreamError::Database(e.to_string()))?;
                for sig in &seed_kg.signatures {
                    let sj = serde_json::to_string(sig)?;
                    stmt_sig
                        .execute(params![sig.id, sj])
                        .map_err(|e| StreamError::Database(e.to_string()))?;
                }
                let rev_hash = hex_sha256(
                    format!("sig-backfill:{now}:{}", seed_kg.signatures.len()).as_bytes(),
                );
                tx.execute(
                    "INSERT INTO kv_meta (key, value) VALUES ('kg_revision', ?1)
                     ON CONFLICT(key) DO UPDATE SET value = excluded.value",
                    params![rev_hash],
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;
            }
            tx.commit().map_err(|e| StreamError::Database(e.to_string()))?;
            return Ok(true);
        }

        Ok(false)
    }

    /// Atomically merges an entire `KnowledgeGraph` delta into SQLite inside a single WAL transaction,
    /// preserving stronger edge weights, existing topological diagrams, and updating the KG revision hash.
    pub fn upsert_knowledge_graph(&mut self, delta: &KnowledgeGraph) -> Result<(), StreamError> {
        let now = unix_now() as i64;
        let tx = self
            .conn
            .transaction()
            .map_err(|e| StreamError::Database(e.to_string()))?;

        {
            let mut stmt_node = tx
                .prepare_cached(
                    "INSERT INTO entities (id, display_name, kind_json, location_json, family_category, external_id, provenance_source)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
                     ON CONFLICT(id) DO UPDATE SET
                       display_name = CASE WHEN entities.display_name = entities.id THEN excluded.display_name ELSE entities.display_name END,
                       provenance_source = excluded.provenance_source",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for node in delta.nodes.values() {
                let kind_j = serde_json::to_string(&node.kind)?;
                let loc_j = serde_json::to_string(&node.location)?;
                stmt_node
                    .execute(params![
                        node.id,
                        node.display_name,
                        kind_j,
                        loc_j,
                        node.family_category,
                        node.external_id,
                        node.provenance_source
                    ])
                    .map_err(|e| StreamError::Database(e.to_string()))?;
            }

            let mut stmt_edge = tx
                .prepare_cached(
                    "INSERT INTO causal_edges (source, target, edge_class_json, sign_i8, weight, mechanism, source_db, pmids_json)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
                     ON CONFLICT(source, target, edge_class_json) DO UPDATE SET
                       weight = MAX(causal_edges.weight, excluded.weight),
                       sign_i8 = CASE WHEN causal_edges.sign_i8 = 0 THEN excluded.sign_i8 ELSE causal_edges.sign_i8 END,
                       source_db = excluded.source_db",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for edge in &delta.edges {
                let cls_j = serde_json::to_string(&edge.edge_class)?;
                let pmids_j = serde_json::to_string(&edge.pmids)?;
                stmt_edge
                    .execute(params![
                        edge.source,
                        edge.target,
                        cls_j,
                        i64::from(edge.sign.as_i8()),
                        edge.weight,
                        edge.mechanism,
                        edge.source_db,
                        pmids_j
                    ])
                    .map_err(|e| StreamError::Database(e.to_string()))?;
            }

            let mut stmt_pw = tx
                .prepare_cached(
                    "INSERT INTO canonical_pathways (
                        id, name, source_db, category, url, members_json, expected_signs_json,
                        topology_json, coords_json, phospho_json, cell_context_label, has_coords, updated_at_unix
                     ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)
                     ON CONFLICT(id) DO UPDATE SET
                        name = excluded.name,
                        members_json = excluded.members_json,
                        topology_json = CASE WHEN excluded.topology_json != '[]' THEN excluded.topology_json ELSE canonical_pathways.topology_json END,
                        coords_json = CASE WHEN excluded.has_coords = 1 THEN excluded.coords_json ELSE canonical_pathways.coords_json END,
                        has_coords = MAX(canonical_pathways.has_coords, excluded.has_coords),
                        updated_at_unix = excluded.updated_at_unix",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for pw in &delta.pathways {
                let members_j = serde_json::to_string(&pw.members)?;
                let signs_j = serde_json::to_string(&pw.member_expected_signs)?;
                let topo_j = serde_json::to_string(&pw.topology_edges)?;
                let coords_j = serde_json::to_string(&pw.node_coords)?;
                let phospho_j = serde_json::to_string(&pw.phospho_nodes)?;
                let has_coords = i64::from(!pw.node_coords.is_empty());
                stmt_pw
                    .execute(params![
                        pw.id,
                        pw.name,
                        pw.source_db,
                        pw.category,
                        pw.url,
                        members_j,
                        signs_j,
                        topo_j,
                        coords_j,
                        phospho_j,
                        pw.cell_context_label,
                        has_coords,
                        now
                    ])
                    .map_err(|e| StreamError::Database(e.to_string()))?;
            }

            let mut stmt_sig = tx
                .prepare_cached(
                    "INSERT INTO perturbation_signatures (id, signature_json)
                     VALUES (?1, ?2)
                     ON CONFLICT(id) DO UPDATE SET signature_json = excluded.signature_json",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for sig in &delta.signatures {
                let sj = serde_json::to_string(sig)?;
                stmt_sig
                    .execute(params![sig.id, sj])
                    .map_err(|e| StreamError::Database(e.to_string()))?;
            }

            let rev_seed = format!(
                "rev:{}:{}:{}:{}",
                now,
                delta.nodes.len(),
                delta.edges.len(),
                delta.pathways.len()
            );
            let rev_hash = hex_sha256(rev_seed.as_bytes());
            tx.execute(
                "INSERT INTO kv_meta (key, value) VALUES ('kg_revision', ?1)
                 ON CONFLICT(key) DO UPDATE SET value = excluded.value",
                params![rev_hash],
            )
            .map_err(|e| StreamError::Database(e.to_string()))?;
        }

        tx.commit().map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(())
    }

    /// Loads the complete `KnowledgeGraph` from SQLite and rebuilds adjacency indices.
    pub fn load_knowledge_graph(&self) -> Result<KnowledgeGraph, StreamError> {
        let mut kg = KnowledgeGraph::new();

        {
            let mut stmt = self
                .conn
                .prepare(
                    "SELECT id, display_name, kind_json, location_json, family_category, external_id, provenance_source FROM entities",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;
            let rows = stmt
                .query_map([], |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, String>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, String>(3)?,
                        r.get::<_, String>(4)?,
                        r.get::<_, Option<String>>(5)?,
                        r.get::<_, String>(6)?,
                    ))
                })
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for item in rows {
                let (id, display_name, kj, lj, fam, ext, prov) =
                    item.map_err(|e| StreamError::Database(e.to_string()))?;
                let kind: EntityKind = serde_json::from_str(&kj).unwrap_or(EntityKind::GeneOther);
                let location: SubcellularLocation =
                    serde_json::from_str(&lj).unwrap_or(SubcellularLocation::Cytoplasm);
                kg.nodes.insert(
                    id.clone(),
                    EntityNode {
                        id,
                        display_name,
                        kind,
                        location,
                        family_category: fam,
                        external_id: ext,
                        provenance_source: prov,
                    },
                );
            }
        }

        {
            let mut stmt = self
                .conn
                .prepare(
                    "SELECT source, target, edge_class_json, sign_i8, weight, mechanism, source_db, pmids_json FROM causal_edges",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;
            let rows = stmt
                .query_map([], |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, String>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, i64>(3)?,
                        r.get::<_, f64>(4)?,
                        r.get::<_, String>(5)?,
                        r.get::<_, String>(6)?,
                        r.get::<_, String>(7)?,
                    ))
                })
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for item in rows {
                let (source, target, cj, si, weight, mechanism, source_db, pj) =
                    item.map_err(|e| StreamError::Database(e.to_string()))?;
                let edge_class: EdgeClass =
                    serde_json::from_str(&cj).unwrap_or(EdgeClass::Activation);
                let pmids: Vec<String> = serde_json::from_str(&pj).unwrap_or_default();
                kg.edges.push(CausalEdge {
                    source,
                    target,
                    edge_class,
                    sign: EdgeSign::from_i8(si.clamp(-1, 1) as i8),
                    weight,
                    mechanism,
                    source_db,
                    pmids,
                });
            }
        }

        {
            let mut stmt = self
                .conn
                .prepare(
                    "SELECT id, name, source_db, category, url, members_json, expected_signs_json,
                            topology_json, coords_json, phospho_json, cell_context_label
                     FROM canonical_pathways",
                )
                .map_err(|e| StreamError::Database(e.to_string()))?;
            let rows = stmt
                .query_map([], |r| {
                    Ok((
                        r.get::<_, String>(0)?,
                        r.get::<_, String>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, String>(3)?,
                        r.get::<_, String>(4)?,
                        r.get::<_, String>(5)?,
                        r.get::<_, String>(6)?,
                        r.get::<_, String>(7)?,
                        r.get::<_, String>(8)?,
                        r.get::<_, String>(9)?,
                        r.get::<_, String>(10)?,
                    ))
                })
                .map_err(|e| StreamError::Database(e.to_string()))?;

            for item in rows {
                let (id, name, source_db, category, url, mj, sj, tj, cj, phj, cell_lbl) =
                    item.map_err(|e| StreamError::Database(e.to_string()))?;
                kg.pathways.push(PathwayDef {
                    id,
                    name,
                    source_db,
                    category,
                    url,
                    members: serde_json::from_str(&mj).unwrap_or_default(),
                    member_expected_signs: serde_json::from_str(&sj).unwrap_or_default(),
                    topology_edges: serde_json::from_str(&tj).unwrap_or_default(),
                    node_coords: serde_json::from_str(&cj).unwrap_or_default(),
                    phospho_nodes: serde_json::from_str(&phj).unwrap_or_default(),
                    cell_context_label: cell_lbl,
                });
            }
        }

        {
            let mut stmt = self
                .conn
                .prepare("SELECT signature_json FROM perturbation_signatures")
                .map_err(|e| StreamError::Database(e.to_string()))?;
            let rows = stmt
                .query_map([], |r| r.get::<_, String>(0))
                .map_err(|e| StreamError::Database(e.to_string()))?;
            for item in rows {
                let sj = item.map_err(|e| StreamError::Database(e.to_string()))?;
                if let Ok(sig) = serde_json::from_str::<PerturbationSignature>(&sj) {
                    kg.signatures.push(sig);
                }
            }
        }

        kg.rebuild_indices();
        Ok(kg)
    }

    /// Retrieves stream HTTP ETag / Last-Modified / SHA-256 cache metadata.
    pub fn get_stream_meta(&self, stream_id: &str) -> Result<Option<StreamCacheMeta>, StreamError> {
        self.conn
            .query_row(
                "SELECT stream_id, etag, last_modified, content_sha256, records_count, synced_at_unix, ttl_seconds
                 FROM stream_sync_meta WHERE stream_id = ?1",
                params![stream_id],
                |r| {
                    Ok(StreamCacheMeta {
                        stream_id: r.get(0)?,
                        etag: r.get(1)?,
                        last_modified: r.get(2)?,
                        content_sha256: r.get(3)?,
                        records_count: r.get::<_, i64>(4)? as usize,
                        synced_at_unix: r.get::<_, i64>(5)? as u64,
                        ttl_seconds: r.get::<_, i64>(6)? as u64,
                    })
                },
            )
            .optional()
            .map_err(|e| StreamError::Database(e.to_string()))
    }

    /// Updates sync metadata after a successful download or 304 Not Modified revalidation.
    pub fn put_stream_meta(&self, meta: &StreamCacheMeta) -> Result<(), StreamError> {
        self.conn
            .execute(
                "INSERT INTO stream_sync_meta (stream_id, etag, last_modified, content_sha256, records_count, synced_at_unix, ttl_seconds)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
                 ON CONFLICT(stream_id) DO UPDATE SET
                   etag = COALESCE(excluded.etag, stream_sync_meta.etag),
                   last_modified = COALESCE(excluded.last_modified, stream_sync_meta.last_modified),
                   content_sha256 = excluded.content_sha256,
                   records_count = excluded.records_count,
                   synced_at_unix = excluded.synced_at_unix,
                   ttl_seconds = excluded.ttl_seconds",
                params![
                    meta.stream_id,
                    meta.etag,
                    meta.last_modified,
                    meta.content_sha256,
                    meta.records_count as i64,
                    meta.synced_at_unix as i64,
                    meta.ttl_seconds as i64
                ],
            )
            .map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(())
    }

    pub fn kg_revision(&self) -> String {
        self.conn
            .query_row(
                "SELECT value FROM kv_meta WHERE key = 'kg_revision'",
                [],
                |r| r.get::<_, String>(0),
            )
            .unwrap_or_else(|_| "seed-v0".to_string())
    }

    /// Computes a deterministic SHA-256 key over the dataset measurements, KG revision, and statistical cutoffs.
    pub fn compute_analysis_cache_key(
        &self,
        dataset: &ExpressionDataset,
        cfg: &CoreAnalysisConfig,
    ) -> String {
        let mut h = Sha256::new();
        h.update(self.kg_revision().as_bytes());
        h.update(b"|ca-core-v2|");
        h.update(dataset.name.as_bytes());
        h.update(format!("|{:.4}|{:.6}|{}|{:.4}|", cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, cfg.cna_max_depth, cfg.cna_weight_cutoff).as_bytes());
        for r in &dataset.rows {
            h.update(r.symbol.as_bytes());
            h.update(r.log2_fc.to_le_bytes());
            h.update(r.fdr.to_le_bytes());
        }
        let d = h.finalize();
        let mut s = String::with_capacity(64);
        for b in d {
            use std::fmt::Write;
            let _ = write!(&mut s, "{b:02x}");
        }
        s
    }

    /// Looks up a previously computed `CoreAnalysisReport` in SQLite (`< 1 ms` cache hit).
    pub fn get_cached_analysis(
        &self,
        dataset: &ExpressionDataset,
        cfg: &CoreAnalysisConfig,
    ) -> Result<Option<CoreAnalysisReport>, StreamError> {
        let key = self.compute_analysis_cache_key(dataset, cfg);
        let blob_opt: Option<Vec<u8>> = self
            .conn
            .query_row(
                "SELECT report_gzip FROM analysis_result_cache WHERE cache_key = ?1",
                params![key],
                |r| r.get(0),
            )
            .optional()
            .map_err(|e| StreamError::Database(e.to_string()))?;

        if let Some(gz_bytes) = blob_opt {
            let json_bytes = gzip_decompress(&gz_bytes)?;
            let rep: CoreAnalysisReport = serde_json::from_slice(&json_bytes)?;
            Ok(Some(rep))
        } else {
            Ok(None)
        }
    }

    /// Stores a `CoreAnalysisReport` compressed in SQLite `analysis_result_cache`.
    pub fn put_cached_analysis(
        &self,
        dataset: &ExpressionDataset,
        cfg: &CoreAnalysisConfig,
        report: &CoreAnalysisReport,
    ) -> Result<(), StreamError> {
        let key = self.compute_analysis_cache_key(dataset, cfg);
        let json_bytes = serde_json::to_vec(report)?;
        let gz_bytes = gzip_compress(&json_bytes)?;
        self.conn
            .execute(
                "INSERT INTO analysis_result_cache (cache_key, dataset_name, computed_at_unix, report_gzip)
                 VALUES (?1, ?2, ?3, ?4)
                 ON CONFLICT(cache_key) DO UPDATE SET
                   computed_at_unix = excluded.computed_at_unix,
                   report_gzip = excluded.report_gzip",
                params![key, dataset.name, unix_now() as i64, gz_bytes],
            )
            .map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(())
    }

    /// Clears cached Core Analysis reports while retaining all downloaded biological knowledge streams.
    pub fn clear_analysis_cache(&self) -> Result<usize, StreamError> {
        let removed = self
            .conn
            .execute("DELETE FROM analysis_result_cache", [])
            .map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(removed)
    }

    /// Runs `PRAGMA wal_checkpoint(TRUNCATE)` and `VACUUM` to compact the SQLite file on disk.
    pub fn vacuum_and_optimize(&self) -> Result<(), StreamError> {
        self.conn
            .execute_batch("PRAGMA wal_checkpoint(TRUNCATE); PRAGMA optimize; VACUUM;")
            .map_err(|e| StreamError::Database(e.to_string()))?;
        Ok(())
    }

    /// Returns live table counts, WAL status, and disk footprint for GUI display.
    pub fn stats(&self) -> Result<SqliteCacheStats, StreamError> {
        let count_of = |sql: &str| -> Result<usize, StreamError> {
            self.conn
                .query_row(sql, [], |r| r.get::<_, i64>(0))
                .map(|v| v.max(0) as usize)
                .map_err(|e| StreamError::Database(e.to_string()))
        };

        let journal: String = self
            .conn
            .query_row("PRAGMA journal_mode", [], |r| r.get(0))
            .unwrap_or_default();

        let file_size_bytes = std::fs::metadata(&self.path).map(|m| m.len()).unwrap_or(0);

        Ok(SqliteCacheStats {
            db_path: self.path.display().to_string(),
            file_size_bytes,
            wal_mode: journal.eq_ignore_ascii_case("wal"),
            cached_entities: count_of("SELECT COUNT(*) FROM entities")?,
            cached_edges: count_of("SELECT COUNT(*) FROM causal_edges")?,
            cached_pathways: count_of("SELECT COUNT(*) FROM canonical_pathways")?,
            cached_diagrams_with_coords: count_of(
                "SELECT COUNT(*) FROM canonical_pathways WHERE has_coords = 1",
            )?,
            cached_signatures: count_of("SELECT COUNT(*) FROM perturbation_signatures")?,
            synced_streams: count_of("SELECT COUNT(*) FROM stream_sync_meta")?,
            cached_analysis_reports: count_of("SELECT COUNT(*) FROM analysis_result_cache")?,
            kg_revision_hash: self.kg_revision().chars().take(12).collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ca_core::{builtin_cohorts, run_core_analysis, seed_kb::build_seed_knowledge_graph};

    #[test]
    fn sqlite_roundtrip_and_analysis_cache_hit() {
        let mut store = SqliteKnowledgeStore::open_in_memory().unwrap();
        let seeded = store.ensure_seeded(build_seed_knowledge_graph).unwrap();
        assert!(seeded);

        let kg = store.load_knowledge_graph().unwrap();
        assert!(kg.pathways.len() >= 6);
        assert!(kg.edges.len() >= 100);

        let ds = &builtin_cohorts()[0];
        let cfg = CoreAnalysisConfig {
            cna_permutations: 150,
            ..CoreAnalysisConfig::default()
        };
        assert!(store.get_cached_analysis(ds, &cfg).unwrap().is_none());

        let rep = run_core_analysis(ds, &kg, &cfg);
        store.put_cached_analysis(ds, &cfg, &rep).unwrap();

        let cached = store
            .get_cached_analysis(ds, &cfg)
            .unwrap()
            .expect("Expected SQLite cache hit");
        assert_eq!(cached.canonical_pathways.len(), rep.canonical_pathways.len());
    }

    #[test]
    fn sqlite_backfill_restores_expanded_seed_signatures() {
        let mut store = SqliteKnowledgeStore::open_in_memory().unwrap();
        assert!(store.ensure_seeded(build_seed_knowledge_graph).unwrap());
        let rev_before = store.kg_revision();

        store
            .conn
            .execute(
                "DELETE FROM perturbation_signatures WHERE id NOT IN (
                    SELECT id FROM perturbation_signatures ORDER BY id LIMIT 7
                 )",
                [],
            )
            .unwrap();
        assert_eq!(store.stats().unwrap().cached_signatures, 7);

        let backfilled = store.ensure_seeded(build_seed_knowledge_graph).unwrap();
        assert!(backfilled, "Expected ensure_seeded to backfill signatures");
        let kg = store.load_knowledge_graph().unwrap();
        assert!(
            kg.signatures.len() >= 25,
            "Expected >= 25 signatures after backfill, got {}",
            kg.signatures.len()
        );
        assert_ne!(store.kg_revision(), rev_before);
    }
}
