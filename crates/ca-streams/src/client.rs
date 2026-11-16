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

use crate::batch_downloader::{sync_reactome_gmt_bundle, sync_wikipathways_gpml_bundle};
use crate::parsers::parse_stream_payload;
use crate::sqlite_store::{unix_now, SqliteKnowledgeStore, StreamCacheMeta};
use crate::{StreamDescriptor, StreamError, StreamFormat};
use ca_core::KnowledgeGraph;
use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::Path;
use std::time::Duration;

const DEFAULT_STREAM_TTL_SECS: u64 = 7 * 24 * 3600;

fn compute_sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let mut out = String::with_capacity(64);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(&mut out, "{b:02x}");
    }
    out
}

fn decode_maybe_gzip(raw_bytes: &[u8]) -> Result<String, StreamError> {
    if raw_bytes.starts_with(&[0x1f, 0x8b]) {
        let mut gz = GzDecoder::new(raw_bytes);
        let mut s = String::new();
        gz.read_to_string(&mut s)?;
        Ok(s)
    } else {
        String::from_utf8(raw_bytes.to_vec())
            .map_err(|e| StreamError::Parse(format!("Invalid UTF-8 in stream payload: {e}")))
    }
}

/// Synchronizes a single stream into `SqliteKnowledgeStore` using TTL, HTTP `ETag` / `Last-Modified`,
/// and SHA-256 content deduplication. Returns `(records_count, was_cache_hit)`.
pub fn fetch_remote_stream_conditional(
    store: &mut SqliteKnowledgeStore,
    desc: &StreamDescriptor,
    force_refresh: bool,
) -> Result<(usize, bool), StreamError> {
    if desc.id == "wikipathways_gmt" {
        return sync_wikipathways_gpml_bundle(store, force_refresh);
    }
    if desc.id == "reactome_content" {
        return sync_reactome_gmt_bundle(store, force_refresh);
    }

    let now = unix_now();
    let prev_meta = store.get_stream_meta(&desc.id)?;

    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.is_fresh(now) {
                return Ok((m.records_count, true));
            }
        }
    }

    let agent = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(12))
        .timeout_read(Duration::from_secs(60))
        .user_agent("causal-atlas/0.2.0 (Open-Source Clean-Room Pathway Workbench)")
        .build();

    let (raw_bytes, etag, last_mod) = if desc.format == StreamFormat::OpenTargetsGraphQlJson {
        let query = serde_json::json!({
            "query": "query DiseaseTargets($efoId: String!) { disease(efoId: $efoId) { id name associatedTargets(page: { index: 0, size: 50 }) { rows { score target { approvedSymbol } } } } }",
            "variables": { "efoId": "EFO_0000305" }
        });
        let resp = agent
            .post(&desc.endpoint_url)
            .set("Content-Type", "application/json")
            .send_json(query)
            .map_err(|e| StreamError::Network(e.to_string()))?;
        let mut buf = Vec::new();
        resp.into_reader().read_to_end(&mut buf)?;
        (buf, None, None)
    } else {
        let mut req = agent.get(&desc.endpoint_url);
        if !force_refresh {
            if let Some(ref m) = prev_meta {
                if let Some(ref et) = m.etag {
                    req = req.set("If-None-Match", et);
                }
                if let Some(ref lm) = m.last_modified {
                    req = req.set("If-Modified-Since", lm);
                }
            }
        }
        let resp = req.call().map_err(|e| StreamError::Network(e.to_string()))?;
        if resp.status() == 304 {
            if let Some(mut m) = prev_meta {
                m.synced_at_unix = now;
                store.put_stream_meta(&m)?;
                return Ok((m.records_count, true));
            }
        }
        let etag = resp.header("ETag").map(String::from);
        let lm = resp.header("Last-Modified").map(String::from);
        let mut buf = Vec::new();
        resp.into_reader().read_to_end(&mut buf)?;
        (buf, etag, lm)
    };

    let sha256 = compute_sha256_hex(&raw_bytes);
    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.content_sha256 == sha256 {
                let mut updated = m.clone();
                updated.synced_at_unix = now;
                store.put_stream_meta(&updated)?;
                return Ok((updated.records_count, true));
            }
        }
    }

    let text = decode_maybe_gzip(&raw_bytes)?;
    let (kg, count) = parse_stream_payload(desc.format, &desc.name, &text)?;
    store.upsert_knowledge_graph(&kg)?;
    store.put_stream_meta(&StreamCacheMeta {
        stream_id: desc.id.clone(),
        etag,
        last_modified: last_mod,
        content_sha256: sha256,
        records_count: count,
        synced_at_unix: now,
        ttl_seconds: DEFAULT_STREAM_TTL_SECS,
    })?;

    Ok((count, false))
}

/// Stateless fetch helper returning `(KnowledgeGraph, records_count, sha256)`.
pub fn fetch_remote_stream(
    desc: &StreamDescriptor,
) -> Result<(KnowledgeGraph, usize, String), StreamError> {
    let mut mem = SqliteKnowledgeStore::open_in_memory()?;
    let (cnt, _) = fetch_remote_stream_conditional(&mut mem, desc, true)?;
    let sha = mem
        .get_stream_meta(&desc.id)?
        .map(|m| m.content_sha256)
        .unwrap_or_default();
    let kg = mem.load_knowledge_graph()?;
    Ok((kg, cnt, sha))
}

/// Imports a local air-gapped open-source stream file (`.gmt`, `.gpml`, `.zip`, `.tsv`, `.gaf.gz`, `.json`).
pub fn import_local_stream_file(
    path: &Path,
    explicit_format: Option<StreamFormat>,
) -> Result<(KnowledgeGraph, usize, String), StreamError> {
    let raw_bytes = std::fs::read(path)?;
    let sha256 = compute_sha256_hex(&raw_bytes);
    let file_name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("local_stream")
        .to_lowercase();

    // Support direct import of downloaded WikiPathways GPML `.zip` or Reactome `.gmt.zip` archives!
    if file_name.ends_with(".zip") {
        use crate::parsers::{parse_gmt_pathways, parse_gpml_pathway};
        use std::io::Cursor;
        let mut archive = zip::ZipArchive::new(Cursor::new(raw_bytes))
            .map_err(|e| StreamError::Parse(format!("Invalid zip archive: {e}")))?;
        let mut combined = KnowledgeGraph::new();
        let mut total_items = 0usize;
        for i in 0..archive.len() {
            let mut entry = archive
                .by_index(i)
                .map_err(|e| StreamError::Parse(e.to_string()))?;
            let ename = entry.name().to_lowercase();
            let mut text = String::new();
            if entry.read_to_string(&mut text).is_ok() {
                if ename.ends_with(".gpml") {
                    if let Ok((sub, _)) = parse_gpml_pathway("WikiPathways", &text) {
                        combined.merge(sub);
                        total_items += 1;
                    }
                } else if ename.ends_with(".gmt") {
                    if let Ok((sub, c)) = parse_gmt_pathways("Imported Zip GMT", &text) {
                        combined.merge(sub);
                        total_items += c;
                    }
                }
            }
        }
        combined.rebuild_indices();
        return Ok((combined, total_items, sha256));
    }

    let inferred = explicit_format.unwrap_or_else(|| {
        if file_name.contains(".gmt") {
            StreamFormat::GmtPathways
        } else if file_name.contains(".gpml") || file_name.ends_with(".xml") {
            StreamFormat::GpmlXml
        } else if file_name.contains("signor") {
            StreamFormat::SignorTsv
        } else if file_name.contains(".gaf") {
            StreamFormat::GoGafTsv
        } else if file_name.contains("string") {
            StreamFormat::StringTsv
        } else if file_name.ends_with(".json") || file_name.ends_with(".json.gz") {
            StreamFormat::LincsSignaturesJson
        } else {
            StreamFormat::OmniPathTsv
        }
    });

    let text = decode_maybe_gzip(&raw_bytes)?;
    let (kg, count) = parse_stream_payload(inferred, &file_name, &text)?;
    Ok((kg, count, sha256))
}
