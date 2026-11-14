#![forbid(unsafe_code)]

use crate::client::fetch_remote_stream_conditional;
use crate::parsers::{parse_gmt_pathways, parse_gpml_pathway};
use crate::sqlite_store::{unix_now, SqliteKnowledgeStore, StreamCacheMeta};
use crate::{default_stream_registry, StreamError};
use ca_core::{EdgeSign, KnowledgeGraph, PathwayDef, PathwayReactionEdge};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::io::{Cursor, Read};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;
use zip::ZipArchive;

/// Progress event emitted during batch downloading and incremental stream updates.
#[derive(Debug, Clone)]
pub struct BatchProgressEvent {
    pub step_index: usize,
    pub total_steps: usize,
    pub source_name: String,
    pub message: String,
    pub cache_hit: bool,
    pub records_processed: usize,
}

#[derive(Debug, Clone, Default)]
pub struct BatchSyncSummary {
    pub streams_checked: usize,
    pub cache_hits: usize,
    pub streams_downloaded: usize,
    pub pathways_ingested: usize,
    pub causal_edges_ingested: usize,
    pub errors: Vec<String>,
}

const WIKIPATHWAYS_GPML_INDEX_URL: &str = "https://data.wikipathways.org/current/gpml/";
const REACTOME_GMT_ZIP_URL: &str = "https://reactome.org/download/current/ReactomePathways.gmt.zip";
const DEFAULT_TTL_SECS: u64 = 7 * 24 * 3600;

fn sha256_hex(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    let d = h.finalize();
    let mut out = String::with_capacity(64);
    for b in d {
        use std::fmt::Write;
        let _ = write!(&mut out, "{b:02x}");
    }
    out
}

fn build_http_agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(15))
        .timeout_read(Duration::from_secs(90))
        .user_agent("causal-atlas/0.2.0 (Open-Source Rust Pathway Workbench)")
        .build()
}

/// Discovers the latest `wikipathways-YYYYMMDD-gpml-Homo_sapiens.zip` archive from
/// `https://data.wikipathways.org/current/gpml/`, downloads it with ETag/SHA-256 caching,
/// extracts all 866+ human `.gpml` XML pathway files in memory, and upserts them into SQLite.
pub fn sync_wikipathways_gpml_bundle(
    store: &mut SqliteKnowledgeStore,
    force_refresh: bool,
) -> Result<(usize, bool), StreamError> {
    let stream_id = "wikipathways_gpml_bundle";
    let now = unix_now();
    let prev_meta = store.get_stream_meta(stream_id)?;

    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.is_fresh(now) {
                return Ok((m.records_count, true));
            }
        }
    }

    let agent = build_http_agent();
    let index_html = agent
        .get(WIKIPATHWAYS_GPML_INDEX_URL)
        .call()
        .map_err(|e| StreamError::Network(e.to_string()))?
        .into_string()
        .map_err(StreamError::Io)?;

    // Locate `wikipathways-YYYYMMDD-gpml-Homo_sapiens.zip` href in directory listing
    let zip_filename = extract_homo_sapiens_gpml_zip_name(&index_html).ok_or_else(|| {
        StreamError::Parse(
            "Could not find Homo_sapiens GPML zip link at data.wikipathways.org/current/gpml/"
                .to_string(),
        )
    })?;

    let zip_url = if zip_filename.starts_with("http") {
        zip_filename
    } else {
        format!("{WIKIPATHWAYS_GPML_INDEX_URL}{zip_filename}")
    };

    let mut req = agent.get(&zip_url);
    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if let Some(ref etag) = m.etag {
                req = req.set("If-None-Match", etag);
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
    let last_mod = resp.header("Last-Modified").map(String::from);
    let mut zip_bytes = Vec::new();
    resp.into_reader().read_to_end(&mut zip_bytes)?;
    let sha = sha256_hex(&zip_bytes);

    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.content_sha256 == sha {
                let mut updated = m.clone();
                updated.synced_at_unix = now;
                store.put_stream_meta(&updated)?;
                return Ok((updated.records_count, true));
            }
        }
    }

    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))
        .map_err(|e| StreamError::Parse(format!("Invalid GPML zip archive: {e}")))?;

    let mut combined = KnowledgeGraph::new();
    let mut pw_count = 0usize;

    for i in 0..archive.len() {
        let mut entry = archive
            .by_index(i)
            .map_err(|e| StreamError::Parse(e.to_string()))?;
        let entry_name = entry.name().to_string();
        if !entry_name.ends_with(".gpml") {
            continue;
        }
        let mut xml = String::new();
        if entry.read_to_string(&mut xml).is_ok() {
            if let Ok((mut sub_kg, _)) = parse_gpml_pathway("WikiPathways", &xml) {
                // Use exact WP identifier from filename (e.g. `..._WP445_...gpml`) when present
                if let Some(wp_id) = extract_wp_id(&entry_name) {
                    if let Some(pw) = sub_kg.pathways.first_mut() {
                        pw.url = format!("https://www.wikipathways.org/pathways/{wp_id}.html");
                        pw.id = wp_id;
                    }
                }
                for (id, n) in sub_kg.nodes {
                    combined.nodes.entry(id).or_insert(n);
                }
                combined.pathways.extend(sub_kg.pathways);
                pw_count += 1;
            }
        }
    }

    combined.rebuild_indices();
    store.upsert_knowledge_graph(&combined)?;
    store.put_stream_meta(&StreamCacheMeta {
        stream_id: stream_id.to_string(),
        etag,
        last_modified: last_mod,
        content_sha256: sha,
        records_count: pw_count,
        synced_at_unix: now,
        ttl_seconds: DEFAULT_TTL_SECS,
    })?;

    Ok((pw_count, false))
}

/// Downloads `ReactomePathways.gmt.zip` (all 2,647+ curated human pathways) with ETag/SHA-256 caching,
/// extracts `ReactomePathways.gmt` in memory, and batch-upserts into SQLite.
pub fn sync_reactome_gmt_bundle(
    store: &mut SqliteKnowledgeStore,
    force_refresh: bool,
) -> Result<(usize, bool), StreamError> {
    let stream_id = "reactome_gmt_bundle";
    let now = unix_now();
    let prev_meta = store.get_stream_meta(stream_id)?;

    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.is_fresh(now) {
                return Ok((m.records_count, true));
            }
        }
    }

    let agent = build_http_agent();
    let mut req = agent.get(REACTOME_GMT_ZIP_URL);
    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if let Some(ref etag) = m.etag {
                req = req.set("If-None-Match", etag);
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
    let last_mod = resp.header("Last-Modified").map(String::from);
    let mut zip_bytes = Vec::new();
    resp.into_reader().read_to_end(&mut zip_bytes)?;
    let sha = sha256_hex(&zip_bytes);

    if !force_refresh {
        if let Some(ref m) = prev_meta {
            if m.content_sha256 == sha {
                let mut updated = m.clone();
                updated.synced_at_unix = now;
                store.put_stream_meta(&updated)?;
                return Ok((updated.records_count, true));
            }
        }
    }

    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))
        .map_err(|e| StreamError::Parse(format!("Invalid Reactome GMT zip: {e}")))?;

    let mut gmt_text = String::new();
    for i in 0..archive.len() {
        let mut entry = archive
            .by_index(i)
            .map_err(|e| StreamError::Parse(e.to_string()))?;
        if entry.name().ends_with(".gmt") {
            entry.read_to_string(&mut gmt_text)?;
            break;
        }
    }

    let (kg, count) = parse_gmt_pathways("Reactome", &gmt_text)?;
    store.upsert_knowledge_graph(&kg)?;
    store.put_stream_meta(&StreamCacheMeta {
        stream_id: stream_id.to_string(),
        etag,
        last_modified: last_mod,
        content_sha256: sha,
        records_count: count,
        synced_at_unix: now,
        ttl_seconds: DEFAULT_TTL_SECS,
    })?;

    Ok((count, false))
}

/// Fetches and caches a Reactome `DiagramJS` layout JSON (`/ContentService/exporter/diagram/{stId}.json`
/// or `/download/current/diagram/{stId}.json`) on demand for any `R-HSA-*` pathway and persists its
/// normalized 2D coordinates and reaction edges into SQLite.
pub fn fetch_and_cache_reactome_diagram(
    store: &mut SqliteKnowledgeStore,
    pathway: &PathwayDef,
) -> Result<PathwayDef, StreamError> {
    if !pathway.id.starts_with("R-HSA-") {
        return Ok(pathway.clone());
    }
    let url = format!(
        "https://reactome.org/download/current/diagram/{}.json",
        pathway.id
    );
    let agent = build_http_agent();
    let resp = agent
        .get(&url)
        .call()
        .map_err(|e| StreamError::Network(e.to_string()))?;
    let val: serde_json::Value = resp
        .into_json()
        .map_err(|e| StreamError::Parse(e.to_string()))?;

    let mut id_to_sym: BTreeMap<u64, String> = BTreeMap::new();
    let mut raw_xy: Vec<(String, f32, f32)> = Vec::new();

    if let Some(nodes) = val.get("nodes").and_then(|v| v.as_array()) {
        for n in nodes {
            let Some(nid) = n.get("id").and_then(|v| v.as_u64()) else {
                continue;
            };
            let Some(disp) = n.get("displayName").and_then(|v| v.as_str()) else {
                continue;
            };
            let sym = disp
                .split([' ', '[', '('])
                .next()
                .unwrap_or(disp)
                .trim()
                .to_uppercase();
            if sym.is_empty() {
                continue;
            }
            id_to_sym.insert(nid, sym.clone());
            if let (Some(x), Some(y)) = (
                n.get("x").and_then(|v| v.as_f64()),
                n.get("y").and_then(|v| v.as_f64()),
            ) {
                raw_xy.push((sym, x as f32, y as f32));
            }
        }
    }

    let mut updated = pathway.clone();
    if raw_xy.len() >= 2 {
        let min_x = raw_xy.iter().map(|p| p.1).fold(f32::INFINITY, f32::min);
        let max_x = raw_xy.iter().map(|p| p.1).fold(f32::NEG_INFINITY, f32::max);
        let min_y = raw_xy.iter().map(|p| p.2).fold(f32::INFINITY, f32::min);
        let max_y = raw_xy.iter().map(|p| p.2).fold(f32::NEG_INFINITY, f32::max);
        let sx = (max_x - min_x).max(1.0);
        let sy = (max_y - min_y).max(1.0);
        for (sym, x, y) in raw_xy {
            updated
                .node_coords
                .insert(sym, (0.08 + ((x - min_x) / sx) * 0.84, 0.08 + ((y - min_y) / sy) * 0.84));
        }
    }

    if let Some(edges) = val.get("edges").and_then(|v| v.as_array()) {
        for e in edges {
            let inputs = e.get("inputs").and_then(|v| v.as_array());
            let outputs = e.get("outputs").and_then(|v| v.as_array());
            let inhibitors = e.get("inhibitors").and_then(|v| v.as_array());

            if let (Some(ins), Some(outs)) = (inputs, outputs) {
                for inp in ins {
                    let Some(in_id) = inp.get("id").and_then(|v| v.as_u64()) else {
                        continue;
                    };
                    let Some(src_sym) = id_to_sym.get(&in_id) else {
                        continue;
                    };
                    for outp in outs {
                        let Some(out_id) = outp.get("id").and_then(|v| v.as_u64()) else {
                            continue;
                        };
                        let Some(dst_sym) = id_to_sym.get(&out_id) else {
                            continue;
                        };
                        if src_sym != dst_sym {
                            updated.topology_edges.push(PathwayReactionEdge {
                                from: src_sym.clone(),
                                to: dst_sym.clone(),
                                sign: EdgeSign::Activating,
                                interaction_label: "reaction".to_string(),
                            });
                        }
                    }
                }
            }
            if let (Some(inhs), Some(outs)) = (inhibitors, outputs) {
                for inh in inhs {
                    let Some(inh_id) = inh.get("id").and_then(|v| v.as_u64()) else {
                        continue;
                    };
                    let Some(src_sym) = id_to_sym.get(&inh_id) else {
                        continue;
                    };
                    for outp in outs {
                        let Some(out_id) = outp.get("id").and_then(|v| v.as_u64()) else {
                            continue;
                        };
                        let Some(dst_sym) = id_to_sym.get(&out_id) else {
                            continue;
                        };
                        if src_sym != dst_sym {
                            updated.topology_edges.push(PathwayReactionEdge {
                                from: src_sym.clone(),
                                to: dst_sym.clone(),
                                sign: EdgeSign::Inhibiting,
                                interaction_label: "inhibition".to_string(),
                            });
                        }
                    }
                }
            }
        }
    }

    let mut delta = KnowledgeGraph::new();
    delta.pathways.push(updated.clone());
    store.upsert_knowledge_graph(&delta)?;
    Ok(updated)
}

/// Orchestrates a full batch download & incremental update across WikiPathways GPML, Reactome GMT,
/// and all enabled causal regulatory streams (`OmniPath`, `CollecTRI`, `SIGNOR 3.0`, `STRING`, `Open Targets`),
/// emitting real-time `BatchProgressEvent` updates for the GUI or CLI progress bar.
pub fn run_batch_download_and_update<F>(
    store: &mut SqliteKnowledgeStore,
    include_pathway_bundles: bool,
    include_regulatory_streams: bool,
    force_refresh: bool,
    cancel_flag: Option<&AtomicBool>,
    mut on_progress: F,
) -> BatchSyncSummary
where
    F: FnMut(BatchProgressEvent),
{
    let mut summary = BatchSyncSummary::default();
    let reg_streams: Vec<_> = if include_regulatory_streams {
        default_stream_registry()
            .into_iter()
            .filter(|d| d.enabled_by_default && d.id != "wikipathways_gmt")
            .collect()
    } else {
        Vec::new()
    };

    let bundle_steps = if include_pathway_bundles { 2 } else { 0 };
    let total_steps = bundle_steps + reg_streams.len();
    let mut step = 0usize;

    if include_pathway_bundles {
        if cancel_flag.is_some_and(|c| c.load(Ordering::Relaxed)) {
            return summary;
        }
        step += 1;
        summary.streams_checked += 1;
        match sync_wikipathways_gpml_bundle(store, force_refresh) {
            Ok((cnt, hit)) => {
                if hit {
                    summary.cache_hits += 1;
                } else {
                    summary.streams_downloaded += 1;
                    summary.pathways_ingested += cnt;
                }
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: "WikiPathways Homo sapiens GPML Zip".to_string(),
                    message: if hit {
                        format!("SQLite cache hit ({cnt} GPML diagrams up to date)")
                    } else {
                        format!("Downloaded & cached {cnt} GPML coordinate diagrams")
                    },
                    cache_hit: hit,
                    records_processed: cnt,
                });
            }
            Err(e) => {
                let msg = format!("WikiPathways GPML: {e}");
                summary.errors.push(msg.clone());
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: "WikiPathways Homo sapiens GPML Zip".to_string(),
                    message: msg,
                    cache_hit: false,
                    records_processed: 0,
                });
            }
        }

        if cancel_flag.is_some_and(|c| c.load(Ordering::Relaxed)) {
            return summary;
        }
        step += 1;
        summary.streams_checked += 1;
        match sync_reactome_gmt_bundle(store, force_refresh) {
            Ok((cnt, hit)) => {
                if hit {
                    summary.cache_hits += 1;
                } else {
                    summary.streams_downloaded += 1;
                    summary.pathways_ingested += cnt;
                }
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: "Reactome Human Pathways GMT Zip".to_string(),
                    message: if hit {
                        format!("SQLite cache hit ({cnt} Reactome pathways up to date)")
                    } else {
                        format!("Downloaded & cached {cnt} Reactome human pathways")
                    },
                    cache_hit: hit,
                    records_processed: cnt,
                });
            }
            Err(e) => {
                let msg = format!("Reactome GMT: {e}");
                summary.errors.push(msg.clone());
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: "Reactome Human Pathways GMT Zip".to_string(),
                    message: msg,
                    cache_hit: false,
                    records_processed: 0,
                });
            }
        }
    }

    for desc in reg_streams {
        if cancel_flag.is_some_and(|c| c.load(Ordering::Relaxed)) {
            break;
        }
        step += 1;
        summary.streams_checked += 1;
        match fetch_remote_stream_conditional(store, &desc, force_refresh) {
            Ok((cnt, hit)) => {
                if hit {
                    summary.cache_hits += 1;
                } else {
                    summary.streams_downloaded += 1;
                    summary.causal_edges_ingested += cnt;
                }
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: desc.name.clone(),
                    message: if hit {
                        format!("SQLite cache hit ({cnt} records valid)")
                    } else {
                        format!("Synced & cached {cnt} records into SQLite")
                    },
                    cache_hit: hit,
                    records_processed: cnt,
                });
            }
            Err(e) => {
                let msg = format!("{}: {e}", desc.name);
                summary.errors.push(msg.clone());
                on_progress(BatchProgressEvent {
                    step_index: step,
                    total_steps,
                    source_name: desc.name,
                    message: msg,
                    cache_hit: false,
                    records_processed: 0,
                });
            }
        }
    }

    summary
}

pub fn extract_homo_sapiens_gpml_zip_name(html: &str) -> Option<String> {
    // Matches e.g. `wikipathways-20260910-gpml-Homo_sapiens.zip` inside href="..."
    for part in html.split('"') {
        let trimmed = part.trim();
        if trimmed.contains("gpml-Homo_sapiens.zip") && !trimmed.contains(".md5") {
            return Some(trimmed.to_string());
        }
    }
    None
}

fn extract_wp_id(filename: &str) -> Option<String> {
    for tok in filename.split(['_', '-', '.', '/']) {
        if tok.starts_with("WP") && tok.len() >= 3 && tok[2..].chars().all(|c| c.is_ascii_digit()) {
            return Some(tok.to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_wikipathways_directory_listing_and_wp_id() {
        let sample_html = r#"<a href="wikipathways-20261010-gpml-Homo_sapiens.zip">Download</a>"#;
        assert_eq!(
            extract_homo_sapiens_gpml_zip_name(sample_html).as_deref(),
            Some("wikipathways-20261010-gpml-Homo_sapiens.zip")
        );
        assert_eq!(
            extract_wp_id("Hs_TNF-alpha_NF-kB_Signaling_Pathway_WP231_12345.gpml").as_deref(),
            Some("WP231")
        );
    }
}
