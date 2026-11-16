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

use crate::svg_export::{render_pathway_map_svg, render_treemap_svg};
use crate::theme::{
    apply_workbench_visuals, draw_directed_causal_edge, draw_dna_double_helix, draw_entity_glyph,
    draw_phospho_badge, ScientificPalette, WorkbenchPaletteMode,
};
use ca_core::{
    build_mechanistic_network, builtin_cohorts, parse_expression_table, run_comparison_analysis,
    run_core_analysis, seed_kb::build_seed_knowledge_graph, simulate_map, ComparisonMetric,
    CoreAnalysisConfig, CoreAnalysisReport, EdgeSign, EntityKind, ExpressionDataset,
    KnowledgeGraph, PredictedState, SubcellularLocation,
};
use ca_streams::{
    client::{fetch_remote_stream_conditional, import_local_stream_file},
    default_stream_registry, fetch_and_cache_reactome_diagram, load_compressed_snapshot,
    run_batch_download_and_update, save_compressed_snapshot, BatchProgressEvent, BatchSyncSummary,
    SqliteCacheStats, SqliteKnowledgeStore, StreamDescriptor,
};
use egui::{Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Ui, Vec2};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc::{self, Receiver, Sender},
    Arc,
};
use std::thread;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkbenchTab {
    DatasetVolcano,
    CanonicalPathwaysMap,
    UpstreamRegulators,
    CausalNetworksMultiHop,
    DownstreamTreeMap,
    RegulatorEffectsFlow,
    ComparisonAndPerturb,
    OpenDataStreams,
}

impl WorkbenchTab {
    pub const ALL: [Self; 8] = [
        Self::DatasetVolcano,
        Self::CanonicalPathwaysMap,
        Self::UpstreamRegulators,
        Self::CausalNetworksMultiHop,
        Self::DownstreamTreeMap,
        Self::RegulatorEffectsFlow,
        Self::ComparisonAndPerturb,
        Self::OpenDataStreams,
    ];

    pub fn title(self) -> &'static str {
        match self {
            Self::DatasetVolcano => "1. Dataset & Volcano",
            Self::CanonicalPathwaysMap => "2. Canonical Pathways & MAP",
            Self::UpstreamRegulators => "3. Upstream Regulators",
            Self::CausalNetworksMultiHop => "4. Causal Networks (Multi-Hop)",
            Self::DownstreamTreeMap => "5. Downstream TreeMap",
            Self::RegulatorEffectsFlow => "6. Regulator Effects",
            Self::ComparisonAndPerturb => "7. Comparison & Reversers",
            Self::OpenDataStreams => "8. Streams & SQLite Cache",
        }
    }
}

enum WorkerMessage {
    BatchProgress(BatchProgressEvent),
    BatchFinished(BatchSyncSummary),
    StreamSynced {
        stream_idx: usize,
        result: Result<(usize, bool), String>,
    },
    LocalStreamImported {
        label: String,
        result: Result<(usize, String), String>,
    },
    ReactomeDiagramFetched {
        pathway_id: String,
        result: Result<usize, String>,
    },
}

enum GuiAction {
    SelectCohort(usize),
    RecomputeActive,
    RecomputeAllComparison,
    ImportDatasetDialog,
    ExportJsonReportDialog,
    ExportPathwaySvgDialog,
    ExportTreeMapSvgDialog,
    SaveSnapshotDialog,
    LoadSnapshotDialog,
    ImportStreamFileDialog,
    SyncRemoteStream { idx: usize, force: bool },
    StartBatchDownload { pathways: bool, regulatory: bool, force: bool },
    CancelBatchDownload,
    FetchReactomeDiagram(String),
    ClearAnalysisCache,
    VacuumSqlite,
    OpenStreamsTab,
    SelectPathway(usize),
    SelectPathwayNode(String),
    SetNodeClamp(String, Option<f64>),
    ClearAllClamps,
    SelectUraIndex(usize),
    SelectCnaIndex(usize),
    SelectDeaTerm(String),
    SelectRegEffectIndex(usize),
    SetComparisonMetric(ComparisonMetric),
    TogglePalette,
    TogglePathwayFullscreen,
    ResetPathwayViewport,
}

pub struct CausalAtlasApp {
    palette_mode: WorkbenchPaletteMode,
    active_tab: WorkbenchTab,
    db_path: PathBuf,
    db: SqliteKnowledgeStore,
    db_stats: SqliteCacheStats,
    analysis_cache_hits: usize,
    analysis_cache_misses: usize,
    kg: KnowledgeGraph,
    streams: Vec<StreamDescriptor>,
    cohorts: Vec<ExpressionDataset>,
    active_cohort_idx: usize,
    config: CoreAnalysisConfig,
    reports: Vec<CoreAnalysisReport>,

    gene_search: String,
    pathway_search: String,
    regulator_search: String,
    selected_pathway_idx: usize,
    selected_pathway_node: Option<String>,
    map_user_clamps: BTreeMap<String, f64>,
    selected_ura_idx: usize,
    mechanistic_node_offsets: BTreeMap<String, Vec2>,
    selected_cna_idx: usize,
    selected_dea_term: Option<String>,
    selected_reg_effect_idx: usize,
    comparison_metric: ComparisonMetric,
    pathway_fullscreen: bool,
    pathway_zoom: f32,
    pathway_pan: Vec2,

    batch_progress_fraction: f32,
    batch_progress_label: String,
    sync_logs: Vec<String>,
    cancel_flag: Arc<AtomicBool>,

    status_banner: String,
    worker_busy: bool,
    tx: Sender<WorkerMessage>,
    rx: Receiver<WorkerMessage>,
}

impl CausalAtlasApp {
    pub fn new(cc: &eframe::CreationContext<'_>, db_path: PathBuf) -> Self {
        Self::new_headless(&cc.egui_ctx, db_path)
    }

    pub fn new_headless(ctx: &egui::Context, db_path: PathBuf) -> Self {
        let palette_mode = WorkbenchPaletteMode::AcademicDark;
        apply_workbench_visuals(ctx, palette_mode);

        let mut db = SqliteKnowledgeStore::open(&db_path)
            .or_else(|_| SqliteKnowledgeStore::open_in_memory())
            .expect("SQLite store initialization must succeed");
        let _ = db.ensure_seeded(build_seed_knowledge_graph);

        let kg = db
            .load_knowledge_graph()
            .unwrap_or_else(|_| build_seed_knowledge_graph());
        let mut streams = default_stream_registry();
        for st in &mut streams {
            if let Ok(Some(meta)) = db.get_stream_meta(&st.id) {
                st.records_ingested = meta.records_count;
                st.last_sha256 = Some(meta.content_sha256.chars().take(12).collect());
                st.last_sync_status = format!("Cached in SQLite ({} records)", meta.records_count);
            }
        }

        let cohorts = builtin_cohorts();
        let config = CoreAnalysisConfig::default();
        let mut hits = 0usize;
        let mut misses = 0usize;

        let reports: Vec<CoreAnalysisReport> = cohorts
            .iter()
            .map(|ds| {
                if let Ok(Some(cached)) = db.get_cached_analysis(ds, &config) {
                    hits += 1;
                    cached
                } else {
                    misses += 1;
                    let rep = run_core_analysis(ds, &kg, &config);
                    let _ = db.put_cached_analysis(ds, &config, &rep);
                    rep
                }
            })
            .collect();

        let db_stats = db.stats().unwrap_or_default();
        let (tx, rx) = mpsc::channel();

        Self {
            palette_mode,
            active_tab: WorkbenchTab::CanonicalPathwaysMap,
            db_path,
            db,
            db_stats: db_stats.clone(),
            analysis_cache_hits: hits,
            analysis_cache_misses: misses,
            kg,
            streams,
            cohorts,
            active_cohort_idx: 0,
            config,
            reports,
            gene_search: String::new(),
            pathway_search: String::new(),
            regulator_search: String::new(),
            selected_pathway_idx: 0,
            selected_pathway_node: None,
            map_user_clamps: BTreeMap::new(),
            selected_ura_idx: 0,
            mechanistic_node_offsets: BTreeMap::new(),
            selected_cna_idx: 0,
            selected_dea_term: None,
            selected_reg_effect_idx: 0,
            comparison_metric: ComparisonMetric::UpstreamRegulatorsZ,
            pathway_fullscreen: false,
            pathway_zoom: 1.0,
            pathway_pan: Vec2::ZERO,
            batch_progress_fraction: 0.0,
            batch_progress_label: "Idle — SQLite WAL cache ready".to_string(),
            sync_logs: vec![format!(
                "Opened SQLite WAL store ({} pathways, {} signed edges, {} cached reports)",
                db_stats.cached_pathways, db_stats.cached_edges, db_stats.cached_analysis_reports
            )],
            cancel_flag: Arc::new(AtomicBool::new(false)),
            status_banner: format!(
                "SQLite WAL active ({}) — {} cache hits, {} computed",
                db_stats.db_path, hits, misses
            ),
            worker_busy: false,
            tx,
            rx,
        }
    }

    fn current_report(&self) -> &CoreAnalysisReport {
        &self.reports[self.active_cohort_idx]
    }

    fn refresh_from_sqlite(&mut self) {
        if let Ok(loaded) = self.db.load_knowledge_graph() {
            self.kg = loaded;
        }
        for st in &mut self.streams {
            if let Ok(Some(meta)) = self.db.get_stream_meta(&st.id) {
                st.records_ingested = meta.records_count;
                st.last_sha256 = Some(meta.content_sha256.chars().take(12).collect());
                st.last_sync_status = format!("Cached in SQLite ({} records)", meta.records_count);
            }
        }
        self.recompute_all();
        if let Ok(st) = self.db.stats() {
            self.db_stats = st;
        }
    }

    fn poll_worker(&mut self) {
        while let Ok(msg) = self.rx.try_recv() {
            match msg {
                WorkerMessage::BatchProgress(ev) => {
                    let denom = ev.total_steps.max(1) as f32;
                    self.batch_progress_fraction = (ev.step_index as f32 / denom).clamp(0.0, 1.0);
                    self.batch_progress_label = format!(
                        "[{}/{}] {} — {}",
                        ev.step_index, ev.total_steps, ev.source_name, ev.message
                    );
                    self.sync_logs.push(self.batch_progress_label.clone());
                }
                WorkerMessage::BatchFinished(sum) => {
                    self.worker_busy = false;
                    self.batch_progress_fraction = 1.0;
                    self.refresh_from_sqlite();
                    let msg = format!(
                        "Batch sync finished: {} checked, {} cache hits, {} downloaded (+{} pathways, +{} edges)",
                        sum.streams_checked,
                        sum.cache_hits,
                        sum.streams_downloaded,
                        sum.pathways_ingested,
                        sum.causal_edges_ingested
                    );
                    self.batch_progress_label = msg.clone();
                    self.status_banner = msg.clone();
                    self.sync_logs.push(msg);
                }
                WorkerMessage::StreamSynced { stream_idx, result } => {
                    self.worker_busy = false;
                    match result {
                        Ok((count, hit)) => {
                            self.refresh_from_sqlite();
                            if let Some(st) = self.streams.get_mut(stream_idx) {
                                st.records_ingested = count;
                                st.last_sync_status = if hit {
                                    format!("Up-to-date in SQLite ({count} records)")
                                } else {
                                    format!("Updated & cached ({count} records)")
                                };
                                let log = format!("{}: {}", st.name, st.last_sync_status);
                                self.status_banner = log.clone();
                                self.sync_logs.push(log);
                            }
                        }
                        Err(e) => {
                            if let Some(st) = self.streams.get_mut(stream_idx) {
                                st.last_sync_status = format!("Sync failed: {e}");
                            }
                            self.status_banner = format!("Stream sync error: {e}");
                            self.sync_logs.push(self.status_banner.clone());
                        }
                    }
                }
                WorkerMessage::LocalStreamImported { label, result } => {
                    self.worker_busy = false;
                    match result {
                        Ok((count, sha)) => {
                            self.refresh_from_sqlite();
                            let short_sha: String = sha.chars().take(12).collect();
                            let log = format!(
                                "Imported {count} records from '{label}' into SQLite (sha256:{short_sha})"
                            );
                            self.status_banner = log.clone();
                            self.sync_logs.push(log);
                        }
                        Err(e) => {
                            self.status_banner = format!("Local import failed: {e}");
                            self.sync_logs.push(self.status_banner.clone());
                        }
                    }
                }
                WorkerMessage::ReactomeDiagramFetched { pathway_id, result } => {
                    self.worker_busy = false;
                    match result {
                        Ok(coord_cnt) => {
                            if let Ok(loaded) = self.db.load_knowledge_graph() {
                                self.kg = loaded;
                            }
                            if let Ok(st) = self.db.stats() {
                                self.db_stats = st;
                            }
                            let log = format!(
                                "Cached Reactome 2D DiagramJS layout for {pathway_id} ({coord_cnt} positioned nodes)"
                            );
                            self.status_banner = log.clone();
                            self.sync_logs.push(log);
                        }
                        Err(e) => {
                            self.status_banner = format!("Diagram fetch note ({pathway_id}): {e}");
                        }
                    }
                }
            }
        }
    }

    fn compute_or_load_cached(&mut self, cohort_idx: usize) -> (CoreAnalysisReport, bool) {
        let ds = &self.cohorts[cohort_idx];
        if let Ok(Some(cached)) = self.db.get_cached_analysis(ds, &self.config) {
            self.analysis_cache_hits += 1;
            return (cached, true);
        }
        self.analysis_cache_misses += 1;
        let rep = run_core_analysis(ds, &self.kg, &self.config);
        let _ = self.db.put_cached_analysis(ds, &self.config, &rep);
        if let Ok(st) = self.db.stats() {
            self.db_stats = st;
        }
        (rep, false)
    }

    fn recompute_active(&mut self) {
        let idx = self.active_cohort_idx;
        let (rep, hit) = self.compute_or_load_cached(idx);
        self.reports[idx] = rep;
        self.status_banner = format!(
            "{} '{}' ({})",
            if hit { "SQLite cache hit (<1 ms) for" } else { "Computed & cached in SQLite:" },
            self.cohorts[idx].name,
            if hit {
                format!("{} hits / {} misses", self.analysis_cache_hits, self.analysis_cache_misses)
            } else {
                format!("rev {}", self.db_stats.kg_revision_hash)
            }
        );
    }

    fn recompute_all(&mut self) {
        let mut new_reports = Vec::with_capacity(self.cohorts.len());
        for idx in 0..self.cohorts.len() {
            let (rep, _) = self.compute_or_load_cached(idx);
            new_reports.push(rep);
        }
        self.reports = new_reports;
        self.selected_pathway_idx = 0;
        self.selected_ura_idx = 0;
        self.selected_cna_idx = 0;
        self.selected_reg_effect_idx = 0;
    }

    fn apply_action(&mut self, ctx: &egui::Context, action: GuiAction) {
        match action {
            GuiAction::OpenStreamsTab => {
                self.active_tab = WorkbenchTab::OpenDataStreams;
            }
            GuiAction::SelectCohort(idx) => {
                if idx < self.cohorts.len() {
                    self.active_cohort_idx = idx;
                    self.selected_pathway_idx = 0;
                    self.selected_ura_idx = 0;
                    self.selected_cna_idx = 0;
                    self.selected_reg_effect_idx = 0;
                    self.map_user_clamps.clear();
                    self.mechanistic_node_offsets.clear();
                }
            }
            GuiAction::RecomputeActive => {
                self.recompute_active();
            }
            GuiAction::RecomputeAllComparison => {
                self.recompute_all();
                self.status_banner = format!(
                    "Synchronized all {} cohorts via SQLite analysis cache",
                    self.cohorts.len()
                );
            }
            GuiAction::StartBatchDownload {
                pathways,
                regulatory,
                force,
            } => {
                if self.worker_busy {
                    return;
                }
                self.worker_busy = true;
                self.cancel_flag.store(false, Ordering::Relaxed);
                self.batch_progress_fraction = 0.02;
                self.batch_progress_label = "Initializing SQLite WAL batch downloader...".to_string();
                let db_path = self.db_path.clone();
                let tx = self.tx.clone();
                let cancel = Arc::clone(&self.cancel_flag);
                thread::spawn(move || {
                    if let Ok(mut bg_store) = SqliteKnowledgeStore::open(&db_path) {
                        let tx_prog = tx.clone();
                        let summary = run_batch_download_and_update(
                            &mut bg_store,
                            pathways,
                            regulatory,
                            force,
                            Some(&cancel),
                            move |ev| {
                                let _ = tx_prog.send(WorkerMessage::BatchProgress(ev));
                            },
                        );
                        let _ = tx.send(WorkerMessage::BatchFinished(summary));
                    } else {
                        let mut sum = BatchSyncSummary::default();
                        sum.errors
                            .push("Could not open background SQLite WAL connection".to_string());
                        let _ = tx.send(WorkerMessage::BatchFinished(sum));
                    }
                });
            }
            GuiAction::CancelBatchDownload => {
                self.cancel_flag.store(true, Ordering::Relaxed);
                self.status_banner = "Cancelling batch download after current request...".to_string();
            }
            GuiAction::FetchReactomeDiagram(pw_id) => {
                if self.worker_busy {
                    return;
                }
                if let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id).cloned() {
                    self.worker_busy = true;
                    self.status_banner =
                        format!("Fetching & caching Reactome 2D DiagramJS JSON for {pw_id}...");
                    let db_path = self.db_path.clone();
                    let tx = self.tx.clone();
                    thread::spawn(move || {
                        let res = SqliteKnowledgeStore::open(&db_path)
                            .and_then(|mut s| fetch_and_cache_reactome_diagram(&mut s, &pw_def))
                            .map(|u| u.node_coords.len())
                            .map_err(|e| e.to_string());
                        let _ = tx.send(WorkerMessage::ReactomeDiagramFetched {
                            pathway_id: pw_id,
                            result: res,
                        });
                    });
                }
            }
            GuiAction::ClearAnalysisCache => {
                match self.db.clear_analysis_cache() {
                    Ok(removed) => {
                        if let Ok(st) = self.db.stats() {
                            self.db_stats = st;
                        }
                        self.status_banner =
                            format!("Cleared {removed} cached CoreAnalysisReport entries from SQLite");
                        self.sync_logs.push(self.status_banner.clone());
                    }
                    Err(e) => self.status_banner = format!("Cache clear error: {e}"),
                }
            }
            GuiAction::VacuumSqlite => {
                match self.db.vacuum_and_optimize() {
                    Ok(()) => {
                        if let Ok(st) = self.db.stats() {
                            self.db_stats = st;
                        }
                        self.status_banner = format!(
                            "SQLite WAL checkpoint & VACUUM complete ({:.1} KB on disk)",
                            self.db_stats.file_size_bytes as f64 / 1024.0
                        );
                        self.sync_logs.push(self.status_banner.clone());
                    }
                    Err(e) => self.status_banner = format!("SQLite VACUUM error: {e}"),
                }
            }
            GuiAction::ImportDatasetDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Expression Table (TSV/CSV)", &["tsv", "csv", "txt"])
                    .pick_file()
                {
                    match std::fs::read_to_string(&path) {
                        Ok(raw) => {
                            let label = path
                                .file_stem()
                                .and_then(|s| s.to_str())
                                .unwrap_or("Custom Cohort");
                            match parse_expression_table(label, &raw) {
                                Ok(ds) => {
                                    self.cohorts.push(ds);
                                    let new_idx = self.cohorts.len() - 1;
                                    let (rep, _) = self.compute_or_load_cached(new_idx);
                                    self.reports.push(rep);
                                    self.active_cohort_idx = new_idx;
                                    self.status_banner =
                                        format!("Loaded & SQLite-cached dataset '{label}'");
                                }
                                Err(e) => {
                                    self.status_banner = format!("Dataset parse error: {e}");
                                }
                            }
                        }
                        Err(e) => {
                            self.status_banner = format!("Could not read file: {e}");
                        }
                    }
                }
            }
            GuiAction::ExportJsonReportDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .set_file_name("causal_atlas_report.json")
                    .save_file()
                {
                    match serde_json::to_string_pretty(self.current_report()) {
                        Ok(j) => match std::fs::write(&path, j) {
                            Ok(()) => {
                                self.status_banner =
                                    format!("Exported Core Analysis JSON to {}", path.display());
                            }
                            Err(e) => self.status_banner = format!("Write failed: {e}"),
                        },
                        Err(e) => self.status_banner = format!("JSON error: {e}"),
                    }
                }
            }
            GuiAction::ExportPathwaySvgDialog => {
                let rep = self.current_report();
                if let Some(pw_hit) = rep.canonical_pathways.get(self.selected_pathway_idx) {
                    if let Some(pw_def) = self
                        .kg
                        .pathways
                        .iter()
                        .find(|p| p.id == pw_hit.pathway_id)
                    {
                        let active = self.cohorts[self.active_cohort_idx].filter_active(
                            self.config.abs_log2fc_cutoff,
                            self.config.fdr_cutoff,
                            &self.kg,
                        );
                        let map_res = simulate_map(
                            pw_def,
                            &active,
                            &self.map_user_clamps,
                            self.config.map_damping,
                            self.config.map_iterations,
                        );
                        let svg = render_pathway_map_svg(pw_def, &map_res, &self.kg);
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name(format!("{}_map.svg", pw_def.id))
                            .save_file()
                        {
                            if let Err(e) = std::fs::write(&path, svg) {
                                self.status_banner = format!("SVG export error: {e}");
                            } else {
                                self.status_banner =
                                    format!("Saved vector pathway MAP SVG to {}", path.display());
                            }
                        }
                    }
                }
            }
            GuiAction::ExportTreeMapSvgDialog => {
                let rep = self.current_report();
                let svg = render_treemap_svg(&rep.dataset_name, &rep.treemap_rects);
                if let Some(path) = rfd::FileDialog::new()
                    .set_file_name("downstream_effects_treemap.svg")
                    .save_file()
                {
                    if let Err(e) = std::fs::write(&path, svg) {
                        self.status_banner = format!("SVG export error: {e}");
                    } else {
                        self.status_banner =
                            format!("Saved Downstream TreeMap SVG to {}", path.display());
                    }
                }
            }
            GuiAction::SaveSnapshotDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .set_file_name("knowledge_bundle.ca-kb.json.gz")
                    .save_file()
                {
                    match save_compressed_snapshot(&self.kg, &path) {
                        Ok(()) => {
                            self.status_banner = format!(
                                "Exported compressed knowledge snapshot ({} pathways) to {}",
                                self.kg.pathways.len(),
                                path.display()
                            );
                        }
                        Err(e) => self.status_banner = format!("Snapshot save error: {e}"),
                    }
                }
            }
            GuiAction::LoadSnapshotDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Knowledge Snapshot", &["gz", "json"])
                    .pick_file()
                {
                    match load_compressed_snapshot(&path) {
                        Ok(loaded) => {
                            let _ = self.db.upsert_knowledge_graph(&loaded);
                            self.refresh_from_sqlite();
                            self.status_banner =
                                format!("Merged snapshot into SQLite from {}", path.display());
                        }
                        Err(e) => self.status_banner = format!("Snapshot load error: {e}"),
                    }
                }
            }
            GuiAction::ImportStreamFileDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter(
                        "Open-Source Stream or Archive",
                        &["zip", "gmt", "gpml", "xml", "tsv", "txt", "gaf", "gz", "json"],
                    )
                    .pick_file()
                {
                    let tx = self.tx.clone();
                    let db_path = self.db_path.clone();
                    let label = path
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or("stream")
                        .to_string();
                    self.worker_busy = true;
                    self.status_banner = format!("Importing '{label}' into SQLite WAL store...");
                    thread::spawn(move || {
                        let res = import_local_stream_file(&path, None)
                            .and_then(|(delta, cnt, sha)| {
                                let mut s = SqliteKnowledgeStore::open(&db_path)?;
                                s.upsert_knowledge_graph(&delta)?;
                                Ok((cnt, sha))
                            })
                            .map_err(|e| e.to_string());
                        let _ = tx.send(WorkerMessage::LocalStreamImported { label, result: res });
                    });
                }
            }
            GuiAction::SyncRemoteStream { idx, force } => {
                if let Some(desc) = self.streams.get(idx).cloned() {
                    let tx = self.tx.clone();
                    let db_path = self.db_path.clone();
                    self.worker_busy = true;
                    self.status_banner = format!("Syncing {} into SQLite...", desc.name);
                    thread::spawn(move || {
                        let res = SqliteKnowledgeStore::open(&db_path)
                            .and_then(|mut s| fetch_remote_stream_conditional(&mut s, &desc, force))
                            .map_err(|e| e.to_string());
                        let _ = tx.send(WorkerMessage::StreamSynced {
                            stream_idx: idx,
                            result: res,
                        });
                    });
                }
            }
            GuiAction::SelectPathway(idx) => {
                self.selected_pathway_idx = idx;
                self.selected_pathway_node = None;
                self.pathway_zoom = 1.0;
                self.pathway_pan = Vec2::ZERO;
            }
            GuiAction::SelectPathwayNode(nid) => {
                self.selected_pathway_node = Some(nid);
            }
            GuiAction::SetNodeClamp(nid, val) => {
                if let Some(v) = val {
                    self.map_user_clamps.insert(nid, v);
                } else {
                    self.map_user_clamps.remove(&nid);
                }
            }
            GuiAction::ClearAllClamps => {
                self.map_user_clamps.clear();
            }
            GuiAction::SelectUraIndex(idx) => {
                self.selected_ura_idx = idx;
                self.mechanistic_node_offsets.clear();
            }
            GuiAction::SelectCnaIndex(idx) => {
                self.selected_cna_idx = idx;
            }
            GuiAction::SelectDeaTerm(tid) => {
                self.selected_dea_term = Some(tid);
            }
            GuiAction::SelectRegEffectIndex(idx) => {
                self.selected_reg_effect_idx = idx;
            }
            GuiAction::SetComparisonMetric(m) => {
                self.comparison_metric = m;
            }
            GuiAction::TogglePathwayFullscreen => {
                self.pathway_fullscreen = !self.pathway_fullscreen;
            }
            GuiAction::ResetPathwayViewport => {
                self.pathway_zoom = 1.0;
                self.pathway_pan = Vec2::ZERO;
            }
            GuiAction::TogglePalette => {
                self.palette_mode = match self.palette_mode {
                    WorkbenchPaletteMode::AcademicDark => WorkbenchPaletteMode::PublicationLight,
                    WorkbenchPaletteMode::PublicationLight => WorkbenchPaletteMode::AcademicDark,
                };
                apply_workbench_visuals(ctx, self.palette_mode);
            }
        }
    }
}

impl eframe::App for CausalAtlasApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.render_ui_frame(ctx);
    }
}

impl CausalAtlasApp {
    pub fn configure_huvec_tnf_screenshot_state(&mut self, ctx: &egui::Context) {
        self.configure_screenshot_view(ctx, 1);
    }

    pub fn configure_screenshot_view(&mut self, ctx: &egui::Context, view: u8) {
        self.apply_action(ctx, GuiAction::SelectCohort(0));
        self.apply_action(ctx, GuiAction::RecomputeActive);
        self.apply_action(ctx, GuiAction::RecomputeAllComparison);
        self.apply_action(
            ctx,
            GuiAction::SetComparisonMetric(ComparisonMetric::UpstreamRegulatorsZ),
        );
        self.active_tab = match view {
            2 => WorkbenchTab::CanonicalPathwaysMap,
            3 => WorkbenchTab::UpstreamRegulators,
            _ => WorkbenchTab::ComparisonAndPerturb,
        };
    }

    pub fn render_ui_frame(&mut self, ctx: &egui::Context) {
        self.poll_worker();
        if self.worker_busy {
            ctx.request_repaint();
        }

        let pal = ScientificPalette::for_mode(self.palette_mode);
        let mut actions: Vec<GuiAction> = Vec::new();

        egui::TopBottomPanel::top("top_nav").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    RichText::new("CAUSAL ATLAS")
                        .strong()
                        .size(16.0)
                        .color(pal.activated_strong),
                );
                ui.separator();
                for tab in WorkbenchTab::ALL {
                    let selected = self.active_tab == tab;
                    if ui.selectable_label(selected, tab.title()).clicked() {
                        self.active_tab = tab;
                    }
                }
                ui.separator();
                let theme_lbl = match self.palette_mode {
                    WorkbenchPaletteMode::AcademicDark => "Light Theme",
                    WorkbenchPaletteMode::PublicationLight => "Dark Theme",
                };
                if ui
                    .add_enabled(!self.worker_busy, egui::Button::new("Sync & Update DB"))
                    .clicked()
                {
                    actions.push(GuiAction::OpenStreamsTab);
                    actions.push(GuiAction::StartBatchDownload {
                        pathways: true,
                        regulatory: true,
                        force: false,
                    });
                }
                if ui.button("Export Report JSON").clicked() {
                    actions.push(GuiAction::ExportJsonReportDialog);
                }
                if ui.button(theme_lbl).clicked() {
                    actions.push(GuiAction::TogglePalette);
                }
            });
            ui.add_space(3.0);
        });

        egui::TopBottomPanel::bottom("bottom_status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(&self.status_banner).size(12.0));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        RichText::new(format!(
                            "SQLite WAL: {:.1} KB | {} pathways ({} w/ 2D coords) | {} edges | Cache Hits: {} / Misses: {}",
                            self.db_stats.file_size_bytes as f64 / 1024.0,
                            self.db_stats.cached_pathways,
                            self.db_stats.cached_diagrams_with_coords,
                            self.db_stats.cached_edges,
                            self.analysis_cache_hits,
                            self.analysis_cache_misses
                        ))
                        .size(11.5)
                        .color(pal.text_muted),
                    );
                });
            });
        });

        egui::SidePanel::left("left_controls")
            .resizable(true)
            .default_width(275.0)
            .min_width(235.0)
            .show(ctx, |ui| {
                ui.heading("Cohort & Cutoffs");
                ui.add_space(4.0);

                egui::ComboBox::from_label("Cohort")
                    .selected_text(&self.cohorts[self.active_cohort_idx].name)
                    .width(210.0)
                    .show_ui(ui, |ui| {
                        for (idx, ds) in self.cohorts.iter().enumerate() {
                            if ui
                                .selectable_label(self.active_cohort_idx == idx, &ds.name)
                                .clicked()
                            {
                                actions.push(GuiAction::SelectCohort(idx));
                            }
                        }
                    });

                ui.label(
                    RichText::new(&self.cohorts[self.active_cohort_idx].description)
                        .size(11.5)
                        .color(pal.text_muted),
                );
                ui.add_space(6.0);

                if ui.button("Import Expression CSV / TSV...").clicked() {
                    actions.push(GuiAction::ImportDatasetDialog);
                }

                ui.separator();
                ui.label(RichText::new("Differential Cutoffs").strong());
                let mut cutoff_changed = false;
                cutoff_changed |= ui
                    .add(
                        egui::Slider::new(&mut self.config.abs_log2fc_cutoff, 0.25..=3.0)
                            .text("|log2FC|"),
                    )
                    .changed();
                cutoff_changed |= ui
                    .add(
                        egui::Slider::new(&mut self.config.fdr_cutoff, 0.0001..=0.25)
                            .logarithmic(true)
                            .text("Max FDR"),
                    )
                    .changed();
                cutoff_changed |= ui
                    .add(
                        egui::Slider::new(&mut self.config.cna_max_depth, 1..=3)
                            .text("CNA Hop K"),
                    )
                    .changed();

                ui.horizontal(|ui| {
                    if ui.button("Run Core Analysis").clicked() || cutoff_changed {
                        actions.push(GuiAction::RecomputeActive);
                    }
                    if ui.button("Sync All Cohorts").clicked() {
                        actions.push(GuiAction::RecomputeAllComparison);
                    }
                });

                ui.separator();
                let sum = &self.reports[self.active_cohort_idx].summary;
                ui.label(RichText::new("Analysis Molecule Summary").strong());
                egui::Grid::new("kpi_grid")
                    .num_columns(2)
                    .spacing([12.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Measured Rows:");
                        ui.label(format!("{}", sum.total_measured));
                        ui.end_row();

                        ui.label("Mapped to KG:");
                        ui.label(format!("{}", sum.mapped_to_graph));
                        ui.end_row();

                        ui.label("Analysis-Ready DEGs:");
                        ui.label(RichText::new(format!("{}", sum.active_genes)).strong());
                        ui.end_row();

                        ui.label("Upregulated (+):");
                        ui.label(
                            RichText::new(format!("{}", sum.upregulated_count)).color(pal.expr_up),
                        );
                        ui.end_row();

                        ui.label("Downregulated (-):");
                        ui.label(
                            RichText::new(format!("{}", sum.downregulated_count))
                                .color(pal.expr_down),
                        );
                        ui.end_row();

                        ui.label("Sign Skew (mu_2):");
                        ui.label(format!("{:+.2}", sum.mean_sign_skew));
                        ui.end_row();
                    });

                ui.separator();
                ui.label(RichText::new("Visual Legend").strong().size(12.5));
                ui.horizontal(|ui| {
                    let (r, p) = ui.allocate_painter(Vec2::new(14.0, 14.0), Sense::hover());
                    p.rect_filled(r.rect, 3.0, pal.activated_strong);
                    ui.label(RichText::new("Predicted Activated (Z >= +2)").size(11.5));
                });
                ui.horizontal(|ui| {
                    let (r, p) = ui.allocate_painter(Vec2::new(14.0, 14.0), Sense::hover());
                    p.rect_filled(r.rect, 3.0, pal.inhibited_strong);
                    ui.label(RichText::new("Predicted Inhibited (Z <= -2)").size(11.5));
                });
                ui.horizontal(|ui| {
                    let (r, p) = ui.allocate_painter(Vec2::new(14.0, 14.0), Sense::hover());
                    p.rect_filled(r.rect, 3.0, pal.expr_up);
                    ui.label(RichText::new("Observed Upregulated mRNA").size(11.5));
                });
                ui.horizontal(|ui| {
                    let (r, p) = ui.allocate_painter(Vec2::new(14.0, 14.0), Sense::hover());
                    p.rect_filled(r.rect, 3.0, pal.expr_down);
                    ui.label(RichText::new("Observed Downregulated mRNA").size(11.5));
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| match self.active_tab {
            WorkbenchTab::DatasetVolcano => {
                self.render_dataset_volcano_tab(ui, &pal);
            }
            WorkbenchTab::CanonicalPathwaysMap => {
                self.render_canonical_map_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::UpstreamRegulators => {
                self.render_ura_mechanistic_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::CausalNetworksMultiHop => {
                self.render_cna_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::DownstreamTreeMap => {
                self.render_dea_treemap_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::RegulatorEffectsFlow => {
                self.render_regulator_effects_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::ComparisonAndPerturb => {
                self.render_comparison_perturb_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::OpenDataStreams => {
                self.render_open_streams_tab(ui, &pal, &mut actions);
            }
        });

        for act in actions {
            self.apply_action(ctx, act);
        }
    }
}

impl CausalAtlasApp {
    fn render_dataset_volcano_tab(&mut self, ui: &mut Ui, pal: &ScientificPalette) {
        let Self {
            cohorts,
            active_cohort_idx,
            config,
            gene_search,
            ..
        } = self;

        let ds = &cohorts[*active_cohort_idx];
        let fc_cut = config.abs_log2fc_cutoff;
        let q_cut = config.fdr_cutoff;

        ui.columns(2, |cols| {
            cols[0].heading("Interactive Volcano Plot (log2FC vs. -log10 FDR)");
            let avail = cols[0].available_size();
            let plot_size = Vec2::new(avail.x.max(280.0), (avail.y - 32.0).max(280.0));
            let (resp, painter) = cols[0].allocate_painter(plot_size, Sense::hover());
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, pal.bg_canvas);
            painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

            let max_x = ds
                .rows
                .iter()
                .map(|r| r.log2_fc.abs())
                .fold(3.5f64, f64::max)
                * 1.1;
            let max_y = ds
                .rows
                .iter()
                .map(|r| -r.fdr.max(1e-25).log10())
                .fold(6.0f64, f64::max)
                * 1.08;

            let pad = 36.0f32;
            let inner = Rect::from_min_max(
                Pos2::new(rect.left() + pad, rect.top() + 18.0),
                Pos2::new(rect.right() - 18.0, rect.bottom() - pad),
            );

            let map_pt = |fc: f64, neg_log_q: f64| -> Pos2 {
                let nx = ((fc + max_x) / (2.0 * max_x)).clamp(0.0, 1.0) as f32;
                let ny = (neg_log_q / max_y).clamp(0.0, 1.0) as f32;
                Pos2::new(
                    inner.left() + nx * inner.width(),
                    inner.bottom() - ny * inner.height(),
                )
            };

            let x_pos = map_pt(fc_cut, 0.0).x;
            let x_neg = map_pt(-fc_cut, 0.0).x;
            let y_q = map_pt(0.0, -q_cut.max(1e-25).log10()).y;
            let guide_stroke = Stroke::new(1.0, pal.grid_line);
            painter.line_segment(
                [Pos2::new(x_pos, inner.top()), Pos2::new(x_pos, inner.bottom())],
                guide_stroke,
            );
            painter.line_segment(
                [Pos2::new(x_neg, inner.top()), Pos2::new(x_neg, inner.bottom())],
                guide_stroke,
            );
            painter.line_segment(
                [Pos2::new(inner.left(), y_q), Pos2::new(inner.right(), y_q)],
                guide_stroke,
            );

            let hover_pos = resp.hover_pos();
            let mut hovered_gene: Option<(&str, f64, f64)> = None;

            for row in &ds.rows {
                let yv = -row.fdr.max(1e-25).log10();
                let pt = map_pt(row.log2_fc, yv);
                let active = row.log2_fc.abs() >= fc_cut && row.fdr <= q_cut;
                let col = if !active {
                    pal.text_muted
                } else if row.log2_fc > 0.0 {
                    pal.expr_up
                } else {
                    pal.expr_down
                };
                let rad = if active { 4.8 } else { 3.0 };
                painter.circle_filled(pt, rad, col);
                if active && yv > max_y * 0.35 {
                    painter.text(
                        Pos2::new(pt.x, pt.y - 9.0),
                        Align2::CENTER_BOTTOM,
                        &row.symbol,
                        FontId::proportional(10.5),
                        pal.text_primary,
                    );
                }
                if let Some(hp) = hover_pos {
                    if hp.distance(pt) <= 8.0 {
                        hovered_gene = Some((&row.symbol, row.log2_fc, row.fdr));
                    }
                }
            }

            if let Some((sym, fc, q)) = hovered_gene {
                resp.on_hover_text(format!("{sym}\nlog2FC: {fc:+.3}\nFDR: {q:.2e}"));
            }

            cols[1].heading("Measured Gene Table");
            cols[1].horizontal(|ui| {
                ui.label("Filter symbol:");
                ui.text_edit_singleline(gene_search);
            });
            let needle = gene_search.trim().to_uppercase();
            cols[1].push_id("dataset_gene_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("gene_table_grid")
                        .striped(true)
                        .num_columns(4)
                        .spacing([14.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Symbol");
                            ui.strong("log2FC");
                            ui.strong("P-value");
                            ui.strong("BH-FDR");
                            ui.end_row();

                            for r in ds
                                .rows
                                .iter()
                                .filter(|r| needle.is_empty() || r.symbol.contains(&needle))
                            {
                                let pass = r.log2_fc.abs() >= fc_cut && r.fdr <= q_cut;
                                let c = if !pass {
                                    pal.text_muted
                                } else if r.log2_fc > 0.0 {
                                    pal.expr_up
                                } else {
                                    pal.expr_down
                                };
                                ui.colored_label(c, &r.symbol);
                                ui.colored_label(c, format!("{:+.3}", r.log2_fc));
                                ui.label(format!("{:.2e}", r.p_value));
                                ui.label(format!("{:.2e}", r.fdr));
                                ui.end_row();
                            }
                        });
                });
            });
        });
    }

    fn render_canonical_map_tab(
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let Self {
            reports,
            cohorts,
            active_cohort_idx,
            config,
            kg,
            pathway_search,
            selected_pathway_idx,
            selected_pathway_node,
            map_user_clamps,
            pathway_fullscreen,
            pathway_zoom,
            pathway_pan,
            ..
        } = self;

        let rep = &reports[*active_cohort_idx];
        let sel_pw_idx = *selected_pathway_idx;
        let is_full = *pathway_fullscreen;

        let render_diagram_pane = |pane_ui: &mut Ui,
                                   actions: &mut Vec<GuiAction>,
                                   pathway_zoom: &mut f32,
                                   pathway_pan: &mut Vec2| {
            pane_ui.horizontal_wrapped(|ui| {
                ui.heading("Canonical Pathway Diagram & MAP");
                let fs_lbl = if is_full { "Exit Full Screen" } else { "Full Screen Map" };
                if ui.button(fs_lbl).clicked() {
                    actions.push(GuiAction::TogglePathwayFullscreen);
                }
                if ui.button(format!("Reset View ({:.0}%)", *pathway_zoom * 100.0)).clicked() {
                    actions.push(GuiAction::ResetPathwayViewport);
                }
                if ui.button("Export SVG").clicked() {
                    actions.push(GuiAction::ExportPathwaySvgDialog);
                }
                if let Some(sh) = rep.canonical_pathways.get(sel_pw_idx) {
                    if sh.pathway_id.starts_with("R-HSA-") {
                        let needs_coords = kg
                            .pathways
                            .iter()
                            .find(|p| p.id == sh.pathway_id)
                            .is_some_and(|p| p.node_coords.is_empty());
                        if needs_coords && ui.button("Fetch Reactome 2D Layout").clicked() {
                            actions.push(GuiAction::FetchReactomeDiagram(sh.pathway_id.clone()));
                        }
                    }
                }
                if !map_user_clamps.is_empty() && ui.button("Reset MAP Clamps").clicked() {
                    actions.push(GuiAction::ClearAllClamps);
                }
            });

            if let Some(sel_hit) = rep.canonical_pathways.get(sel_pw_idx) {
                if let Some(pw_def) = kg.pathways.iter().find(|p| p.id == sel_hit.pathway_id) {
                    let active = cohorts[*active_cohort_idx].filter_active(
                        config.abs_log2fc_cutoff,
                        config.fdr_cutoff,
                        kg,
                    );
                    let map_res = simulate_map(
                        pw_def,
                        &active,
                        map_user_clamps,
                        config.map_damping,
                        config.map_iterations,
                    );

                    if let Some(ref sel_nid) = *selected_pathway_node {
                        pane_ui.horizontal(|ui| {
                            ui.label(RichText::new(format!("Node: {sel_nid}")).strong());
                            if ui.button("Activate (+1)").clicked() {
                                actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), Some(1.0)));
                            }
                            if ui.button("Inhibit (-1)").clicked() {
                                actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), Some(-1.0)));
                            }
                            if ui.button("Unclamp").clicked() {
                                actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), None));
                            }
                        });
                    } else {
                        pane_ui.label(
                            RichText::new(
                                "Drag to pan, scroll wheel to zoom, or click any molecule to clamp in-silico MAP activity.",
                            )
                            .size(11.5)
                            .color(pal.text_muted),
                        );
                    }

                    let canvas_sz = Vec2::new(
                        pane_ui.available_width().max(320.0),
                        (pane_ui.available_height() - 8.0).max(340.0),
                    );
                    let (resp, painter) =
                        pane_ui.allocate_painter(canvas_sz, Sense::click_and_drag());
                    let rect = resp.rect;
                    let painter = painter.with_clip_rect(rect);
                    painter.rect_filled(rect, 6.0, pal.bg_canvas);
                    painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

                    if resp.dragged() {
                        *pathway_pan += resp.drag_delta();
                    }
                    if resp.hovered() {
                        let scroll_y = pane_ui.input(|i| i.smooth_scroll_delta.y);
                        if scroll_y.abs() > 0.5 {
                            let factor = if scroll_y > 0.0 { 1.08 } else { 0.92 };
                            *pathway_zoom = (*pathway_zoom * factor).clamp(0.55, 2.4);
                        }
                    }

                    let zoom = *pathway_zoom;
                    let pan = *pathway_pan;
                    let vw = rect.width().max(320.0);
                    let vh = rect.height().max(320.0);
                    let center = rect.center();

                    let to_screen = |nx: f32, ny: f32| -> Pos2 {
                        let wx = (nx - 0.5) * (vw - 96.0).max(200.0);
                        let wy = (ny - 0.5) * (vh - 56.0).max(200.0);
                        Pos2::new(center.x + wx * zoom + pan.x, center.y + wy * zoom + pan.y)
                    };

                    // 1. Solid architectural Plasma Membrane (y=0.18) & Nuclear Envelope (y=0.72) bars
                    let pm_l = to_screen(0.02, 0.18);
                    let pm_r = to_screen(0.98, 0.18);
                    painter.line_segment([pm_l, pm_r], Stroke::new(6.5 * zoom, pal.border_default));

                    let nm_l = to_screen(0.06, 0.72);
                    let nm_r = to_screen(0.94, 0.72);
                    painter.line_segment([nm_l, nm_r], Stroke::new(6.0 * zoom, pal.border_default));

                    painter.text(
                        Pos2::new(rect.left() + 12.0, rect.top() + 8.0),
                        Align2::LEFT_TOP,
                        SubcellularLocation::Extracellular.label(),
                        FontId::proportional(11.0 * zoom.clamp(0.8, 1.3)),
                        pal.text_muted,
                    );
                    painter.text(
                        Pos2::new(pm_l.x + 4.0, pm_l.y + 22.0 * zoom),
                        Align2::LEFT_TOP,
                        SubcellularLocation::Cytoplasm.label(),
                        FontId::proportional(11.0 * zoom.clamp(0.8, 1.3)),
                        pal.text_muted,
                    );
                    if !pw_def.cell_context_label.is_empty() {
                        painter.text(
                            Pos2::new(rect.right() - 12.0, rect.top() + 8.0),
                            Align2::RIGHT_TOP,
                            &pw_def.cell_context_label,
                            FontId::proportional(11.0 * zoom.clamp(0.8, 1.3)),
                            pal.text_muted,
                        );
                    }
                    painter.text(
                        Pos2::new(nm_l.x + 4.0, nm_l.y + 11.0 * zoom),
                        Align2::LEFT_TOP,
                        SubcellularLocation::Nucleus.label(),
                        FontId::proportional(11.0 * zoom.clamp(0.8, 1.3)),
                        pal.text_muted,
                    );

                    // 2. Nuclear DNA Double Helix across promoter band (y=0.855)
                    let dna_l = to_screen(0.11, 0.855);
                    let dna_r = to_screen(0.89, 0.855);
                    draw_dna_double_helix(
                        &painter,
                        dna_l.x,
                        dna_r.x,
                        dna_l.y,
                        6.0 * zoom,
                        pal.grid_line,
                    );

                    // 3. Resolve vertical-lane or GPML coordinates
                    let norm_coords = pw_def.resolve_layout_coords(kg);
                    let positions: BTreeMap<String, Pos2> = norm_coords
                        .into_iter()
                        .map(|(k, (nx, ny))| (k, to_screen(nx, ny)))
                        .collect();

                    // 4. Draw directed pathway edges
                    for (e_idx, edge) in pw_def.topology_edges.iter().enumerate() {
                        if let (Some(&p1), Some(&p2)) =
                            (positions.get(&edge.from), positions.get(&edge.to))
                        {
                            let (consistent, flux) =
                                map_res.edge_flux.get(e_idx).copied().unwrap_or((true, 0.0));
                            let edge_col = if !consistent {
                                pal.inconsistent_amber
                            } else if flux > 0.12 {
                                pal.activated_strong
                            } else if flux < -0.12 {
                                pal.inhibited_strong
                            } else {
                                pal.border_default
                            };
                            draw_directed_causal_edge(
                                &painter,
                                p1,
                                p2,
                                22.0 * zoom,
                                edge.sign,
                                Stroke::new(1.65 * zoom.clamp(0.8, 1.4), edge_col),
                            );
                        }
                    }

                    // 5. Draw classic entity glyphs & (P) phosphorylation badges
                    let click_pt = if resp.clicked() {
                        resp.interact_pointer_pos()
                    } else {
                        None
                    };
                    let hover_pt = resp.hover_pos();

                    for (nid, st) in &map_res.states {
                        if let Some(&pos) = positions.get(nid) {
                            let kind = kg
                                .nodes
                                .get(nid)
                                .map(|n| n.kind)
                                .unwrap_or(EntityKind::GeneOther);
                            let fill = if let Some(fc) = st.observed_log2fc {
                                if st.user_clamp.is_some() {
                                    pal.color_for_activity(st.predicted_activity)
                                } else {
                                    pal.color_for_log2fc(fc)
                                }
                            } else {
                                pal.color_for_activity(st.predicted_activity)
                            };

                            let is_sel = selected_pathway_node.as_deref() == Some(nid.as_str());
                            let border = if is_sel || st.user_clamp.is_some() {
                                Stroke::new(2.5, Color32::from_rgb(168, 85, 247))
                            } else if st.observed_log2fc.is_some() {
                                Stroke::new(1.8, pal.text_primary)
                            } else {
                                Stroke::new(1.2, pal.border_default)
                            };

                            draw_entity_glyph(
                                &painter,
                                pos,
                                Vec2::new(29.0 * zoom, 15.0 * zoom),
                                kind,
                                fill,
                                border,
                            );
                            if pw_def.phospho_nodes.contains(nid) {
                                draw_phospho_badge(
                                    &painter,
                                    Pos2::new(pos.x + 23.0 * zoom, pos.y - 13.0 * zoom),
                                    pal.bg_canvas,
                                    pal.border_default,
                                );
                            }
                            painter.text(
                                Pos2::new(
                                    pos.x,
                                    pos.y + if kind == EntityKind::Kinase { -2.0 * zoom } else { 0.0 },
                                ),
                                Align2::CENTER_CENTER,
                                nid,
                                FontId::proportional(10.0 * zoom.clamp(0.85, 1.35)),
                                pal.text_primary,
                            );

                            if let Some(cp) = click_pt {
                                if cp.distance(pos) <= 26.0 * zoom {
                                    actions.push(GuiAction::SelectPathwayNode(nid.clone()));
                                }
                            }
                            if let Some(hp) = hover_pt {
                                if hp.distance(pos) <= 24.0 * zoom {
                                    let obs_str = st
                                        .observed_log2fc
                                        .map(|v| format!("{v:+.2}"))
                                        .unwrap_or_else(|| "unmeasured".to_string());
                                    resp.clone().on_hover_text(format!(
                                        "{nid} ({})
Observed log2FC: {obs_str}
MAP Activity: {:+.2}",
                                        kind.label(),
                                        st.predicted_activity
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        };

        if is_full {
            render_diagram_pane(ui, actions, pathway_zoom, pathway_pan);
            return;
        }

        ui.columns(2, |cols| {
            cols[0].horizontal(|ui| {
                ui.heading("Canonical Pathways");
                ui.add_space(8.0);
                ui.text_edit_singleline(pathway_search);
            });

            let needle = pathway_search.trim().to_lowercase();
            cols[0].push_id("canonical_pw_scroll", |ui| {
                egui::ScrollArea::vertical()
                    .max_height(260.0)
                    .show(ui, |ui| {
                        let max_bar = rep
                            .canonical_pathways
                            .iter()
                            .map(|p| p.neg_log10_fdr)
                            .fold(3.0f64, f64::max);

                        for (idx, hit) in rep.canonical_pathways.iter().enumerate() {
                            if !needle.is_empty()
                                && !hit.pathway_name.to_lowercase().contains(&needle)
                                && !hit.pathway_id.to_lowercase().contains(&needle)
                            {
                                continue;
                            }
                            let is_sel = sel_pw_idx == idx;
                            ui.horizontal_wrapped(|ui| {
                                if ui
                                    .selectable_label(
                                        is_sel,
                                        format!(
                                            "{} ({}/{})",
                                            hit.pathway_name, hit.overlap_count, hit.pathway_size
                                        ),
                                    )
                                    .clicked()
                                {
                                    actions.push(GuiAction::SelectPathway(idx));
                                }
                            });
                            let bar_w = ((hit.neg_log10_fdr / max_bar).clamp(0.04, 1.0) as f32)
                                * (ui.available_width() - 168.0).max(60.0);
                            ui.horizontal(|ui| {
                                let (r, p) =
                                    ui.allocate_painter(Vec2::new(bar_w, 12.0), Sense::click());
                                let fill = hit
                                    .activation_z
                                    .map(|z| pal.color_for_z(z))
                                    .unwrap_or(pal.neutral_node);
                                p.rect_filled(r.rect, 3.0, fill);
                                if r.clicked() {
                                    actions.push(GuiAction::SelectPathway(idx));
                                }
                                let z_txt = hit
                                    .activation_z
                                    .map(|z| format!("Z={z:+.2}"))
                                    .unwrap_or_else(|| "Z=n/a".to_string());
                                ui.label(
                                    RichText::new(format!(
                                        "-log10(FDR)={:.1} | {z_txt}",
                                        hit.neg_log10_fdr
                                    ))
                                    .size(11.0)
                                    .color(pal.text_muted),
                                );
                            });
                            ui.add_space(3.0);
                        }
                    });
            });

            cols[0].separator();
            if let Some(sel_hit) = rep.canonical_pathways.get(sel_pw_idx) {
                cols[0].label(
                    RichText::new(format!(
                        "{} [{}] — Source: {}",
                        sel_hit.pathway_name, sel_hit.pathway_id, sel_hit.source_db
                    ))
                    .strong(),
                );
                cols[0].label(format!(
                    "Overlap: {}/{} ({:.1}%) | P = {:.2e} | BH-FDR = {:.2e}",
                    sel_hit.overlap_count,
                    sel_hit.pathway_size,
                    sel_hit.overlap_ratio * 100.0,
                    sel_hit.p_value,
                    sel_hit.fdr
                ));
                if !sel_hit.upregulated_genes.is_empty() {
                    cols[0].colored_label(
                        pal.expr_up,
                        format!("Up DEGs: {}", sel_hit.upregulated_genes.join(", ")),
                    );
                }
                if !sel_hit.downregulated_genes.is_empty() {
                    cols[0].colored_label(
                        pal.expr_down,
                        format!("Down DEGs: {}", sel_hit.downregulated_genes.join(", ")),
                    );
                }
            }

            render_diagram_pane(&mut cols[1], actions, pathway_zoom, pathway_pan);
        });
    }

    fn render_ura_mechanistic_tab(
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let Self {
            reports,
            cohorts,
            active_cohort_idx,
            config,
            kg,
            regulator_search,
            selected_ura_idx,
            mechanistic_node_offsets,
            ..
        } = self;

        let rep = &reports[*active_cohort_idx];
        let sel_u_idx = *selected_ura_idx;

        ui.columns(2, |cols| {
            cols[0].horizontal(|ui| {
                ui.heading("Upstream Regulators (URA)");
                ui.text_edit_singleline(regulator_search);
            });
            let needle = regulator_search.trim().to_uppercase();

            cols[0].push_id("ura_table_scroll", |ui| {
                egui::ScrollArea::vertical()
                    .max_height(310.0)
                    .show(ui, |ui| {
                        egui::Grid::new("ura_grid")
                            .striped(true)
                            .num_columns(7)
                            .spacing([10.0, 5.0])
                            .show(ui, |ui| {
                                ui.strong("Regulator");
                                ui.strong("Type");
                                ui.strong("State");
                                ui.strong("Act Z");
                                ui.strong("Bias Z");
                                ui.strong("Overlap P");
                                ui.strong("Targets");
                                ui.end_row();

                                for (idx, h) in rep.upstream_regulators.iter().enumerate() {
                                    if !needle.is_empty() && !h.regulator_id.contains(&needle) {
                                        continue;
                                    }
                                    let sel = sel_u_idx == idx;
                                    if ui.selectable_label(sel, &h.regulator_id).clicked() {
                                        actions.push(GuiAction::SelectUraIndex(idx));
                                    }
                                    ui.label(h.kind.label());
                                    let st_col = match h.predicted_state {
                                        PredictedState::Activated => pal.activated_strong,
                                        PredictedState::Inhibited => pal.inhibited_strong,
                                        PredictedState::Neutral => pal.text_muted,
                                    };
                                    ui.colored_label(st_col, h.predicted_state.label());
                                    let skew_mark = if h.skew_flagged { "*" } else { "" };
                                    ui.colored_label(
                                        st_col,
                                        format!("{:+.2}{skew_mark}", h.activation_z),
                                    );
                                    ui.label(format!("{:+.2}", h.bias_corrected_z));
                                    ui.label(format!("{:.2e}", h.p_value));
                                    ui.label(format!("{}/{}", h.overlap_count, h.regulon_size));
                                    ui.end_row();
                                }
                            });
                    });
            });

            if let Some(sel_u) = rep.upstream_regulators.get(sel_u_idx) {
                cols[0].separator();
                cols[0].label(
                    RichText::new(format!(
                        "Target Evidence for {} ({})",
                        sel_u.regulator_name,
                        sel_u.mechanistic_summary.as_deref().unwrap_or("Direct regulon")
                    ))
                    .strong(),
                );
                cols[0].push_id("ura_targets_scroll", |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        egui::Grid::new("ura_ev_grid")
                            .striped(true)
                            .num_columns(5)
                            .spacing([12.0, 4.0])
                            .show(ui, |ui| {
                                ui.strong("Target");
                                ui.strong("Known Edge");
                                ui.strong("Dataset log2FC");
                                ui.strong("Match");
                                ui.strong("Stream");
                                ui.end_row();
                                for ev in &sel_u.target_evidences {
                                    ui.label(&ev.target_gene);
                                    ui.label(match ev.edge_sign {
                                        EdgeSign::Activating => "Activates (+1)",
                                        EdgeSign::Inhibiting => "Inhibits (-1)",
                                        EdgeSign::Ambiguous => "Binding (0)",
                                    });
                                    let c = if ev.observed_log2fc > 0.0 {
                                        pal.expr_up
                                    } else {
                                        pal.expr_down
                                    };
                                    ui.colored_label(c, format!("{:+.2}", ev.observed_log2fc));
                                    ui.label(if ev.consistent {
                                        "Consistent (+)"
                                    } else {
                                        "Opposite (-)"
                                    });
                                    ui.label(&ev.source_db);
                                    ui.end_row();
                                }
                            });
                    });
                });

                cols[1].heading(format!("Mechanistic Network — {}", sel_u.regulator_id));
                cols[1].label(
                    RichText::new(
                        "Transitive causal cascade connecting upstream co-regulators (drag nodes to rearrange).",
                    )
                    .size(11.5)
                    .color(pal.text_muted),
                );

                let active = cohorts[*active_cohort_idx].filter_active(
                    config.abs_log2fc_cutoff,
                    config.fdr_cutoff,
                    kg,
                );
                let mn = build_mechanistic_network(
                    &sel_u.regulator_id,
                    &rep.upstream_regulators,
                    &active,
                    kg,
                    config.mechanistic_breadth,
                    config.mechanistic_depth,
                );

                let sz = Vec2::new(
                    cols[1].available_width().max(280.0),
                    (cols[1].available_height() - 8.0).max(320.0),
                );
                let (resp, painter) = cols[1].allocate_painter(sz, Sense::click_and_drag());
                let rect = resp.rect;
                painter.rect_filled(rect, 6.0, pal.bg_canvas);
                painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

                let center = rect.center();
                let mut node_pos: BTreeMap<String, Pos2> = BTreeMap::new();

                let root_base = Pos2::new(center.x, rect.top() + 58.0);
                let root_off = mechanistic_node_offsets
                    .get(&mn.root_regulator)
                    .copied()
                    .unwrap_or(Vec2::ZERO);
                node_pos.insert(mn.root_regulator.clone(), root_base + root_off);

                let other_regs: Vec<&String> = mn
                    .regulator_ids
                    .iter()
                    .filter(|r| *r != &mn.root_regulator)
                    .collect();
                for (i, r) in other_regs.iter().enumerate() {
                    let frac = (i as f32 + 0.5) / (other_regs.len().max(1) as f32);
                    let base = Pos2::new(
                        rect.left() + 65.0 + frac * (rect.width() - 130.0),
                        rect.top() + rect.height() * 0.42,
                    );
                    let off = mechanistic_node_offsets
                        .get(*r)
                        .copied()
                        .unwrap_or(Vec2::ZERO);
                    node_pos.insert((*r).clone(), base + off);
                }

                let shown_targets: Vec<TargetGeneRef<'_>> = sel_u
                    .target_evidences
                    .iter()
                    .take(10)
                    .map(|t| TargetGeneRef {
                        sym: &t.target_gene,
                        fc: t.observed_log2fc,
                        sign: t.edge_sign,
                    })
                    .collect();

                for (i, tg) in shown_targets.iter().enumerate() {
                    let frac = (i as f32 + 0.5) / (shown_targets.len().max(1) as f32);
                    let base = Pos2::new(
                        rect.left() + 45.0 + frac * (rect.width() - 90.0),
                        rect.bottom() - 56.0,
                    );
                    let off = mechanistic_node_offsets
                        .get(tg.sym)
                        .copied()
                        .unwrap_or(Vec2::ZERO);
                    node_pos.insert(tg.sym.to_string(), base + off);
                }

                if resp.dragged() {
                    if let Some(ptr) = resp.interact_pointer_pos() {
                        let delta = resp.drag_delta();
                        if let Some((best_id, _)) =
                            node_pos.iter().find(|(_, &p)| p.distance(ptr) <= 32.0)
                        {
                            let entry = mechanistic_node_offsets
                                .entry(best_id.clone())
                                .or_insert(Vec2::ZERO);
                            *entry += delta;
                        }
                    }
                }

                for e in &mn.edges {
                    if let (Some(&p1), Some(&p2)) = (
                        node_pos.get(&e.from_regulator),
                        node_pos.get(&e.to_regulator),
                    ) {
                        let c = match e.sign {
                            EdgeSign::Activating => pal.activated_strong,
                            EdgeSign::Inhibiting => pal.inhibited_strong,
                            EdgeSign::Ambiguous => pal.border_default,
                        };
                        draw_directed_causal_edge(
                            &painter,
                            p1,
                            p2,
                            24.0,
                            e.sign,
                            Stroke::new(2.2, c),
                        );
                    }
                }

                if let Some(&rp) = node_pos.get(&mn.root_regulator) {
                    for tg in &shown_targets {
                        if let Some(&tp) = node_pos.get(tg.sym) {
                            draw_directed_causal_edge(
                                &painter,
                                rp,
                                tp,
                                20.0,
                                tg.sign,
                                Stroke::new(1.1, pal.grid_line),
                            );
                        }
                    }
                }

                for rid in &mn.regulator_ids {
                    if let Some(&p) = node_pos.get(rid) {
                        let z = mn.regulator_z_scores.get(rid).copied().unwrap_or(0.0);
                        let kind = kg
                            .nodes
                            .get(rid)
                            .map(|n| n.kind)
                            .unwrap_or(EntityKind::TranscriptionRegulator);
                        draw_entity_glyph(
                            &painter,
                            p,
                            Vec2::new(32.0, 16.0),
                            kind,
                            pal.color_for_z(z),
                            Stroke::new(2.0, pal.text_primary),
                        );
                        painter.text(
                            p,
                            Align2::CENTER_CENTER,
                            rid,
                            FontId::proportional(11.0),
                            pal.text_primary,
                        );
                    }
                }

                for tg in &shown_targets {
                    if let Some(&p) = node_pos.get(tg.sym) {
                        draw_entity_glyph(
                            &painter,
                            p,
                            Vec2::new(25.0, 13.0),
                            EntityKind::GeneOther,
                            pal.color_for_log2fc(tg.fc),
                            Stroke::new(1.3, pal.border_default),
                        );
                        painter.text(
                            p,
                            Align2::CENTER_CENTER,
                            tg.sym,
                            FontId::proportional(10.0),
                            pal.text_primary,
                        );
                    }
                }
            }
        });
    }

    fn render_cna_tab(
        &self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = &self.reports[self.active_cohort_idx];
        ui.columns(2, |cols| {
            cols[0].heading("Multi-Hop Causal Networks (CNA, Depth K=1..3)");
            cols[0].label(
                RichText::new(
                    "Identifies root master regulators (even when transcriptionally unchanged) via multi-step signaling cascades and degree-binned Monte Carlo p-values.",
                )
                .size(11.5)
                .color(pal.text_muted),
            );

            cols[0].push_id("cna_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("cna_grid")
                        .striped(true)
                        .num_columns(6)
                        .spacing([10.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Master Root");
                            ui.strong("Depth");
                            ui.strong("State");
                            ui.strong("Act Z");
                            ui.strong("Overlap P");
                            ui.strong("Net-Bias P");
                            ui.end_row();

                            for (idx, h) in rep.causal_networks.iter().enumerate() {
                                let sel = self.selected_cna_idx == idx;
                                if ui.selectable_label(sel, &h.root_regulator).clicked() {
                                    actions.push(GuiAction::SelectCnaIndex(idx));
                                }
                                ui.label(format!("K={}", h.max_depth_k));
                                let c = match h.predicted_state {
                                    PredictedState::Activated => pal.activated_strong,
                                    PredictedState::Inhibited => pal.inhibited_strong,
                                    PredictedState::Neutral => pal.text_muted,
                                };
                                ui.colored_label(c, h.predicted_state.label());
                                ui.colored_label(c, format!("{:+.2}", h.activation_z));
                                ui.label(format!("{:.2e}", h.overlap_p_value));
                                ui.label(format!("{:.3e}", h.network_bias_p_value));
                                ui.end_row();
                            }
                        });
                });
            });

            if let Some(hit) = rep.causal_networks.get(self.selected_cna_idx) {
                cols[1].heading(format!(
                    "Root Hypothesis: {} ({}, K={})",
                    hit.root_name,
                    hit.predicted_state.label(),
                    hit.max_depth_k
                ));
                if !hit.intermediate_regulators.is_empty() {
                    cols[1].label(format!(
                        "Intermediate Signaling / TF Relays: {}",
                        hit.intermediate_regulators.join(" -> ")
                    ));
                }
                cols[1].separator();
                cols[1].label(RichText::new("Multi-Hop Causal Paths to Dataset Genes").strong());
                cols[1].push_id("cna_paths_scroll", |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for p in &hit.sample_paths {
                            let chain = p.nodes.join("  ->  ");
                            let c = if p.consistent {
                                pal.activated_strong
                            } else {
                                pal.inconsistent_amber
                            };
                            ui.horizontal_wrapped(|ui| {
                                ui.colored_label(c, if p.net_sign > 0 { "[+]" } else { "[-]" });
                                ui.label(&chain);
                                ui.label(
                                    RichText::new(format!("(log2FC = {:+.2})", p.observed_log2fc))
                                        .color(pal.text_muted),
                                );
                            });
                        }
                    });
                });
            }
        });
    }

    fn render_dea_treemap_tab(
        &self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = &self.reports[self.active_cohort_idx];
        ui.horizontal(|ui| {
            ui.heading("Downstream Effects Analysis — Biological Process & Disease TreeMap");
            if ui.button("Export TreeMap SVG").clicked() {
                actions.push(GuiAction::ExportTreeMapSvgDialog);
            }
        });

        ui.columns(2, |cols| {
            let sz = Vec2::new(
                cols[0].available_width().max(300.0),
                (cols[0].available_height() - 10.0).max(340.0),
            );
            let (resp, painter) = cols[0].allocate_painter(sz, Sense::click());
            let rect = resp.rect;
            painter.rect_filled(rect, 4.0, pal.bg_canvas);

            let click_pt = if resp.clicked() {
                resp.interact_pointer_pos()
            } else {
                None
            };
            let hover_pt = resp.hover_pos();

            for r in &rep.treemap_rects {
                let cell = Rect::from_min_size(
                    Pos2::new(
                        rect.left() + r.x * rect.width() + 1.5,
                        rect.top() + r.y * rect.height() + 1.5,
                    ),
                    Vec2::new(
                        (r.w * rect.width() - 3.0).max(2.0),
                        (r.h * rect.height() - 3.0).max(2.0),
                    ),
                );
                let fill = pal.color_for_z(r.activation_z);
                let is_sel = self.selected_dea_term.as_deref() == Some(r.term_id.as_str());
                painter.rect_filled(cell, 3.0, fill);
                painter.rect_stroke(
                    cell,
                    3.0,
                    Stroke::new(if is_sel { 2.4 } else { 1.0 }, pal.text_primary),
                );

                if cell.width() > 68.0 && cell.height() > 26.0 {
                    painter.text(
                        Pos2::new(cell.left() + 6.0, cell.top() + 12.0),
                        Align2::LEFT_CENTER,
                        &r.term_name,
                        FontId::proportional(11.0),
                        pal.text_primary,
                    );
                    if cell.height() > 42.0 {
                        painter.text(
                            Pos2::new(cell.left() + 6.0, cell.top() + 27.0),
                            Align2::LEFT_CENTER,
                            format!("Z = {:+.2} ({} genes)", r.activation_z, r.overlap_count),
                            FontId::proportional(10.0),
                            pal.text_primary,
                        );
                    }
                }

                if let Some(cp) = click_pt {
                    if cell.contains(cp) {
                        actions.push(GuiAction::SelectDeaTerm(r.term_id.clone()));
                    }
                }
                if let Some(hp) = hover_pt {
                    if cell.contains(hp) {
                        resp.clone().on_hover_text(format!(
                            "{}\nCategory: {}\nActivation Z: {:+.2}\n-log10(P): {:.2}",
                            r.term_name, r.category, r.activation_z, r.neg_log10_p
                        ));
                    }
                }
            }

            cols[1].heading("Enriched Processes & Disease Phenotypes");
            cols[1].push_id("dea_list_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("dea_grid")
                        .striped(true)
                        .num_columns(5)
                        .spacing([10.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Function / Disease");
                            ui.strong("Prediction");
                            ui.strong("Act Z");
                            ui.strong("BH-FDR");
                            ui.strong("DEGs");
                            ui.end_row();

                            for d in &rep.downstream_effects {
                                let is_sel =
                                    self.selected_dea_term.as_deref() == Some(d.term_id.as_str());
                                if ui.selectable_label(is_sel, &d.term_name).clicked() {
                                    actions.push(GuiAction::SelectDeaTerm(d.term_id.clone()));
                                }
                                let c = pal.color_for_z(d.activation_z);
                                ui.colored_label(c, &d.predicted_direction);
                                ui.colored_label(c, format!("{:+.2}", d.activation_z));
                                ui.label(format!("{:.2e}", d.fdr));
                                ui.label(d.overlap_genes.join(", "));
                                ui.end_row();
                            }
                        });
                });
            });
        });
    }

    fn render_regulator_effects_tab(
        &self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = &self.reports[self.active_cohort_idx];
        ui.heading("Regulator Effects — 3-Tier Causal Mechanism Explorer");
        ui.horizontal(|ui| {
            ui.label("Select Causal Flow Hypothesis:");
            for (idx, re) in rep.regulator_effects.iter().take(8).enumerate() {
                let title = re
                    .regulators
                    .first()
                    .map(|(r, _)| format!("{r} (CS={:.1})", re.consistency_score))
                    .unwrap_or_else(|| format!("Network #{}", idx + 1));
                if ui
                    .selectable_label(self.selected_reg_effect_idx == idx, title)
                    .clicked()
                {
                    actions.push(GuiAction::SelectRegEffectIndex(idx));
                }
            }
        });

        if let Some(hit) = rep.regulator_effects.get(self.selected_reg_effect_idx) {
            let sz = Vec2::new(
                ui.available_width().max(320.0),
                (ui.available_height() - 12.0).max(340.0),
            );
            let (resp, painter) = ui.allocate_painter(sz, Sense::hover());
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, pal.bg_canvas);
            painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

            let y_top = rect.top() + 58.0;
            let y_mid = rect.center().y;
            let y_bot = rect.bottom() - 62.0;

            let mut reg_pts = Vec::new();
            for (i, (rid, z)) in hit.regulators.iter().enumerate() {
                let frac = (i as f32 + 0.5) / (hit.regulators.len().max(1) as f32);
                let p = Pos2::new(rect.left() + 120.0 + frac * (rect.width() - 240.0), y_top);
                reg_pts.push((rid.as_str(), *z, p));
            }

            let mut gene_pts = Vec::new();
            for (i, (gid, fc)) in hit.shared_dataset_targets.iter().take(12).enumerate() {
                let n = hit.shared_dataset_targets.len().min(12).max(1) as f32;
                let frac = (i as f32 + 0.5) / n;
                let p = Pos2::new(rect.left() + 45.0 + frac * (rect.width() - 90.0), y_mid);
                gene_pts.push((gid.as_str(), *fc, p));
            }

            let mut pheno_pts = Vec::new();
            for (i, (_pid, pname, pz)) in hit.downstream_phenotypes.iter().enumerate() {
                let n = hit.downstream_phenotypes.len().max(1) as f32;
                let frac = (i as f32 + 0.5) / n;
                let p = Pos2::new(rect.left() + 90.0 + frac * (rect.width() - 180.0), y_bot);
                pheno_pts.push((pname.as_str(), *pz, p));
            }

            for &(_, z, rp) in &reg_pts {
                let s = if z >= 0.0 {
                    EdgeSign::Activating
                } else {
                    EdgeSign::Inhibiting
                };
                for &(_, _, gp) in &gene_pts {
                    draw_directed_causal_edge(
                        &painter,
                        rp,
                        gp,
                        20.0,
                        s,
                        Stroke::new(1.5, pal.color_for_z(z)),
                    );
                }
            }

            for &(_, _, gp) in &gene_pts {
                for &(_, pz, pp) in &pheno_pts {
                    let s = if pz >= 0.0 {
                        EdgeSign::Activating
                    } else {
                        EdgeSign::Inhibiting
                    };
                    draw_directed_causal_edge(
                        &painter,
                        gp,
                        pp,
                        24.0,
                        s,
                        Stroke::new(1.5, pal.color_for_z(pz)),
                    );
                }
            }

            for (rid, z, p) in reg_pts {
                draw_entity_glyph(
                    &painter,
                    p,
                    Vec2::new(42.0, 18.0),
                    EntityKind::TranscriptionRegulator,
                    pal.color_for_z(z),
                    Stroke::new(2.0, pal.text_primary),
                );
                painter.text(
                    p,
                    Align2::CENTER_CENTER,
                    rid,
                    FontId::proportional(11.5),
                    pal.text_primary,
                );
            }

            for (gid, fc, p) in gene_pts {
                draw_entity_glyph(
                    &painter,
                    p,
                    Vec2::new(28.0, 14.0),
                    EntityKind::GeneOther,
                    pal.color_for_log2fc(fc),
                    Stroke::new(1.5, pal.text_primary),
                );
                painter.text(
                    p,
                    Align2::CENTER_CENTER,
                    gid,
                    FontId::proportional(10.5),
                    pal.text_primary,
                );
            }

            for (pname, pz, p) in pheno_pts {
                draw_entity_glyph(
                    &painter,
                    p,
                    Vec2::new(92.0, 19.0),
                    EntityKind::BiologicalProcess,
                    pal.color_for_z(pz),
                    Stroke::new(2.0, pal.text_primary),
                );
                painter.text(
                    p,
                    Align2::CENTER_CENTER,
                    pname,
                    FontId::proportional(10.5),
                    pal.text_primary,
                );
            }
        } else {
            ui.label("No multi-tier Regulator Effects pass current significance thresholds.");
        }
    }

    fn render_comparison_perturb_tab(
        &self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        ui.columns(2, |cols| {
            cols[0].heading("Multi-Cohort Comparison Heatmap");
            cols[0].horizontal_wrapped(|ui| {
                for m in [
                    ComparisonMetric::UpstreamRegulatorsZ,
                    ComparisonMetric::CanonicalPathwaysZ,
                    ComparisonMetric::DownstreamEffectsZ,
                ] {
                    if ui
                        .selectable_label(self.comparison_metric == m, m.label())
                        .clicked()
                    {
                        actions.push(GuiAction::SetComparisonMetric(m));
                    }
                }
            });

            let comp = run_comparison_analysis(&self.reports, self.comparison_metric, 18);
            cols[0].push_id("comp_heatmap_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("comp_matrix_grid")
                        .spacing([6.0, 4.0])
                        .show(ui, |ui| {
                            ui.strong("Entity");
                            for cname in &comp.cohort_names {
                                let short: String = cname.chars().take(18).collect();
                                ui.strong(short);
                            }
                            ui.end_row();

                            for (r_idx, label) in comp.entity_labels.iter().enumerate() {
                                ui.label(label);
                                for &z in &comp.values[r_idx] {
                                    let (cell_r, cell_p) =
                                        ui.allocate_painter(Vec2::new(66.0, 20.0), Sense::hover());
                                    cell_p.rect_filled(cell_r.rect, 3.0, pal.color_for_z(z));
                                    cell_p.text(
                                        cell_r.rect.center(),
                                        Align2::CENTER_CENTER,
                                        format!("{z:+.2}"),
                                        FontId::proportional(10.5),
                                        pal.text_primary,
                                    );
                                }
                                ui.end_row();
                            }
                        });
                });
            });

            let rep = &self.reports[self.active_cohort_idx];
            cols[1].heading("Perturbation Connectivity (SigCom LINCS Open Signatures)");
            cols[1].label(
                RichText::new(
                    "Negative scores (blue) identify candidate small-molecule Reversers that oppose the disease/stimulus expression state; positive scores (orange) identify concordant Mimickers.",
                )
                .size(11.5)
                .color(pal.text_muted),
            );

            cols[1].push_id("perturb_scroll", |ui| {
                egui::ScrollArea::both().show(ui, |ui| {
                    egui::Grid::new("perturb_grid")
                        .striped(true)
                        .num_columns(5)
                        .spacing([8.0, 6.0])
                        .max_col_width(150.0)
                        .show(ui, |ui| {
                            ui.strong("Perturbagen");
                            ui.strong("Class");
                            ui.strong("Connectivity");
                            ui.strong("Role");
                            ui.strong("Mechanism of Action");
                            ui.end_row();

                            for p in &rep.perturbation_matches {
                                ui.add(egui::Label::new(&p.perturbagen).wrap());
                                ui.add(egui::Label::new(p.kind.label()).wrap());
                                let c = pal.color_for_z(p.connectivity_score / 25.0);
                                ui.colored_label(c, format!("{:+.1}", p.connectivity_score));
                                ui.add(
                                    egui::Label::new(RichText::new(&p.mode).color(c)).wrap(),
                                );
                                ui.add(egui::Label::new(&p.mechanism_of_action).wrap());
                                ui.end_row();
                            }
                        });
                });
            });
        });
    }

    fn render_open_streams_tab(
        &self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.heading("SQLite 3 WAL Cache & Open-Source Batch Downloader");
            ui.add_space(10.0);
            if ui
                .add_enabled(
                    !self.worker_busy,
                    egui::Button::new("Batch Download All Pathways (GPML + Reactome)"),
                )
                .clicked()
            {
                actions.push(GuiAction::StartBatchDownload {
                    pathways: true,
                    regulatory: false,
                    force: false,
                });
            }
            if ui
                .add_enabled(
                    !self.worker_busy,
                    egui::Button::new("Check & Sync All Updates (ETag / TTL)"),
                )
                .clicked()
            {
                actions.push(GuiAction::StartBatchDownload {
                    pathways: true,
                    regulatory: true,
                    force: false,
                });
            }
            if ui
                .add_enabled(!self.worker_busy, egui::Button::new("Force Full Re-Sync"))
                .clicked()
            {
                actions.push(GuiAction::StartBatchDownload {
                    pathways: true,
                    regulatory: true,
                    force: true,
                });
            }
            if self.worker_busy && ui.button("Cancel Download").clicked() {
                actions.push(GuiAction::CancelBatchDownload);
            }
        });

        ui.add_space(4.0);
        ui.add(
            egui::ProgressBar::new(self.batch_progress_fraction)
                .text(&self.batch_progress_label)
                .animate(self.worker_busy),
        );
        ui.add_space(6.0);

        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.label(RichText::new("SQLite WAL Database Telemetry").strong());
                egui::Grid::new("sqlite_stats_grid")
                    .num_columns(2)
                    .spacing([14.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Database File:");
                        ui.monospace(&self.db_stats.db_path);
                        ui.end_row();

                        ui.label("Journal Mode & Footprint:");
                        ui.label(format!(
                            "{} | {:.1} KB on disk | KG rev {}",
                            if self.db_stats.wal_mode { "WAL" } else { "MEMORY/DELETE" },
                            self.db_stats.file_size_bytes as f64 / 1024.0,
                            self.db_stats.kg_revision_hash
                        ));
                        ui.end_row();

                        ui.label("Cached Canonical Pathways:");
                        ui.label(format!(
                            "{} total ({} with 2D diagram coordinates)",
                            self.db_stats.cached_pathways, self.db_stats.cached_diagrams_with_coords
                        ));
                        ui.end_row();

                        ui.label("Cached Entities & Edges:");
                        ui.label(format!(
                            "{} molecules/phenotypes | {} signed causal edges",
                            self.db_stats.cached_entities, self.db_stats.cached_edges
                        ));
                        ui.end_row();

                        ui.label("Analysis Result Cache:");
                        ui.label(format!(
                            "{} stored reports ({} hits / {} misses)",
                            self.db_stats.cached_analysis_reports,
                            self.analysis_cache_hits,
                            self.analysis_cache_misses
                        ));
                        ui.end_row();
                    });

                ui.horizontal_wrapped(|ui| {
                    if ui
                        .add_enabled(!self.worker_busy, egui::Button::new("Import Local .zip / .gpml / .gmt..."))
                        .clicked()
                    {
                        actions.push(GuiAction::ImportStreamFileDialog);
                    }
                    if ui.button("Clear Analysis Cache").clicked() {
                        actions.push(GuiAction::ClearAnalysisCache);
                    }
                    if ui.button("SQLite VACUUM").clicked() {
                        actions.push(GuiAction::VacuumSqlite);
                    }
                    if ui.button("Export .json.gz").clicked() {
                        actions.push(GuiAction::SaveSnapshotDialog);
                    }
                });
            });

            cols[1].group(|ui| {
                ui.label(RichText::new("Live Stream Sync & SQLite WAL Log").strong());
                ui.push_id("sqlite_sync_log_scroll", |ui| {
                    egui::ScrollArea::vertical()
                        .max_height(135.0)
                        .stick_to_bottom(true)
                        .show(ui, |ui| {
                            for line in &self.sync_logs {
                                ui.label(
                                    RichText::new(line)
                                        .monospace()
                                        .size(11.0)
                                        .color(pal.text_muted),
                                );
                            }
                        });
                });
            });
        });

        ui.separator();
        ui.push_id("streams_registry_scroll", |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("streams_grid")
                    .striped(true)
                    .num_columns(6)
                    .spacing([12.0, 7.0])
                    .show(ui, |ui| {
                        ui.strong("Open Data Stream");
                        ui.strong("Analytical Role");
                        ui.strong("Format");
                        ui.strong("License");
                        ui.strong("SQLite Cache Status");
                        ui.strong("Per-Stream Sync");
                        ui.end_row();

                        for (idx, st) in self.streams.iter().enumerate() {
                            ui.label(RichText::new(&st.name).strong());
                            ui.label(&st.role_in_analysis);
                            ui.label(st.format.label());
                            ui.label(&st.license);
                            let sha_note = st
                                .last_sha256
                                .as_deref()
                                .map(|h| format!(" [sha256:{h}]"))
                                .unwrap_or_default();
                            ui.label(format!("{}{sha_note}", st.last_sync_status));
                            ui.horizontal(|ui| {
                                if ui
                                    .add_enabled(!self.worker_busy, egui::Button::new("Sync"))
                                    .clicked()
                                {
                                    actions.push(GuiAction::SyncRemoteStream { idx, force: false });
                                }
                                if ui
                                    .add_enabled(!self.worker_busy, egui::Button::new("Force"))
                                    .clicked()
                                {
                                    actions.push(GuiAction::SyncRemoteStream { idx, force: true });
                                }
                            });
                            ui.end_row();
                        }
                    });
            });
        });
    }
}

struct TargetGeneRef<'a> {
    sym: &'a str,
    fc: f64,
    sign: EdgeSign,
}
