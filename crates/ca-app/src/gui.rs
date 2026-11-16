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
    apply_workbench_visuals, draw_bezier_causal_edge, draw_dna_double_helix, draw_entity_glyph,
    draw_gold_pulse_halo, draw_phospho_badge, ScientificPalette, WorkbenchPaletteMode,
};
use ca_core::{
    build_mechanistic_network, builtin_cohorts, compute_pathway_world_layout_cached,
    compute_tiered_causal_layout, parse_expression_table, recompute_and_cache_pathway_layout,
    run_comparison_analysis, run_core_analysis, seed_kb::build_seed_knowledge_graph,
    simulate_map_in_kg, ComparisonMetric, CoreAnalysisConfig, CoreAnalysisReport, EdgeSign,
    EntityKind, ExpressionDataset, KnowledgeGraph, PredictedState, SubcellularLocation,
    TieredNodeSpec, MIN_NODE_DX, MIN_NODE_DY,
};
use ca_streams::{
    client::{fetch_remote_stream_conditional, import_local_stream_file},
    default_stream_registry, fetch_and_cache_reactome_diagram, load_compressed_snapshot,
    run_batch_download_and_update, save_compressed_snapshot, BatchProgressEvent, BatchSyncSummary,
    SqliteCacheStats, SqliteKnowledgeStore, StreamDescriptor,
};
use egui::{Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Ui, Vec2};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::PathBuf;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc::{self, Receiver, Sender},
    Arc,
};
use std::thread;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkbenchTab {
    CohortVolcano,
    PathwayAtlas,
    CausalRegulators,
    PhenotypeDrugMatch,
    KnowledgeDatabases,
}

impl WorkbenchTab {
    pub const ALL: [Self; 5] = [
        Self::CohortVolcano,
        Self::PathwayAtlas,
        Self::CausalRegulators,
        Self::PhenotypeDrugMatch,
        Self::KnowledgeDatabases,
    ];

    pub fn title(self) -> &'static str {
        match self {
            Self::CohortVolcano => "1. Cohort & Volcano",
            Self::PathwayAtlas => "2. Pathway Atlas",
            Self::CausalRegulators => "3. Causal Regulators",
            Self::PhenotypeDrugMatch => "4. Phenotype & Drug Match",
            Self::KnowledgeDatabases => "5. Knowledge Databases",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CausalRegulatorsSubTab {
    UpstreamRegulators,
    MultiHopNetworks,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhenotypeDrugSubTab {
    ComparisonAndReversers,
    DownstreamTreeMap,
    RegulatorEffectsFlow,
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
    ToggleSidebarHidden,
    FitPathwayViewport,
    CenterSelectedPathwayNode,
    RelayoutActivePathway,
    TogglePathwayEgoFilter,
    UnpinPathwayNodes(String),
    JumpToGene {
        symbol: String,
        pathway_id: Option<String>,
    },
}

pub struct CausalAtlasApp {
    palette_mode: WorkbenchPaletteMode,
    active_tab: WorkbenchTab,
    causal_subtab: CausalRegulatorsSubTab,
    phenotype_subtab: PhenotypeDrugSubTab,
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

    // Global Gene/Protein Omnibox (Ctrl+F) & Context Strip
    omnibox_query: String,
    omnibox_dropdown_open: bool,
    omnibox_focus_requested: bool,
    omnibox_selected_gene: Option<String>,
    highlight_pulse_gene: Option<String>,
    highlight_pulse_timer: f32,

    gene_search: String,
    pathway_search: String,
    regulator_search: String,
    selected_pathway_idx: usize,
    override_pathway_id: Option<String>,
    selected_pathway_node: Option<String>,
    map_user_clamps: BTreeMap<String, f64>,
    selected_ura_idx: usize,
    mechanistic_node_offsets: BTreeMap<String, Vec2>,
    dragging_ura_node: Option<String>,
    ura_zoom: f32,
    ura_pan: Vec2,
    ura_fit_pending: bool,
    selected_cna_idx: usize,
    selected_dea_term: Option<String>,
    selected_reg_effect_idx: usize,
    re_node_offsets: BTreeMap<String, Vec2>,
    dragging_re_node: Option<String>,
    re_zoom: f32,
    re_pan: Vec2,
    re_fit_pending: bool,
    comparison_metric: ComparisonMetric,
    pathway_fullscreen: bool,
    sidebar_hidden: bool,
    pathway_zoom: f32,
    pathway_pan: Vec2,
    pathway_fit_pending: bool,
    pathway_ego_filter: bool,
    dragging_pathway_node: Option<String>,
    pinned_pathway_nodes: BTreeMap<String, BTreeMap<String, Pos2>>,
    camera_pan_pending_node: Option<String>,
    camera_anim_start_pan: Vec2,
    camera_anim_target_pan: Vec2,
    camera_anim_start_zoom: f32,
    camera_anim_target_zoom: f32,
    camera_anim_progress: f32,

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
            active_tab: WorkbenchTab::PathwayAtlas,
            causal_subtab: CausalRegulatorsSubTab::UpstreamRegulators,
            phenotype_subtab: PhenotypeDrugSubTab::ComparisonAndReversers,
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
            omnibox_query: String::new(),
            omnibox_dropdown_open: false,
            omnibox_focus_requested: false,
            omnibox_selected_gene: None,
            highlight_pulse_gene: None,
            highlight_pulse_timer: 0.0,
            gene_search: String::new(),
            pathway_search: String::new(),
            regulator_search: String::new(),
            selected_pathway_idx: 0,
            override_pathway_id: None,
            selected_pathway_node: None,
            map_user_clamps: BTreeMap::new(),
            selected_ura_idx: 0,
            mechanistic_node_offsets: BTreeMap::new(),
            dragging_ura_node: None,
            ura_zoom: 1.0,
            ura_pan: Vec2::ZERO,
            ura_fit_pending: true,
            selected_cna_idx: 0,
            selected_dea_term: None,
            selected_reg_effect_idx: 0,
            re_node_offsets: BTreeMap::new(),
            dragging_re_node: None,
            re_zoom: 1.0,
            re_pan: Vec2::ZERO,
            re_fit_pending: true,
            comparison_metric: ComparisonMetric::UpstreamRegulatorsZ,
            pathway_fullscreen: false,
            sidebar_hidden: false,
            pathway_zoom: 1.0,
            pathway_pan: Vec2::ZERO,
            pathway_fit_pending: true,
            pathway_ego_filter: false,
            dragging_pathway_node: None,
            pinned_pathway_nodes: BTreeMap::new(),
            camera_pan_pending_node: None,
            camera_anim_start_pan: Vec2::ZERO,
            camera_anim_target_pan: Vec2::ZERO,
            camera_anim_start_zoom: 1.0,
            camera_anim_target_zoom: 1.0,
            camera_anim_progress: 1.0,
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

    fn active_pathway_id(&self) -> Option<String> {
        if let Some(ref oid) = self.override_pathway_id {
            if self.kg.pathways.iter().any(|p| &p.id == oid) {
                return Some(oid.clone());
            }
        }
        let rep = &self.reports[self.active_cohort_idx];
        rep.canonical_pathways
            .get(self.selected_pathway_idx)
            .map(|h| h.pathway_id.clone())
            .or_else(|| self.kg.pathways.first().map(|p| p.id.clone()))
    }

    fn apply_action(&mut self, ctx: &egui::Context, action: GuiAction) {
        match action {
            GuiAction::OpenStreamsTab => {
                self.active_tab = WorkbenchTab::KnowledgeDatabases;
            }
            GuiAction::SelectCohort(idx) => {
                if idx < self.cohorts.len() {
                    self.active_cohort_idx = idx;
                    self.selected_pathway_idx = 0;
                    self.override_pathway_id = None;
                    self.selected_ura_idx = 0;
                    self.selected_cna_idx = 0;
                    self.selected_reg_effect_idx = 0;
                    self.map_user_clamps.clear();
                    self.mechanistic_node_offsets.clear();
                    self.re_node_offsets.clear();
                    self.pathway_fit_pending = true;
                    self.ura_fit_pending = true;
                    self.re_fit_pending = true;
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
                if let Some(pw_id) = self.active_pathway_id() {
                    if let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id) {
                        let active = self.cohorts[self.active_cohort_idx].filter_active(
                            self.config.abs_log2fc_cutoff,
                            self.config.fdr_cutoff,
                            &self.kg,
                        );
                        let map_res = simulate_map_in_kg(
                            pw_def,
                            &active,
                            &self.kg,
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
                self.override_pathway_id = None;
                self.selected_pathway_node = None;
                self.pathway_fit_pending = true;
                self.camera_anim_progress = 1.0;
                self.camera_pan_pending_node = None;
            }
            GuiAction::SelectPathwayNode(nid) => {
                self.selected_pathway_node = Some(nid.clone());
                self.omnibox_selected_gene = Some(nid);
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
                self.ura_fit_pending = true;
            }
            GuiAction::SelectCnaIndex(idx) => {
                self.selected_cna_idx = idx;
            }
            GuiAction::SelectDeaTerm(tid) => {
                self.selected_dea_term = Some(tid);
            }
            GuiAction::SelectRegEffectIndex(idx) => {
                self.selected_reg_effect_idx = idx;
                self.re_node_offsets.clear();
                self.re_fit_pending = true;
            }
            GuiAction::SetComparisonMetric(m) => {
                self.comparison_metric = m;
            }
            GuiAction::TogglePathwayFullscreen => {
                self.pathway_fullscreen = !self.pathway_fullscreen;
                self.pathway_fit_pending = true;
            }
            GuiAction::ToggleSidebarHidden => {
                self.sidebar_hidden = !self.sidebar_hidden;
            }
            GuiAction::FitPathwayViewport => {
                self.pathway_fit_pending = true;
                self.camera_anim_progress = 1.0;
                self.camera_pan_pending_node = None;
            }
            GuiAction::CenterSelectedPathwayNode => {
                let target = self
                    .selected_pathway_node
                    .clone()
                    .or_else(|| self.omnibox_selected_gene.clone());
                if let Some(nid) = target {
                    self.camera_pan_pending_node = Some(nid.clone());
                    self.highlight_pulse_gene = Some(nid);
                    self.highlight_pulse_timer = 2.8;
                } else {
                    self.pathway_fit_pending = true;
                }
            }
            GuiAction::RelayoutActivePathway => {
                if let Some(pw_id) = self.active_pathway_id() {
                    self.pinned_pathway_nodes.remove(&pw_id);
                    if let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id) {
                        recompute_and_cache_pathway_layout(pw_def, &self.kg);
                    }
                    self.pathway_fit_pending = true;
                    self.status_banner = format!(
                        "Recomputed 120-iteration radiative compartment layout for {pw_id}"
                    );
                }
            }
            GuiAction::TogglePathwayEgoFilter => {
                self.pathway_ego_filter = !self.pathway_ego_filter;
                if self.pathway_ego_filter && self.selected_pathway_node.is_none() {
                    if let Some(pw_id) = self.active_pathway_id() {
                        if let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id) {
                            let wl = compute_pathway_world_layout_cached(pw_def, &self.kg);
                            let mut deg: BTreeMap<&str, usize> = BTreeMap::new();
                            for e in &wl.edges {
                                *deg.entry(e.from.as_str()).or_insert(0) += 1;
                                *deg.entry(e.to.as_str()).or_insert(0) += 1;
                            }
                            if let Some((best, _)) = deg
                                .into_iter()
                                .max_by_key(|&(_, d)| d)
                                .or_else(|| wl.positions.keys().next().map(|k| (k.as_str(), 0)))
                            {
                                self.selected_pathway_node = Some(best.to_string());
                                self.omnibox_selected_gene = Some(best.to_string());
                            }
                        }
                    }
                }
            }
            GuiAction::UnpinPathwayNodes(pw_id) => {
                self.pinned_pathway_nodes.remove(&pw_id);
            }
            GuiAction::JumpToGene { symbol, pathway_id } => {
                let sym = symbol.trim().to_uppercase();
                self.omnibox_selected_gene = Some(sym.clone());
                self.omnibox_dropdown_open = false;
                self.highlight_pulse_gene = Some(sym.clone());
                self.highlight_pulse_timer = 3.2;

                let chosen_pw = pathway_id.or_else(|| {
                    let rep = &self.reports[self.active_cohort_idx];
                    for hit in &rep.canonical_pathways {
                        if let Some(pw) = self.kg.pathways.iter().find(|p| p.id == hit.pathway_id) {
                            if pw.members.iter().any(|m| m == &sym)
                                || pw.topology_edges.iter().any(|e| e.from == sym || e.to == sym)
                            {
                                return Some(pw.id.clone());
                            }
                        }
                    }
                    self.kg
                        .pathways
                        .iter()
                        .find(|pw| {
                            pw.members.iter().any(|m| m == &sym)
                                || pw.topology_edges.iter().any(|e| e.from == sym || e.to == sym)
                        })
                        .map(|p| p.id.clone())
                });

                if let Some(pw_id) = chosen_pw {
                    let rep_idx = self.reports[self.active_cohort_idx]
                        .canonical_pathways
                        .iter()
                        .position(|p| p.pathway_id == pw_id);
                    if let Some(idx) = rep_idx {
                        self.selected_pathway_idx = idx;
                        self.override_pathway_id = None;
                    } else {
                        self.override_pathway_id = Some(pw_id.clone());
                    }
                    self.active_tab = WorkbenchTab::PathwayAtlas;
                    self.selected_pathway_node = Some(sym.clone());
                    self.camera_pan_pending_node = Some(sym.clone());
                    self.status_banner = format!("Focused {sym} inside Pathway {pw_id}");
                } else {
                    self.active_tab = WorkbenchTab::CohortVolcano;
                    self.gene_search = sym.clone();
                    self.status_banner =
                        format!("Highlighted {sym} in Cohort Differential Expression table");
                }
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
        self.sidebar_hidden = false;
        self.pathway_fullscreen = false;
        self.pathway_ego_filter = false;
        self.pathway_fit_pending = true;
        self.ura_fit_pending = true;
        self.re_fit_pending = true;
        match view {
            2 => {
                self.active_tab = WorkbenchTab::PathwayAtlas;
                self.selected_pathway_idx = 0;
                self.override_pathway_id = None;
                self.selected_pathway_node = Some("VCAM1".to_string());
                self.omnibox_selected_gene = Some("VCAM1".to_string());
                self.highlight_pulse_gene = Some("VCAM1".to_string());
                self.highlight_pulse_timer = 2.0;
            }
            3 => {
                self.active_tab = WorkbenchTab::CausalRegulators;
                self.causal_subtab = CausalRegulatorsSubTab::UpstreamRegulators;
                self.omnibox_selected_gene = Some("TNF".to_string());
            }
            _ => {
                self.active_tab = WorkbenchTab::PhenotypeDrugMatch;
                self.phenotype_subtab = PhenotypeDrugSubTab::ComparisonAndReversers;
            }
        }
    }

    pub fn render_ui_frame(&mut self, ctx: &egui::Context) {
        self.poll_worker();
        let dt = ctx.input(|i| i.unstable_dt).clamp(0.001, 0.05);
        if self.highlight_pulse_timer > 0.0 {
            self.highlight_pulse_timer = (self.highlight_pulse_timer - dt).max(0.0);
            if self.highlight_pulse_timer > 0.0 {
                ctx.request_repaint();
            }
        }
        if self.camera_anim_progress < 1.0 {
            ctx.request_repaint();
        }
        if self.worker_busy {
            ctx.request_repaint();
        }

        // Global Ctrl+F / Cmd+F shortcut to focus Omnibox
        if ctx.input(|i| (i.modifiers.ctrl || i.modifiers.command) && i.key_pressed(egui::Key::F)) {
            self.omnibox_focus_requested = true;
            self.omnibox_dropdown_open = true;
        }
        if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            self.omnibox_dropdown_open = false;
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
            ui.separator();
            self.render_global_omnibox_and_context_strip(ctx, ui, &pal, &mut actions);
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

        if !self.sidebar_hidden && !self.pathway_fullscreen {
            egui::SidePanel::left("left_controls")
                .resizable(true)
                .default_width(250.0)
                .min_width(220.0)
                .max_width(320.0)
                .show(ctx, |ui| {
                    ui.heading("Cohort & Cutoffs");
                    ui.add_space(4.0);

                    egui::ComboBox::from_label("Cohort")
                        .selected_text(&self.cohorts[self.active_cohort_idx].name)
                        .width(190.0)
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
                                RichText::new(format!("{}", sum.upregulated_count))
                                    .color(pal.expr_up),
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
        }

        egui::CentralPanel::default().show(ctx, |ui| match self.active_tab {
            WorkbenchTab::CohortVolcano => {
                self.render_dataset_volcano_tab(ui, &pal);
            }
            WorkbenchTab::PathwayAtlas => {
                self.render_canonical_map_tab(ui, &pal, &mut actions);
            }
            WorkbenchTab::CausalRegulators => {
                ui.horizontal(|ui| {
                    if ui
                        .selectable_label(
                            self.causal_subtab == CausalRegulatorsSubTab::UpstreamRegulators,
                            "Upstream Regulators & Mechanistic Network",
                        )
                        .clicked()
                    {
                        self.causal_subtab = CausalRegulatorsSubTab::UpstreamRegulators;
                    }
                    if ui
                        .selectable_label(
                            self.causal_subtab == CausalRegulatorsSubTab::MultiHopNetworks,
                            "Multi-Hop Causal Networks (CNA)",
                        )
                        .clicked()
                    {
                        self.causal_subtab = CausalRegulatorsSubTab::MultiHopNetworks;
                    }
                });
                ui.separator();
                match self.causal_subtab {
                    CausalRegulatorsSubTab::UpstreamRegulators => {
                        self.render_ura_mechanistic_tab(ui, &pal, &mut actions);
                    }
                    CausalRegulatorsSubTab::MultiHopNetworks => {
                        self.render_cna_tab(ui, &pal, &mut actions);
                    }
                }
            }
            WorkbenchTab::PhenotypeDrugMatch => {
                ui.horizontal(|ui| {
                    if ui
                        .selectable_label(
                            self.phenotype_subtab == PhenotypeDrugSubTab::ComparisonAndReversers,
                            "Multi-Cohort Comparison & Drug Reversers",
                        )
                        .clicked()
                    {
                        self.phenotype_subtab = PhenotypeDrugSubTab::ComparisonAndReversers;
                    }
                    if ui
                        .selectable_label(
                            self.phenotype_subtab == PhenotypeDrugSubTab::DownstreamTreeMap,
                            "Downstream Effects TreeMap",
                        )
                        .clicked()
                    {
                        self.phenotype_subtab = PhenotypeDrugSubTab::DownstreamTreeMap;
                    }
                    if ui
                        .selectable_label(
                            self.phenotype_subtab == PhenotypeDrugSubTab::RegulatorEffectsFlow,
                            "Regulator Effects (3-Tier Flow)",
                        )
                        .clicked()
                    {
                        self.phenotype_subtab = PhenotypeDrugSubTab::RegulatorEffectsFlow;
                    }
                });
                ui.separator();
                match self.phenotype_subtab {
                    PhenotypeDrugSubTab::ComparisonAndReversers => {
                        self.render_comparison_perturb_tab(ui, &pal, &mut actions);
                    }
                    PhenotypeDrugSubTab::DownstreamTreeMap => {
                        self.render_dea_treemap_tab(ui, &pal, &mut actions);
                    }
                    PhenotypeDrugSubTab::RegulatorEffectsFlow => {
                        self.render_regulator_effects_tab(ui, &pal, &mut actions);
                    }
                }
            }
            WorkbenchTab::KnowledgeDatabases => {
                self.render_open_streams_tab(ui, &pal, &mut actions);
            }
        });

        for act in actions {
            self.apply_action(ctx, act);
        }
    }

    fn render_global_omnibox_and_context_strip(
        &mut self,
        ctx: &egui::Context,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let ds = &self.cohorts[self.active_cohort_idx];
        let mut meas_by_sym: HashMap<&str, (f64, f64)> = HashMap::new();
        for r in &ds.rows {
            meas_by_sym.insert(r.symbol.as_str(), (r.log2_fc, r.fdr));
        }

        let mut sym_pathways: BTreeMap<&str, Vec<(&str, &str)>> = BTreeMap::new();
        for pw in &self.kg.pathways {
            for m in &pw.members {
                let list = sym_pathways.entry(m.as_str()).or_default();
                if !list.iter().any(|(id, _)| *id == pw.id.as_str()) {
                    list.push((pw.id.as_str(), pw.name.as_str()));
                }
            }
            for e in &pw.topology_edges {
                for ep in [e.from.as_str(), e.to.as_str()] {
                    let list = sym_pathways.entry(ep).or_default();
                    if !list.iter().any(|(id, _)| *id == pw.id.as_str()) {
                        list.push((pw.id.as_str(), pw.name.as_str()));
                    }
                }
            }
        }

        let mut edit_rect = Rect::NOTHING;
        let mut enter_pressed = false;

        ui.horizontal_wrapped(|ui| {
            ui.label(
                RichText::new("Gene / Protein Omnibox (Ctrl+F):")
                    .strong()
                    .size(12.0),
            );
            let resp = ui.add(
                egui::TextEdit::singleline(&mut self.omnibox_query)
                    .hint_text("Jump to symbol, receptor, TF, or drug...")
                    .desired_width(250.0),
            );
            edit_rect = resp.rect;

            if self.omnibox_focus_requested {
                resp.request_focus();
                self.omnibox_focus_requested = false;
            }
            if resp.gained_focus() || resp.changed() {
                self.omnibox_dropdown_open = !self.omnibox_query.trim().is_empty();
            }
            if resp.has_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                enter_pressed = true;
            }
            if !self.omnibox_query.is_empty() && ui.small_button("✕").clicked() {
                self.omnibox_query.clear();
                self.omnibox_dropdown_open = false;
            }

            // Context Strip for the active inspected gene/protein
            if let Some(ref active_sym) = self.omnibox_selected_gene {
                ui.separator();
                let (kind_lbl, comp_lbl) = self
                    .kg
                    .nodes
                    .get(active_sym)
                    .map(|n| (n.kind.label(), n.location.label()))
                    .unwrap_or(("Gene / Transcript", "Cytoplasm"));
                ui.label(
                    RichText::new(active_sym)
                        .strong()
                        .color(Color32::from_rgb(250, 204, 21)),
                );
                ui.label(
                    RichText::new(format!("{kind_lbl} | {comp_lbl}"))
                        .size(11.5)
                        .color(pal.text_muted),
                );
                if let Some(&(fc, fdr)) = meas_by_sym.get(active_sym.as_str()) {
                    let c = if fc > 0.0 { pal.expr_up } else { pal.expr_down };
                    ui.colored_label(
                        c,
                        RichText::new(format!("log2FC {fc:+.2} (FDR {fdr:.2e})")).size(11.5),
                    );
                }
                if let Some(pws) = sym_pathways.get(active_sym.as_str()) {
                    ui.label(RichText::new("Pathways:").size(11.5).color(pal.text_muted));
                    for &(pw_id, pw_name) in pws.iter().take(4) {
                        if ui
                            .small_button(pw_id)
                            .on_hover_text(pw_name)
                            .clicked()
                        {
                            actions.push(GuiAction::JumpToGene {
                                symbol: active_sym.clone(),
                                pathway_id: Some(pw_id.to_string()),
                            });
                        }
                    }
                }
            } else {
                ui.separator();
                ui.label(
                    RichText::new("Tip: Press Ctrl+F to locate any molecule across pathways and cohorts.")
                        .size(11.5)
                        .color(pal.text_muted),
                );
            }
        });

        let needle = self.omnibox_query.trim().to_ascii_uppercase();
        if !needle.is_empty() && (self.omnibox_dropdown_open || enter_pressed) {
            struct MatchItem<'a> {
                score: u8,
                symbol: String,
                kind: &'static str,
                compartment: &'static str,
                log2_fc: Option<f64>,
                pathways: Vec<(&'a str, &'a str)>,
            }

            let mut seen_syms: BTreeSet<String> = BTreeSet::new();
            let mut matches: Vec<MatchItem<'_>> = Vec::new();

            for (sym, node) in &self.kg.nodes {
                let id_up = sym.to_ascii_uppercase();
                let name_up = node.display_name.to_ascii_uppercase();
                let score = if id_up == needle {
                    4
                } else if id_up.starts_with(&needle) {
                    3
                } else if id_up.contains(&needle) {
                    2
                } else if name_up.contains(&needle) {
                    1
                } else {
                    0
                };
                if score > 0 {
                    seen_syms.insert(sym.clone());
                    matches.push(MatchItem {
                        score,
                        symbol: sym.clone(),
                        kind: node.kind.label(),
                        compartment: node.location.label(),
                        log2_fc: meas_by_sym.get(sym.as_str()).map(|v| v.0),
                        pathways: sym_pathways
                            .get(sym.as_str())
                            .cloned()
                            .unwrap_or_default(),
                    });
                }
            }
            for r in &ds.rows {
                if !seen_syms.contains(&r.symbol) {
                    let id_up = r.symbol.to_ascii_uppercase();
                    let score = if id_up == needle {
                        4
                    } else if id_up.starts_with(&needle) {
                        3
                    } else if id_up.contains(&needle) {
                        2
                    } else {
                        0
                    };
                    if score > 0 {
                        seen_syms.insert(r.symbol.clone());
                        matches.push(MatchItem {
                            score,
                            symbol: r.symbol.clone(),
                            kind: "Measured Transcript",
                            compartment: "Cytoplasm",
                            log2_fc: Some(r.log2_fc),
                            pathways: sym_pathways
                                .get(r.symbol.as_str())
                                .cloned()
                                .unwrap_or_default(),
                        });
                    }
                }
            }

            matches.sort_by(|a, b| {
                b.score
                    .cmp(&a.score)
                    .then_with(|| b.pathways.len().cmp(&a.pathways.len()))
                    .then_with(|| {
                        b.log2_fc
                            .unwrap_or(0.0)
                            .abs()
                            .total_cmp(&a.log2_fc.unwrap_or(0.0).abs())
                    })
                    .then_with(|| a.symbol.cmp(&b.symbol))
            });

            if enter_pressed {
                if let Some(top) = matches.first() {
                    actions.push(GuiAction::JumpToGene {
                        symbol: top.symbol.clone(),
                        pathway_id: top.pathways.first().map(|p| p.0.to_string()),
                    });
                }
            } else if self.omnibox_dropdown_open && !matches.is_empty() {
                let popup_pos = Pos2::new(edit_rect.left(), edit_rect.bottom() + 4.0);
                let mut close_dropdown = false;
                egui::Area::new(egui::Id::new("causal_atlas_omnibox_dropdown"))
                    .fixed_pos(popup_pos)
                    .order(egui::Order::Foreground)
                    .show(ctx, |ui| {
                        egui::Frame::popup(ui.style()).show(ui, |ui| {
                            ui.set_min_width(640.0);
                            egui::Grid::new("omnibox_autocomplete_grid")
                                .striped(true)
                                .num_columns(5)
                                .spacing([14.0, 5.0])
                                .show(ui, |ui| {
                                    ui.strong("Symbol");
                                    ui.strong("Kind");
                                    ui.strong("Compartment");
                                    ui.strong("log2FC");
                                    ui.strong("Pathways");
                                    ui.end_row();

                                    for item in matches.iter().take(12) {
                                        if ui
                                            .selectable_label(
                                                false,
                                                RichText::new(&item.symbol).strong(),
                                            )
                                            .clicked()
                                        {
                                            actions.push(GuiAction::JumpToGene {
                                                symbol: item.symbol.clone(),
                                                pathway_id: item
                                                    .pathways
                                                    .first()
                                                    .map(|p| p.0.to_string()),
                                            });
                                            close_dropdown = true;
                                        }
                                        ui.label(item.kind);
                                        ui.label(item.compartment);
                                        if let Some(fc) = item.log2_fc {
                                            let c = if fc > 0.0 {
                                                pal.expr_up
                                            } else {
                                                pal.expr_down
                                            };
                                            ui.colored_label(c, format!("{fc:+.2}"));
                                        } else {
                                            ui.label(
                                                RichText::new("—").color(pal.text_muted),
                                            );
                                        }
                                        if item.pathways.is_empty() {
                                            ui.label(
                                                RichText::new("Cohort DEG").color(pal.text_muted),
                                            );
                                        } else {
                                            let pw_str: String = item
                                                .pathways
                                                .iter()
                                                .take(3)
                                                .map(|p| p.0)
                                                .collect::<Vec<_>>()
                                                .join(", ");
                                            ui.label(pw_str);
                                        }
                                        ui.end_row();
                                    }
                                });
                        });
                    });
                if close_dropdown {
                    self.omnibox_dropdown_open = false;
                }
            }
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
        let rep = self.reports[self.active_cohort_idx].clone();
        let active_pw_id = self.active_pathway_id();
        let is_full = self.pathway_fullscreen;
        let show_drawer = !is_full && !self.sidebar_hidden;

        let total_avail_h = ui.available_height().max(380.0);

        ui.horizontal_top(|ui| {
            if show_drawer {
                ui.allocate_ui_with_layout(
                    Vec2::new(310.0, total_avail_h),
                    egui::Layout::top_down(egui::Align::Min),
                    |drawer_ui| {
                        drawer_ui.set_max_width(310.0);
                        drawer_ui.horizontal(|ui| {
                            ui.heading("Pathways");
                            ui.add_space(6.0);
                            ui.add(
                                egui::TextEdit::singleline(&mut self.pathway_search)
                                    .hint_text("Filter...")
                                    .desired_width(175.0),
                            );
                        });

                        let needle = self.pathway_search.trim().to_lowercase();
                        drawer_ui.push_id("canonical_pw_scroll", |ui| {
                            egui::ScrollArea::vertical()
                                .max_height((total_avail_h * 0.56).max(200.0))
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
                                        let is_sel =
                                            active_pw_id.as_deref() == Some(hit.pathway_id.as_str());
                                        ui.horizontal_wrapped(|ui| {
                                            if ui
                                                .selectable_label(
                                                    is_sel,
                                                    format!(
                                                        "{} ({}/{})",
                                                        hit.pathway_name,
                                                        hit.overlap_count,
                                                        hit.pathway_size
                                                    ),
                                                )
                                                .clicked()
                                            {
                                                actions.push(GuiAction::SelectPathway(idx));
                                            }
                                        });
                                        let bar_w = ((hit.neg_log10_fdr / max_bar).clamp(0.04, 1.0)
                                            as f32)
                                            * 145.0;
                                        ui.horizontal(|ui| {
                                            let (r, p) = ui.allocate_painter(
                                                Vec2::new(bar_w, 11.0),
                                                Sense::click(),
                                            );
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
                                                    "-log10(q)={:.1} | {z_txt}",
                                                    hit.neg_log10_fdr
                                                ))
                                                .size(10.5)
                                                .color(pal.text_muted),
                                            );
                                        });
                                        ui.add_space(3.0);
                                    }
                                });
                        });

                        drawer_ui.separator();
                        if let Some(ref pw_id) = active_pw_id {
                            if let Some(hit) =
                                rep.canonical_pathways.iter().find(|p| &p.pathway_id == pw_id)
                            {
                                drawer_ui.label(
                                    RichText::new(format!(
                                        "{} [{}]",
                                        hit.pathway_name, hit.pathway_id
                                    ))
                                    .strong()
                                    .size(12.0),
                                );
                                drawer_ui.label(
                                    RichText::new(format!(
                                        "{} | Overlap {}/{} ({:.0}%) | FDR {:.2e}",
                                        hit.source_db,
                                        hit.overlap_count,
                                        hit.pathway_size,
                                        hit.overlap_ratio * 100.0,
                                        hit.fdr
                                    ))
                                    .size(11.0)
                                    .color(pal.text_muted),
                                );
                                if !hit.upregulated_genes.is_empty() {
                                    drawer_ui.colored_label(
                                        pal.expr_up,
                                        RichText::new(format!(
                                            "Up: {}",
                                            hit.upregulated_genes.join(", ")
                                        ))
                                        .size(11.0),
                                    );
                                }
                                if !hit.downregulated_genes.is_empty() {
                                    drawer_ui.colored_label(
                                        pal.expr_down,
                                        RichText::new(format!(
                                            "Down: {}",
                                            hit.downregulated_genes.join(", ")
                                        ))
                                        .size(11.0),
                                    );
                                }
                            } else if let Some(pw) =
                                self.kg.pathways.iter().find(|p| &p.id == pw_id)
                            {
                                drawer_ui.label(
                                    RichText::new(format!("{} [{}]", pw.name, pw.id))
                                        .strong()
                                        .size(12.0),
                                );
                                drawer_ui.label(
                                    RichText::new(format!(
                                        "{} | {} members (unmatched in active DEG cutoff)",
                                        pw.source_db,
                                        pw.members.len()
                                    ))
                                    .size(11.0)
                                    .color(pal.text_muted),
                                );
                            }
                        }
                    },
                );
                ui.separator();
            }

            ui.vertical(|pane_ui| {
                let Some(pw_id) = active_pw_id.clone() else {
                    pane_ui.label("No canonical pathways available.");
                    return;
                };
                let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id).cloned() else {
                    pane_ui.label("Selected pathway not found in knowledge graph.");
                    return;
                };

                let wl = compute_pathway_world_layout_cached(&pw_def, &self.kg);
                let mut world_positions: BTreeMap<String, Pos2> = wl
                    .positions
                    .iter()
                    .map(|(k, &(x, y))| (k.clone(), Pos2::new(x, y)))
                    .collect();
                if let Some(pinned_map) = self.pinned_pathway_nodes.get(&pw_def.id) {
                    for (k, &p) in pinned_map {
                        if world_positions.contains_key(k) {
                            world_positions.insert(k.clone(), p);
                        }
                    }
                }
                let pinned_count = self
                    .pinned_pathway_nodes
                    .get(&pw_def.id)
                    .map(|m| m.len())
                    .unwrap_or(0);

                // Toolbar row
                pane_ui.horizontal_wrapped(|ui| {
                    ui.heading("Pathway Atlas & MAP");
                    ui.label(
                        RichText::new(format!(
                            "[{} | {} nodes | {} edges | {:.0}×{:.0}px world | {:.0}%]",
                            pw_def.id,
                            world_positions.len(),
                            wl.edges.len(),
                            wl.world_w,
                            wl.world_h,
                            self.pathway_zoom * 100.0
                        ))
                        .size(11.5)
                        .color(pal.text_muted),
                    );
                    if ui.button("Fit All").clicked() {
                        actions.push(GuiAction::FitPathwayViewport);
                    }
                    if ui.button("Center Selected").clicked() {
                        actions.push(GuiAction::CenterSelectedPathwayNode);
                    }
                    if ui.button("Re-Layout").clicked() {
                        actions.push(GuiAction::RelayoutActivePathway);
                    }
                    if ui
                        .selectable_label(self.pathway_ego_filter, "1-Hop Ego Filter")
                        .clicked()
                    {
                        actions.push(GuiAction::TogglePathwayEgoFilter);
                    }
                    let sb_lbl = if self.sidebar_hidden {
                        "Show Sidebar ▶"
                    } else {
                        "◀ Hide Sidebar"
                    };
                    if ui.button(sb_lbl).clicked() {
                        actions.push(GuiAction::ToggleSidebarHidden);
                    }
                    let fs_lbl = if is_full {
                        "Exit Full Screen"
                    } else {
                        "Full Screen"
                    };
                    if ui.button(fs_lbl).clicked() {
                        actions.push(GuiAction::TogglePathwayFullscreen);
                    }
                    if ui.button("Export SVG").clicked() {
                        actions.push(GuiAction::ExportPathwaySvgDialog);
                    }
                    if pinned_count > 0
                        && ui.button(format!("Unpin ({pinned_count})")).clicked()
                    {
                        actions.push(GuiAction::UnpinPathwayNodes(pw_def.id.clone()));
                    }
                    if pw_def.id.starts_with("R-HSA-") && pw_def.node_coords.is_empty() {
                        if ui.button("Fetch Reactome 2D").clicked() {
                            actions.push(GuiAction::FetchReactomeDiagram(pw_def.id.clone()));
                        }
                    }
                    if !self.map_user_clamps.is_empty() && ui.button("Reset MAP Clamps").clicked() {
                        actions.push(GuiAction::ClearAllClamps);
                    }
                });

                let active = self.cohorts[self.active_cohort_idx].filter_active(
                    self.config.abs_log2fc_cutoff,
                    self.config.fdr_cutoff,
                    &self.kg,
                );
                let map_res = simulate_map_in_kg(
                    &pw_def,
                    &active,
                    &self.kg,
                    &self.map_user_clamps,
                    self.config.map_damping,
                    self.config.map_iterations,
                );

                // 1-Hop Ego Filter computation
                let ego_visible: Option<BTreeSet<String>> = if self.pathway_ego_filter {
                    let root = self
                        .selected_pathway_node
                        .as_ref()
                        .filter(|s| world_positions.contains_key(*s))
                        .cloned()
                        .or_else(|| world_positions.keys().next().cloned());
                    root.map(|r| {
                        let mut set = BTreeSet::new();
                        set.insert(r.clone());
                        for e in &wl.edges {
                            if e.from == r {
                                set.insert(e.to.clone());
                            } else if e.to == r {
                                set.insert(e.from.clone());
                            }
                        }
                        set
                    })
                } else {
                    None
                };

                // Selected node MAP clamp bar
                pane_ui.horizontal_wrapped(|ui| {
                    if let Some(ref sel_nid) = self.selected_pathway_node {
                        ui.label(
                            RichText::new(format!("Selected: {sel_nid}"))
                                .strong()
                                .color(Color32::from_rgb(250, 204, 21)),
                        );
                        if ui.small_button("Activate (+1)").clicked() {
                            actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), Some(1.0)));
                        }
                        if ui.small_button("Inhibit (-1)").clicked() {
                            actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), Some(-1.0)));
                        }
                        if ui.small_button("Unclamp").clicked() {
                            actions.push(GuiAction::SetNodeClamp(sel_nid.clone(), None));
                        }
                    } else {
                        ui.label(
                            RichText::new(
                                "Scroll to zoom toward cursor (15%-350%), drag background to pan, drag any molecule to pin, click to inspect & clamp MAP.",
                            )
                            .size(11.5)
                            .color(pal.text_muted),
                        );
                    }
                    if let Some(ref ego) = ego_visible {
                        ui.separator();
                        ui.label(
                            RichText::new(format!("1-Hop Ego Subgraph: {} nodes", ego.len()))
                                .size(11.5)
                                .color(pal.activated_strong),
                        );
                    }
                });

                let canvas_sz = Vec2::new(
                    pane_ui.available_width().max(320.0),
                    (pane_ui.available_height() - 6.0).max(340.0),
                );
                let (resp, painter) =
                    pane_ui.allocate_painter(canvas_sz, Sense::click_and_drag());
                let rect = resp.rect;
                let painter = painter.with_clip_rect(rect);
                painter.rect_filled(rect, 6.0, pal.bg_canvas);
                painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

                let world_center = Pos2::new(wl.world_w * 0.5, wl.world_h * 0.5);

                if self.pathway_fit_pending {
                    self.pathway_zoom = ((rect.width() - 36.0) / wl.world_w)
                        .min((rect.height() - 36.0) / wl.world_h)
                        .clamp(0.15, 3.5);
                    self.pathway_pan = Vec2::ZERO;
                    self.pathway_fit_pending = false;
                    self.camera_anim_progress = 1.0;
                }

                if let Some(target_nid) = self.camera_pan_pending_node.take() {
                    if let Some(&wpos) = world_positions.get(&target_nid) {
                        self.camera_anim_start_pan = self.pathway_pan;
                        self.camera_anim_start_zoom = self.pathway_zoom;
                        let target_zoom = self.pathway_zoom.max(1.18).clamp(0.15, 3.5);
                        self.camera_anim_target_zoom = target_zoom;
                        self.camera_anim_target_pan = Vec2::new(
                            (world_center.x - wpos.x) * target_zoom,
                            (world_center.y - wpos.y) * target_zoom,
                        );
                        self.camera_anim_progress = 0.0;
                    }
                }

                if self.camera_anim_progress < 1.0 {
                    let dt = pane_ui.input(|i| i.unstable_dt).clamp(0.008, 0.05);
                    self.camera_anim_progress = (self.camera_anim_progress + dt * 3.8).min(1.0);
                    let ease = 1.0 - (1.0 - self.camera_anim_progress).powi(3);
                    self.pathway_zoom = self.camera_anim_start_zoom
                        + (self.camera_anim_target_zoom - self.camera_anim_start_zoom) * ease;
                    self.pathway_pan = self.camera_anim_start_pan
                        + (self.camera_anim_target_pan - self.camera_anim_start_pan) * ease;
                }

                // Cursor-anchored scroll zoom in 0.15..=3.5
                if resp.hovered() {
                    let scroll_y = pane_ui.input(|i| i.smooth_scroll_delta.y);
                    if scroll_y.abs() > 0.1 {
                        self.camera_anim_progress = 1.0;
                        let old_zoom = self.pathway_zoom;
                        let factor = (1.0 + scroll_y * 0.0022).clamp(0.80, 1.25);
                        let new_zoom = (old_zoom * factor).clamp(0.15, 3.5);
                        if let Some(cursor) = pane_ui.input(|i| i.pointer.hover_pos()) {
                            let w_under = Pos2::new(
                                world_center.x
                                    + (cursor.x - rect.center().x - self.pathway_pan.x)
                                        / old_zoom.max(0.05),
                                world_center.y
                                    + (cursor.y - rect.center().y - self.pathway_pan.y)
                                        / old_zoom.max(0.05),
                            );
                            self.pathway_pan = Vec2::new(
                                cursor.x - rect.center().x - (w_under.x - world_center.x) * new_zoom,
                                cursor.y - rect.center().y - (w_under.y - world_center.y) * new_zoom,
                            );
                        }
                        self.pathway_zoom = new_zoom;
                    }
                }

                let minimap_rect = Rect::from_min_size(
                    Pos2::new(rect.right() - 162.0, rect.bottom() - 122.0),
                    Vec2::new(150.0, 110.0),
                );
                let minimap_inner = minimap_rect.shrink(6.0);

                // Handle Mini-Map interaction OR node drag-to-pin OR background pan
                let pointer_pos = resp
                    .interact_pointer_pos()
                    .or_else(|| pane_ui.input(|i| i.pointer.hover_pos()));

                let zoom_now = self.pathway_zoom;
                let pan_now = self.pathway_pan;
                let world_to_screen = |wx: f32, wy: f32| -> Pos2 {
                    Pos2::new(
                        rect.center().x + pan_now.x + (wx - world_center.x) * zoom_now,
                        rect.center().y + pan_now.y + (wy - world_center.y) * zoom_now,
                    )
                };
                let screen_to_world = |sp: Pos2| -> Pos2 {
                    Pos2::new(
                        world_center.x + (sp.x - rect.center().x - pan_now.x) / zoom_now.max(0.05),
                        world_center.y + (sp.y - rect.center().y - pan_now.y) / zoom_now.max(0.05),
                    )
                };

                let pointer_in_minimap = pointer_pos.is_some_and(|pp| minimap_rect.contains(pp));

                if resp.drag_started() && !pointer_in_minimap {
                    self.camera_anim_progress = 1.0;
                    if let Some(pp) = pointer_pos {
                        let hit_nid = world_positions.iter().find_map(|(nid, &wpos)| {
                            if let Some(ref ego) = ego_visible {
                                if !ego.contains(nid) {
                                    return None;
                                }
                            }
                            let sp = world_to_screen(wpos.x, wpos.y);
                            let hw = (36.0 * zoom_now).max(16.0);
                            let hh = (18.0 * zoom_now).max(10.0);
                            if (pp.x - sp.x).abs() <= hw && (pp.y - sp.y).abs() <= hh {
                                Some(nid.clone())
                            } else {
                                None
                            }
                        });
                        self.dragging_pathway_node = hit_nid.clone();
                        if let Some(nid) = hit_nid {
                            self.selected_pathway_node = Some(nid.clone());
                            self.omnibox_selected_gene = Some(nid);
                        }
                    }
                }

                if resp.dragged() {
                    self.camera_anim_progress = 1.0;
                    if pointer_in_minimap {
                        if let Some(pp) = pointer_pos {
                            let u = ((pp.x - minimap_inner.left()) / minimap_inner.width().max(1.0))
                                .clamp(0.0, 1.0);
                            let v = ((pp.y - minimap_inner.top()) / minimap_inner.height().max(1.0))
                                .clamp(0.0, 1.0);
                            let wx = u * wl.world_w;
                            let wy = v * wl.world_h;
                            self.pathway_pan = Vec2::new(
                                (world_center.x - wx) * self.pathway_zoom,
                                (world_center.y - wy) * self.pathway_zoom,
                            );
                        }
                    } else if let Some(ref drag_nid) = self.dragging_pathway_node.clone() {
                        let delta_world = resp.drag_delta() / self.pathway_zoom.max(0.05);
                        if let Some(&cur_w) = world_positions.get(drag_nid) {
                            let new_w = Pos2::new(
                                (cur_w.x + delta_world.x).clamp(40.0, wl.world_w - 40.0),
                                (cur_w.y + delta_world.y).clamp(32.0, wl.world_h - 32.0),
                            );
                            world_positions.insert(drag_nid.clone(), new_w);
                            self.pinned_pathway_nodes
                                .entry(pw_def.id.clone())
                                .or_default()
                                .insert(drag_nid.clone(), new_w);
                        }
                    } else {
                        self.pathway_pan += resp.drag_delta();
                    }
                }

                if resp.drag_stopped() {
                    self.dragging_pathway_node = None;
                }

                if resp.clicked() && pointer_in_minimap {
                    if let Some(pp) = resp.interact_pointer_pos() {
                        let u = ((pp.x - minimap_inner.left()) / minimap_inner.width().max(1.0))
                            .clamp(0.0, 1.0);
                        let v = ((pp.y - minimap_inner.top()) / minimap_inner.height().max(1.0))
                            .clamp(0.0, 1.0);
                        let wx = u * wl.world_w;
                        let wy = v * wl.world_h;
                        self.pathway_pan = Vec2::new(
                            (world_center.x - wx) * self.pathway_zoom,
                            (world_center.y - wy) * self.pathway_zoom,
                        );
                    }
                }

                // Refresh closures after any pan/pin update
                let zoom = self.pathway_zoom;
                let pan = self.pathway_pan;
                let to_screen = |wx: f32, wy: f32| -> Pos2 {
                    Pos2::new(
                        rect.center().x + pan.x + (wx - world_center.x) * zoom,
                        rect.center().y + pan.y + (wy - world_center.y) * zoom,
                    )
                };

                // 1. 5-Compartment Subcellular Architectural Bands in Virtual World Space
                let pm_l = to_screen(24.0, wl.bands.plasma_membrane_y);
                let pm_r = to_screen(wl.world_w - 24.0, wl.bands.plasma_membrane_y);
                painter.line_segment(
                    [pm_l, pm_r],
                    Stroke::new((6.0 * zoom).clamp(2.0, 8.0), pal.border_default),
                );

                let org_l = to_screen(48.0, wl.bands.organelle_divider_y);
                let org_r = to_screen(wl.world_w - 48.0, wl.bands.organelle_divider_y);
                painter.line_segment([org_l, org_r], Stroke::new(1.2, pal.grid_line));

                let nm_l = to_screen(52.0, wl.bands.nuclear_envelope_y);
                let nm_r = to_screen(wl.world_w - 52.0, wl.bands.nuclear_envelope_y);
                painter.line_segment(
                    [nm_l, nm_r],
                    Stroke::new((5.5 * zoom).clamp(2.0, 7.5), pal.border_default),
                );

                let label_sz = (11.5 * zoom.clamp(0.75, 1.2)).max(9.5);
                painter.text(
                    Pos2::new(rect.left() + 12.0, rect.top() + 8.0),
                    Align2::LEFT_TOP,
                    SubcellularLocation::Extracellular.label(),
                    FontId::proportional(label_sz),
                    pal.text_muted,
                );
                painter.text(
                    Pos2::new(pm_l.x.max(rect.left() + 12.0), pm_l.y + 14.0 * zoom),
                    Align2::LEFT_TOP,
                    SubcellularLocation::Cytoplasm.label(),
                    FontId::proportional(label_sz),
                    pal.text_muted,
                );
                painter.text(
                    Pos2::new(org_l.x.max(rect.left() + 12.0), org_l.y + 6.0 * zoom),
                    Align2::LEFT_TOP,
                    SubcellularLocation::Organelle.label(),
                    FontId::proportional(label_sz * 0.92),
                    pal.text_muted,
                );
                painter.text(
                    Pos2::new(nm_l.x.max(rect.left() + 12.0), nm_l.y + 10.0 * zoom),
                    Align2::LEFT_TOP,
                    SubcellularLocation::Nucleus.label(),
                    FontId::proportional(label_sz),
                    pal.text_muted,
                );
                if !pw_def.cell_context_label.is_empty() {
                    painter.text(
                        Pos2::new(rect.right() - 12.0, rect.top() + 8.0),
                        Align2::RIGHT_TOP,
                        &pw_def.cell_context_label,
                        FontId::proportional(label_sz),
                        pal.text_muted,
                    );
                }

                // 2. Nuclear DNA Double Helix across promoter band
                let dna_l = to_screen(wl.world_w * 0.11, wl.bands.dna_helix_y);
                let dna_r = to_screen(wl.world_w * 0.89, wl.bands.dna_helix_y);
                draw_dna_double_helix(
                    &painter,
                    dna_l.x,
                    dna_r.x,
                    dna_l.y,
                    5.5 * zoom.clamp(0.4, 1.4),
                    pal.grid_line,
                );

                // 3. Draw cubic Bezier curved causal edges
                let mut pair_counts: HashMap<(String, String), usize> = HashMap::new();
                let mut edge_seen_rev: HashSet<(String, String)> = HashSet::new();
                for e in &wl.edges {
                    edge_seen_rev.insert((e.from.clone(), e.to.clone()));
                }

                for (e_idx, edge) in wl.edges.iter().enumerate() {
                    if let Some(ref ego) = ego_visible {
                        if !ego.contains(&edge.from) || !ego.contains(&edge.to) {
                            continue;
                        }
                    }
                    if let (Some(&wp1), Some(&wp2)) =
                        (world_positions.get(&edge.from), world_positions.get(&edge.to))
                    {
                        let p1 = to_screen(wp1.x, wp1.y);
                        let p2 = to_screen(wp2.x, wp2.y);
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
                        let key = (edge.from.clone(), edge.to.clone());
                        let dup_i = pair_counts.entry(key).or_insert(0);
                        let has_rev =
                            edge_seen_rev.contains(&(edge.to.clone(), edge.from.clone()));
                        let bow = if has_rev {
                            22.0 * zoom
                        } else {
                            (12.0 + (*dup_i as f32) * 14.0) * zoom
                        };
                        *dup_i += 1;

                        draw_bezier_causal_edge(
                            &painter,
                            p1,
                            p2,
                            (20.0 * zoom).clamp(8.0, 36.0),
                            (24.0 * zoom).clamp(10.0, 42.0),
                            bow,
                            edge.sign,
                            Stroke::new((1.65 * zoom).clamp(0.9, 2.2), edge_col),
                        );
                    }
                }

                // 4. Draw entity glyphs, phosphorylation badges, pin indicators & gold pulse halo
                let click_pt = if resp.clicked() && !pointer_in_minimap {
                    resp.interact_pointer_pos()
                } else {
                    None
                };
                let hover_pt = if pointer_in_minimap { None } else { resp.hover_pos() };
                let pulse_phase = pane_ui.input(|i| i.time as f32);

                for (nid, st) in &map_res.states {
                    if let Some(ref ego) = ego_visible {
                        if !ego.contains(nid) {
                            continue;
                        }
                    }
                    if let Some(&wpos) = world_positions.get(nid) {
                        let pos = to_screen(wpos.x, wpos.y);
                        let kind = self
                            .kg
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

                        let half_ext = Vec2::new((32.0 * zoom).max(16.0), (16.0 * zoom).max(9.0));
                        let is_sel = self.selected_pathway_node.as_deref() == Some(nid.as_str());
                        let is_pulsed = self.highlight_pulse_gene.as_deref() == Some(nid.as_str())
                            && self.highlight_pulse_timer > 0.0;

                        if is_pulsed || is_sel {
                            draw_gold_pulse_halo(&painter, pos, half_ext, pulse_phase);
                        }

                        let border = if is_sel || st.user_clamp.is_some() {
                            Stroke::new(2.4, Color32::from_rgb(250, 204, 21))
                        } else if st.observed_log2fc.is_some() {
                            Stroke::new(1.8, pal.text_primary)
                        } else {
                            Stroke::new(1.2, pal.border_default)
                        };

                        draw_entity_glyph(&painter, pos, half_ext, kind, fill, border);

                        if pw_def.phospho_nodes.contains(nid) {
                            draw_phospho_badge(
                                &painter,
                                Pos2::new(pos.x + half_ext.x * 0.80, pos.y - half_ext.y * 0.85),
                                pal.bg_canvas,
                                pal.border_default,
                            );
                        }
                        if self
                            .pinned_pathway_nodes
                            .get(&pw_def.id)
                            .is_some_and(|m| m.contains_key(nid))
                        {
                            painter.circle_filled(
                                Pos2::new(pos.x - half_ext.x * 0.78, pos.y - half_ext.y * 0.78),
                                3.0 * zoom.clamp(0.7, 1.3),
                                Color32::from_rgb(56, 189, 248),
                            );
                        }

                        painter.text(
                            Pos2::new(
                                pos.x,
                                pos.y
                                    + if kind == EntityKind::Kinase {
                                        -2.0 * zoom
                                    } else {
                                        0.0
                                    },
                            ),
                            Align2::CENTER_CENTER,
                            nid,
                            FontId::proportional((10.5 * zoom).clamp(8.0, 14.5)),
                            pal.text_primary,
                        );

                        if let Some(cp) = click_pt {
                            if (cp.x - pos.x).abs() <= half_ext.x
                                && (cp.y - pos.y).abs() <= half_ext.y
                            {
                                actions.push(GuiAction::SelectPathwayNode(nid.clone()));
                            }
                        }
                        if let Some(hp) = hover_pt {
                            if (hp.x - pos.x).abs() <= half_ext.x
                                && (hp.y - pos.y).abs() <= half_ext.y
                            {
                                let obs_str = st
                                    .observed_log2fc
                                    .map(|v| format!("{v:+.2}"))
                                    .unwrap_or_else(|| "unmeasured".to_string());
                                resp.clone().on_hover_text(format!(
                                    "{nid} ({})\nCompartment: {}\nObserved log2FC: {obs_str}\nMAP Activity: {:+.2}",
                                    kind.label(),
                                    self.kg
                                        .nodes
                                        .get(nid)
                                        .map(|n| n.location.label())
                                        .unwrap_or("Cytoplasm"),
                                    st.predicted_activity
                                ));
                            }
                        }
                    }
                }

                // 5. Bottom-right 150x110 px Mini-Map with live viewport box
                painter.rect_filled(
                    minimap_rect,
                    5.0,
                    Color32::from_rgba_unmultiplied(12, 16, 24, 224),
                );
                painter.rect_stroke(minimap_rect, 5.0, Stroke::new(1.2, pal.border_default));

                let map_to_mini = |wx: f32, wy: f32| -> Pos2 {
                    Pos2::new(
                        minimap_inner.left()
                            + (wx / wl.world_w.max(1.0)).clamp(0.0, 1.0) * minimap_inner.width(),
                        minimap_inner.top()
                            + (wy / wl.world_h.max(1.0)).clamp(0.0, 1.0) * minimap_inner.height(),
                    )
                };

                let mm_pm_y = map_to_mini(0.0, wl.bands.plasma_membrane_y).y;
                let mm_nm_y = map_to_mini(0.0, wl.bands.nuclear_envelope_y).y;
                painter.line_segment(
                    [
                        Pos2::new(minimap_inner.left(), mm_pm_y),
                        Pos2::new(minimap_inner.right(), mm_pm_y),
                    ],
                    Stroke::new(0.8, pal.grid_line),
                );
                painter.line_segment(
                    [
                        Pos2::new(minimap_inner.left(), mm_nm_y),
                        Pos2::new(minimap_inner.right(), mm_nm_y),
                    ],
                    Stroke::new(0.8, pal.grid_line),
                );

                for edge in &wl.edges {
                    if let (Some(&wp1), Some(&wp2)) =
                        (world_positions.get(&edge.from), world_positions.get(&edge.to))
                    {
                        painter.line_segment(
                            [map_to_mini(wp1.x, wp1.y), map_to_mini(wp2.x, wp2.y)],
                            Stroke::new(0.7, Color32::from_rgba_unmultiplied(120, 136, 160, 110)),
                        );
                    }
                }

                for (nid, &wpos) in &world_positions {
                    let is_sel = self.selected_pathway_node.as_deref() == Some(nid.as_str());
                    let dot_col = if is_sel {
                        Color32::from_rgb(250, 204, 21)
                    } else if let Some(st) = map_res.states.get(nid) {
                        if let Some(fc) = st.observed_log2fc {
                            pal.color_for_log2fc(fc)
                        } else {
                            pal.color_for_activity(st.predicted_activity)
                        }
                    } else {
                        pal.neutral_node
                    };
                    painter.circle_filled(
                        map_to_mini(wpos.x, wpos.y),
                        if is_sel { 3.2 } else { 2.1 },
                        dot_col,
                    );
                }

                // Live viewport box inside Mini-Map
                let w_tl = screen_to_world(rect.min);
                let w_br = screen_to_world(rect.max);
                let v_tl = map_to_mini(w_tl.x, w_tl.y);
                let v_br = map_to_mini(w_br.x, w_br.y);
                let vp_rect = Rect::from_min_max(v_tl, v_br).intersect(minimap_inner);
                painter.rect_filled(
                    vp_rect,
                    2.0,
                    Color32::from_rgba_unmultiplied(56, 189, 248, 34),
                );
                painter.rect_stroke(
                    vp_rect,
                    2.0,
                    Stroke::new(1.4, Color32::from_rgb(56, 189, 248)),
                );
                painter.text(
                    Pos2::new(minimap_rect.left() + 5.0, minimap_rect.top() + 3.0),
                    Align2::LEFT_TOP,
                    "MINI-MAP",
                    FontId::proportional(8.5),
                    pal.text_muted,
                );
            });
        });
    }

    fn render_ura_mechanistic_tab(
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = self.reports[self.active_cohort_idx].clone();
        let sel_u_idx = self.selected_ura_idx;
        let total_avail_h = ui.available_height().max(380.0);

        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(
                Vec2::new(560.0, total_avail_h),
                egui::Layout::top_down(egui::Align::Min),
                |left_ui| {
                    left_ui.horizontal(|ui| {
                        ui.heading("Upstream Regulators (URA)");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.regulator_search)
                                .hint_text("Filter...")
                                .desired_width(150.0),
                        );
                    });
                    let needle = self.regulator_search.trim().to_uppercase();

                    left_ui.push_id("ura_table_scroll", |ui| {
                        egui::ScrollArea::vertical()
                            .max_height((total_avail_h * 0.44).max(210.0))
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
                                            if !needle.is_empty()
                                                && !h.regulator_id.contains(&needle)
                                            {
                                                continue;
                                            }
                                            let sel = sel_u_idx == idx;
                                            if ui.selectable_label(sel, &h.regulator_id).clicked()
                                            {
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
                                            ui.label(format!(
                                                "{}/{}",
                                                h.overlap_count, h.regulon_size
                                            ));
                                            ui.end_row();
                                        }
                                    });
                            });
                    });

                    if let Some(sel_u) = rep.upstream_regulators.get(sel_u_idx) {
                        left_ui.separator();
                        left_ui.label(
                            RichText::new(format!(
                                "Target Evidence for {} ({})",
                                sel_u.regulator_name,
                                sel_u.mechanistic_summary.as_deref().unwrap_or("Direct regulon")
                            ))
                            .strong(),
                        );
                        left_ui.push_id("ura_targets_scroll", |ui| {
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
                                            ui.colored_label(
                                                c,
                                                format!("{:+.2}", ev.observed_log2fc),
                                            );
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
                    }
                },
            );
            ui.separator();

            ui.vertical(|right_ui| {
                let Some(sel_u) = rep.upstream_regulators.get(sel_u_idx) else {
                    return;
                };
                right_ui.horizontal_wrapped(|ui| {
                    ui.heading(format!("Mechanistic Network — {}", sel_u.regulator_id));
                    if ui.button("Fit All").clicked() {
                        self.ura_fit_pending = true;
                        self.mechanistic_node_offsets.clear();
                    }
                });
                right_ui.label(
                    RichText::new(
                        "Radiative-forcing layered cascade connecting upstream co-regulators & targets (scroll to zoom, drag background to pan, drag nodes to pin).",
                    )
                    .size(11.5)
                    .color(pal.text_muted),
                );

                let active = self.cohorts[self.active_cohort_idx].filter_active(
                    self.config.abs_log2fc_cutoff,
                    self.config.fdr_cutoff,
                    &self.kg,
                );
                let mn = build_mechanistic_network(
                    &sel_u.regulator_id,
                    &rep.upstream_regulators,
                    &active,
                    &self.kg,
                    self.config.mechanistic_breadth,
                    self.config.mechanistic_depth,
                );

                let shown_targets: Vec<TargetGeneRef<'_>> = sel_u
                    .target_evidences
                    .iter()
                    .take(14)
                    .map(|t| TargetGeneRef {
                        sym: &t.target_gene,
                        fc: t.observed_log2fc,
                        sign: t.edge_sign,
                    })
                    .collect();

                let mut tier_nodes = vec![TieredNodeSpec {
                    id: mn.root_regulator.clone(),
                    tier: 0,
                }];
                let mut layout_edges = Vec::new();
                for r in &mn.regulator_ids {
                    if r != &mn.root_regulator {
                        tier_nodes.push(TieredNodeSpec {
                            id: r.clone(),
                            tier: 1,
                        });
                    }
                }
                for e in &mn.edges {
                    layout_edges.push((e.from_regulator.clone(), e.to_regulator.clone()));
                }
                for tg in &shown_targets {
                    if !tier_nodes.iter().any(|n| n.id == tg.sym) {
                        tier_nodes.push(TieredNodeSpec {
                            id: tg.sym.to_string(),
                            tier: 2,
                        });
                    }
                    layout_edges.push((mn.root_regulator.clone(), tg.sym.to_string()));
                }

                let tiered = compute_tiered_causal_layout(
                    &format!("ura_{}_{}", self.active_cohort_idx, sel_u.regulator_id),
                    &tier_nodes,
                    &layout_edges,
                    MIN_NODE_DX,
                    MIN_NODE_DY,
                );

                let sz = Vec2::new(
                    right_ui.available_width().max(280.0),
                    (right_ui.available_height() - 8.0).max(320.0),
                );
                let (resp, painter) = right_ui.allocate_painter(sz, Sense::click_and_drag());
                let rect = resp.rect;
                let painter = painter.with_clip_rect(rect);
                painter.rect_filled(rect, 6.0, pal.bg_canvas);
                painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

                let world_center = Pos2::new(tiered.world_w * 0.5, tiered.world_h * 0.5);
                if self.ura_fit_pending {
                    self.ura_zoom = ((rect.width() - 32.0) / tiered.world_w)
                        .min((rect.height() - 32.0) / tiered.world_h)
                        .clamp(0.15, 3.5);
                    self.ura_pan = Vec2::ZERO;
                    self.ura_fit_pending = false;
                }
                if resp.hovered() {
                    let scroll_y = right_ui.input(|i| i.smooth_scroll_delta.y);
                    if scroll_y.abs() > 0.1 {
                        let old_zoom = self.ura_zoom;
                        let factor = (1.0 + scroll_y * 0.0022).clamp(0.80, 1.25);
                        let new_zoom = (old_zoom * factor).clamp(0.15, 3.5);
                        if let Some(cursor) = right_ui.input(|i| i.pointer.hover_pos()) {
                            let w_under = Pos2::new(
                                world_center.x
                                    + (cursor.x - rect.center().x - self.ura_pan.x)
                                        / old_zoom.max(0.05),
                                world_center.y
                                    + (cursor.y - rect.center().y - self.ura_pan.y)
                                        / old_zoom.max(0.05),
                            );
                            self.ura_pan = Vec2::new(
                                cursor.x - rect.center().x - (w_under.x - world_center.x) * new_zoom,
                                cursor.y - rect.center().y - (w_under.y - world_center.y) * new_zoom,
                            );
                        }
                        self.ura_zoom = new_zoom;
                    }
                }

                let zoom = self.ura_zoom;
                let pan = self.ura_pan;
                let mut node_pos: BTreeMap<String, Pos2> = BTreeMap::new();
                for (nid, &(wx, wy)) in &tiered.positions {
                    let off = self
                        .mechanistic_node_offsets
                        .get(nid)
                        .copied()
                        .unwrap_or(Vec2::ZERO);
                    let sx = rect.center().x + pan.x + (wx + off.x - world_center.x) * zoom;
                    let sy = rect.center().y + pan.y + (wy + off.y - world_center.y) * zoom;
                    node_pos.insert(nid.clone(), Pos2::new(sx, sy));
                }

                if resp.drag_started() {
                    if let Some(ptr) = resp.interact_pointer_pos() {
                        self.dragging_ura_node = node_pos
                            .iter()
                            .find(|(_, &p)| p.distance(ptr) <= (32.0 * zoom).max(18.0))
                            .map(|(k, _)| k.clone());
                    }
                }
                if resp.dragged() {
                    if let Some(ref drag_id) = self.dragging_ura_node {
                        let entry = self
                            .mechanistic_node_offsets
                            .entry(drag_id.clone())
                            .or_insert(Vec2::ZERO);
                        *entry += resp.drag_delta() / self.ura_zoom.max(0.05);
                    } else {
                        self.ura_pan += resp.drag_delta();
                    }
                }
                if resp.drag_stopped() {
                    self.dragging_ura_node = None;
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
                        draw_bezier_causal_edge(
                            &painter,
                            p1,
                            p2,
                            (22.0 * zoom).clamp(10.0, 32.0),
                            (24.0 * zoom).clamp(12.0, 34.0),
                            18.0 * zoom,
                            e.sign,
                            Stroke::new((2.2 * zoom).clamp(1.2, 2.8), c),
                        );
                    }
                }

                if let Some(&rp) = node_pos.get(&mn.root_regulator) {
                    for (idx, tg) in shown_targets.iter().enumerate() {
                        if let Some(&tp) = node_pos.get(tg.sym) {
                            let bow = ((idx as f32) - (shown_targets.len() as f32) * 0.5) * 4.5 * zoom;
                            draw_bezier_causal_edge(
                                &painter,
                                rp,
                                tp,
                                (18.0 * zoom).clamp(8.0, 28.0),
                                (20.0 * zoom).clamp(10.0, 30.0),
                                bow,
                                tg.sign,
                                Stroke::new((1.25 * zoom).clamp(0.8, 1.8), pal.grid_line),
                            );
                        }
                    }
                }

                for rid in &mn.regulator_ids {
                    if let Some(&p) = node_pos.get(rid) {
                        let z = mn.regulator_z_scores.get(rid).copied().unwrap_or(0.0);
                        let kind = self
                            .kg
                            .nodes
                            .get(rid)
                            .map(|n| n.kind)
                            .unwrap_or(EntityKind::TranscriptionRegulator);
                        draw_entity_glyph(
                            &painter,
                            p,
                            Vec2::new((34.0 * zoom).max(20.0), (17.0 * zoom).max(11.0)),
                            kind,
                            pal.color_for_z(z),
                            Stroke::new(2.0, pal.text_primary),
                        );
                        painter.text(
                            p,
                            Align2::CENTER_CENTER,
                            rid,
                            FontId::proportional((11.0 * zoom).clamp(8.5, 14.0)),
                            pal.text_primary,
                        );
                    }
                }

                for tg in &shown_targets {
                    if let Some(&p) = node_pos.get(tg.sym) {
                        draw_entity_glyph(
                            &painter,
                            p,
                            Vec2::new((28.0 * zoom).max(18.0), (14.0 * zoom).max(10.0)),
                            EntityKind::GeneOther,
                            pal.color_for_log2fc(tg.fc),
                            Stroke::new(1.3, pal.border_default),
                        );
                        painter.text(
                            p,
                            Align2::CENTER_CENTER,
                            tg.sym,
                            FontId::proportional((10.0 * zoom).clamp(8.0, 13.0)),
                            pal.text_primary,
                        );
                    }
                }
            });
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
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = self.reports[self.active_cohort_idx].clone();
        ui.horizontal_wrapped(|ui| {
            ui.heading("Regulator Effects — 3-Tier Causal Mechanism Explorer");
            if ui.button("Fit All").clicked() {
                self.re_fit_pending = true;
                self.re_node_offsets.clear();
            }
        });
        ui.horizontal_wrapped(|ui| {
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
            let mut tier_nodes = Vec::new();
            let mut layout_edges = Vec::new();

            for (rid, _) in &hit.regulators {
                tier_nodes.push(TieredNodeSpec {
                    id: rid.clone(),
                    tier: 0,
                });
            }
            for (gid, _) in hit.shared_dataset_targets.iter().take(24) {
                if !tier_nodes.iter().any(|n| n.id == *gid) {
                    tier_nodes.push(TieredNodeSpec {
                        id: gid.clone(),
                        tier: 1,
                    });
                }
                for (rid, _) in &hit.regulators {
                    layout_edges.push((rid.clone(), gid.clone()));
                }
            }
            for (_pid, pname, _) in &hit.downstream_phenotypes {
                let key = format!("PHENO::{pname}");
                tier_nodes.push(TieredNodeSpec {
                    id: key.clone(),
                    tier: 2,
                });
                for (gid, _) in hit.shared_dataset_targets.iter().take(24) {
                    layout_edges.push((gid.clone(), key.clone()));
                }
            }

            let tiered = compute_tiered_causal_layout(
                &format!("reg_eff_{}_{}", self.active_cohort_idx, self.selected_reg_effect_idx),
                &tier_nodes,
                &layout_edges,
                196.0,
                52.0,
            );

            let sz = Vec2::new(
                ui.available_width().max(320.0),
                (ui.available_height() - 12.0).max(340.0),
            );
            let (resp, painter) = ui.allocate_painter(sz, Sense::click_and_drag());
            let rect = resp.rect;
            let painter = painter.with_clip_rect(rect);
            painter.rect_filled(rect, 6.0, pal.bg_canvas);
            painter.rect_stroke(rect, 6.0, Stroke::new(1.0, pal.grid_line));

            let world_center = Pos2::new(tiered.world_w * 0.5, tiered.world_h * 0.5);
            if self.re_fit_pending {
                self.re_zoom = ((rect.width() - 32.0) / tiered.world_w)
                    .min((rect.height() - 32.0) / tiered.world_h)
                    .clamp(0.15, 3.5);
                self.re_pan = Vec2::ZERO;
                self.re_fit_pending = false;
            }
            if resp.hovered() {
                let scroll_y = ui.input(|i| i.smooth_scroll_delta.y);
                if scroll_y.abs() > 0.1 {
                    let old_zoom = self.re_zoom;
                    let factor = (1.0 + scroll_y * 0.0022).clamp(0.80, 1.25);
                    let new_zoom = (old_zoom * factor).clamp(0.15, 3.5);
                    if let Some(cursor) = ui.input(|i| i.pointer.hover_pos()) {
                        let w_under = Pos2::new(
                            world_center.x
                                + (cursor.x - rect.center().x - self.re_pan.x)
                                    / old_zoom.max(0.05),
                            world_center.y
                                + (cursor.y - rect.center().y - self.re_pan.y)
                                    / old_zoom.max(0.05),
                        );
                        self.re_pan = Vec2::new(
                            cursor.x - rect.center().x - (w_under.x - world_center.x) * new_zoom,
                            cursor.y - rect.center().y - (w_under.y - world_center.y) * new_zoom,
                        );
                    }
                    self.re_zoom = new_zoom;
                }
            }

            let zoom = self.re_zoom;
            let pan = self.re_pan;
            let mut node_pos: BTreeMap<String, Pos2> = BTreeMap::new();
            for (nid, &(wx, wy)) in &tiered.positions {
                let off = self.re_node_offsets.get(nid).copied().unwrap_or(Vec2::ZERO);
                let sx = rect.center().x + pan.x + (wx + off.x - world_center.x) * zoom;
                let sy = rect.center().y + pan.y + (wy + off.y - world_center.y) * zoom;
                node_pos.insert(nid.clone(), Pos2::new(sx, sy));
            }

            if resp.drag_started() {
                if let Some(ptr) = resp.interact_pointer_pos() {
                    self.dragging_re_node = node_pos
                        .iter()
                        .find(|(_, &p)| p.distance(ptr) <= (38.0 * zoom).max(18.0))
                        .map(|(k, _)| k.clone());
                }
            }
            if resp.dragged() {
                if let Some(ref drag_id) = self.dragging_re_node {
                    let entry = self
                        .re_node_offsets
                        .entry(drag_id.clone())
                        .or_insert(Vec2::ZERO);
                    *entry += resp.drag_delta() / self.re_zoom.max(0.05);
                } else {
                    self.re_pan += resp.drag_delta();
                }
            }
            if resp.drag_stopped() {
                self.dragging_re_node = None;
            }

            for (rid, z) in &hit.regulators {
                let Some(&rp) = node_pos.get(rid) else {
                    continue;
                };
                let s = if *z >= 0.0 {
                    EdgeSign::Activating
                } else {
                    EdgeSign::Inhibiting
                };
                for (gid, _) in hit.shared_dataset_targets.iter().take(24) {
                    if let Some(&gp) = node_pos.get(gid) {
                        draw_bezier_causal_edge(
                            &painter,
                            rp,
                            gp,
                            (18.0 * zoom).clamp(8.0, 28.0),
                            (20.0 * zoom).clamp(10.0, 30.0),
                            14.0 * zoom,
                            s,
                            Stroke::new((1.5 * zoom).clamp(0.9, 2.2), pal.color_for_z(*z)),
                        );
                    }
                }
            }

            for (gid, _) in hit.shared_dataset_targets.iter().take(24) {
                let Some(&gp) = node_pos.get(gid) else {
                    continue;
                };
                for (_pid, pname, pz) in &hit.downstream_phenotypes {
                    let key = format!("PHENO::{pname}");
                    if let Some(&pp) = node_pos.get(&key) {
                        let s = if *pz >= 0.0 {
                            EdgeSign::Activating
                        } else {
                            EdgeSign::Inhibiting
                        };
                        draw_bezier_causal_edge(
                            &painter,
                            gp,
                            pp,
                            (16.0 * zoom).clamp(8.0, 26.0),
                            (24.0 * zoom).clamp(10.0, 34.0),
                            -14.0 * zoom,
                            s,
                            Stroke::new((1.5 * zoom).clamp(0.9, 2.2), pal.color_for_z(*pz)),
                        );
                    }
                }
            }

            for (rid, z) in &hit.regulators {
                if let Some(&p) = node_pos.get(rid) {
                    draw_entity_glyph(
                        &painter,
                        p,
                        Vec2::new((42.0 * zoom).max(24.0), (18.0 * zoom).max(11.0)),
                        EntityKind::TranscriptionRegulator,
                        pal.color_for_z(*z),
                        Stroke::new(2.0, pal.text_primary),
                    );
                    painter.text(
                        p,
                        Align2::CENTER_CENTER,
                        rid,
                        FontId::proportional((11.5 * zoom).clamp(8.5, 14.0)),
                        pal.text_primary,
                    );
                }
            }

            for (gid, fc) in hit.shared_dataset_targets.iter().take(24) {
                if let Some(&p) = node_pos.get(gid) {
                    draw_entity_glyph(
                        &painter,
                        p,
                        Vec2::new((28.0 * zoom).max(18.0), (14.0 * zoom).max(9.0)),
                        EntityKind::GeneOther,
                        pal.color_for_log2fc(*fc),
                        Stroke::new(1.5, pal.text_primary),
                    );
                    painter.text(
                        p,
                        Align2::CENTER_CENTER,
                        gid,
                        FontId::proportional((10.5 * zoom).clamp(8.0, 13.0)),
                        pal.text_primary,
                    );
                }
            }

            for (_pid, pname, pz) in &hit.downstream_phenotypes {
                let key = format!("PHENO::{pname}");
                if let Some(&p) = node_pos.get(&key) {
                    draw_entity_glyph(
                        &painter,
                        p,
                        Vec2::new((92.0 * zoom).max(46.0), (19.0 * zoom).max(12.0)),
                        EntityKind::BiologicalProcess,
                        pal.color_for_z(*pz),
                        Stroke::new(2.0, pal.text_primary),
                    );
                    painter.text(
                        p,
                        Align2::CENTER_CENTER,
                        pname,
                        FontId::proportional((10.5 * zoom).clamp(8.0, 13.0)),
                        pal.text_primary,
                    );
                }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn omnibox_jump_and_five_pillars_render_headless() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "causal_atlas_gui_test_{}.sqlite3",
            std::process::id()
        ));
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());
        assert_eq!(WorkbenchTab::ALL.len(), 5);

        app.apply_action(
            &ctx,
            GuiAction::JumpToGene {
                symbol: "VCAM1".to_string(),
                pathway_id: Some("WP445".to_string()),
            },
        );
        assert_eq!(app.active_tab, WorkbenchTab::PathwayAtlas);
        assert_eq!(app.selected_pathway_node.as_deref(), Some("VCAM1"));
        assert_eq!(app.omnibox_selected_gene.as_deref(), Some("VCAM1"));
        assert_eq!(app.highlight_pulse_gene.as_deref(), Some("VCAM1"));

        for view in [1u8, 2, 3] {
            app.configure_screenshot_view(&ctx, view);
            let input = egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                ..Default::default()
            };
            let _ = ctx.run(input, |c| app.render_ui_frame(c));
        }
        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }
}
