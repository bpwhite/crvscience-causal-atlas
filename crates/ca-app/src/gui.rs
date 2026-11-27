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
    apply_workbench_visuals, draw_bezier_causal_edge, draw_deflected_bezier_causal_edge,
    draw_dna_double_helix, draw_entity_glyph, draw_gold_pulse_halo, draw_phospho_badge,
    lerp_color, EdgeObstacle, ScientificPalette, WorkbenchPaletteMode,
};
use ca_core::{
    build_mechanistic_network, builtin_cohorts, compute_pathway_world_layout_cached,
    compute_regulator_rank_diff, compute_tiered_causal_layout, format_regulator_rank_diff_csv,
    format_regulator_rank_diff_json, parse_expression_table, recompute_and_cache_pathway_layout,
    run_comparison_analysis, run_core_analysis, seed_kb::build_seed_knowledge_graph,
    select_stratified_mechanistic_targets, simulate_map_in_kg, ComparisonMetric,
    CoreAnalysisConfig, CoreAnalysisReport, EdgeSign, EntityKind, ExpressionDataset,
    KnowledgeGraph, PredictedState, SubcellularLocation, TieredNodeSpec, MIN_NODE_DX, MIN_NODE_DY,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SortDir {
    Asc,
    Desc,
}

impl SortDir {
    fn toggle(self) -> Self {
        match self {
            Self::Asc => Self::Desc,
            Self::Desc => Self::Asc,
        }
    }

    fn glyph(self) -> &'static str {
        match self {
            Self::Asc => "▲",
            Self::Desc => "▼",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TableSort<C: Copy + Eq> {
    col: C,
    dir: SortDir,
}

impl<C: Copy + Eq> TableSort<C> {
    fn toggle_on(slot: &mut Option<Self>, col: C, default_dir: SortDir) {
        if let Some(cur) = slot {
            if cur.col == col {
                cur.dir = cur.dir.toggle();
            } else {
                *cur = Self {
                    col,
                    dir: default_dir,
                };
            }
        } else {
            *slot = Some(Self {
                col,
                dir: default_dir,
            });
        }
    }

    fn apply(self, ord: std::cmp::Ordering) -> std::cmp::Ordering {
        match self.dir {
            SortDir::Asc => ord,
            SortDir::Desc => ord.reverse(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GeneSortCol {
    Symbol,
    Log2Fc,
    PValue,
    Fdr,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GeneDirectionFilter {
    All,
    DegsOnly,
    UpOnly,
    DownOnly,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RegStateFilter {
    All,
    Activated,
    Inhibited,
}

impl RegStateFilter {
    fn matches(self, state: PredictedState) -> bool {
        match self {
            Self::All => true,
            Self::Activated => state == PredictedState::Activated,
            Self::Inhibited => state == PredictedState::Inhibited,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UraSortCol {
    Regulator,
    Type,
    State,
    ActZ,
    BiasZ,
    OverlapP,
    Targets,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CnaSortCol {
    Root,
    Depth,
    State,
    ActZ,
    OverlapP,
    NetBiasP,
}

const CUTOFF_PRESETS: [(&str, f64, f64, &str); 3] = [
    (
        "Stringent",
        1.5,
        0.01,
        "Stringent (|log2FC| >= 1.50, FDR <= 0.01): High-specificity core drivers with minimal false positives.",
    ),
    (
        "Standard",
        1.0,
        0.05,
        "Standard (|log2FC| >= 1.00, FDR <= 0.05): Balanced 2-fold expression and 5% BH-FDR default.",
    ),
    (
        "Exploratory",
        0.5,
        0.10,
        "Exploratory (|log2FC| >= 0.50, FDR <= 0.10): High-sensitivity sweep for subtle transcriptional cascades.",
    ),
];

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
    ExportGenesCsvDialog,
    ExportUraCsvDialog,
    ExportCnaCsvDialog,
    SetDiffCohortA(usize),
    SetDiffCohortB(usize),
    ExportRankDiffJsonDialog,
    ExportRankDiffCsvDialog,
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
    TogglePathwayMinimap,
    FitPathwayViewport,
    CenterSelectedPathwayNode,
    RelayoutActivePathway,
    TogglePathwayEgoFilter,
    UnpinPathwayNodes(String),
    JumpToGene {
        symbol: String,
        pathway_id: Option<String>,
    },
    ToggleAbout,
    SelectPerturbIndex(usize),
    ApplyDrugRescue {
        hit_idx: usize,
        switch_tab: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CaPaletteAction {
    TabCohortVolcano,
    TabPathwayAtlas,
    TabCausalRegulators,
    TabPhenotypeComparison,
    TabKnowledgeDatabases,
    FocusOmnibox,
    RunCoreAnalysis,
    SyncAllCohorts,
    FitPathwayCamera,
    ClearMapClamps,
    ToggleTheme,
    ExportReportJson,
    ExportPathwaySvg,
    ExportRankDiffJson,
    ExportRankDiffCsv,
    ToggleAbout,
}

struct CaPaletteCommand {
    label: &'static str,
    category: &'static str,
    hotkey: &'static str,
    keywords: &'static str,
    action: CaPaletteAction,
}

const CA_PALETTE_COMMANDS: &[CaPaletteCommand] = &[
    CaPaletteCommand {
        label: "Switch to Cohort & Volcano Table",
        category: "View",
        hotkey: "1",
        keywords: "cohort volcano deg table expression tab",
        action: CaPaletteAction::TabCohortVolcano,
    },
    CaPaletteCommand {
        label: "Switch to Pathway Atlas (MAP)",
        category: "View",
        hotkey: "2",
        keywords: "pathway atlas map diagram wikipathways reactome tab",
        action: CaPaletteAction::TabPathwayAtlas,
    },
    CaPaletteCommand {
        label: "Switch to Causal Regulators (URA / CNA)",
        category: "View",
        hotkey: "3",
        keywords: "causal regulators upstream ura cna mechanistic tab",
        action: CaPaletteAction::TabCausalRegulators,
    },
    CaPaletteCommand {
        label: "Switch to Multi-Cohort Comparison & Reversers",
        category: "View",
        hotkey: "4",
        keywords: "comparison heatmap reversers mimickers lincs phenotype tab",
        action: CaPaletteAction::TabPhenotypeComparison,
    },
    CaPaletteCommand {
        label: "Switch to Knowledge Databases (SQLite Streams)",
        category: "View",
        hotkey: "5",
        keywords: "knowledge databases streams sqlite sync tab",
        action: CaPaletteAction::TabKnowledgeDatabases,
    },
    CaPaletteCommand {
        label: "Focus Gene / Protein Omnibox",
        category: "Navigate",
        hotkey: "Ctrl+F",
        keywords: "omnibox search gene protein jump focus",
        action: CaPaletteAction::FocusOmnibox,
    },
    CaPaletteCommand {
        label: "Run Core Analysis on Active Cohort",
        category: "Analysis",
        hotkey: "",
        keywords: "run compute core analysis ura cna canonical",
        action: CaPaletteAction::RunCoreAnalysis,
    },
    CaPaletteCommand {
        label: "Sync & Recompute All Comparison Cohorts",
        category: "Analysis",
        hotkey: "",
        keywords: "sync all cohorts comparison cache",
        action: CaPaletteAction::SyncAllCohorts,
    },
    CaPaletteCommand {
        label: "Fit Pathway Atlas Camera (Fit All)",
        category: "Viewport",
        hotkey: "F / 0",
        keywords: "fit camera zoom reset pathway",
        action: CaPaletteAction::FitPathwayCamera,
    },
    CaPaletteCommand {
        label: "Clear All MAP Node Clamps",
        category: "Simulation",
        hotkey: "R",
        keywords: "clear reset map clamps drug rescue",
        action: CaPaletteAction::ClearMapClamps,
    },
    CaPaletteCommand {
        label: "Toggle Dark / Publication Light Theme",
        category: "View",
        hotkey: "",
        keywords: "theme palette dark light publication",
        action: CaPaletteAction::ToggleTheme,
    },
    CaPaletteCommand {
        label: "Export Full Core Analysis Report (JSON)",
        category: "Export",
        hotkey: "Ctrl+S",
        keywords: "export report json core analysis",
        action: CaPaletteAction::ExportReportJson,
    },
    CaPaletteCommand {
        label: "Export Active Pathway Diagram (SVG)",
        category: "Export",
        hotkey: "Ctrl+E",
        keywords: "export pathway svg vector diagram",
        action: CaPaletteAction::ExportPathwaySvg,
    },
    CaPaletteCommand {
        label: "Export Regulator Rank Diff Digest (JSON)",
        category: "Export",
        hotkey: "",
        keywords: "export regulator rank diff comparison digest json",
        action: CaPaletteAction::ExportRankDiffJson,
    },
    CaPaletteCommand {
        label: "Export Regulator Rank Diff Digest (CSV)",
        category: "Export",
        hotkey: "",
        keywords: "export regulator rank diff comparison digest csv",
        action: CaPaletteAction::ExportRankDiffCsv,
    },
    CaPaletteCommand {
        label: "Toggle About Causal Atlas Dialog",
        category: "Help",
        hotkey: "F1",
        keywords: "about help version license citation",
        action: CaPaletteAction::ToggleAbout,
    },
];

fn ca_fuzzy_palette_score(query: &str, label: &str, keywords: &str) -> Option<i32> {
    let q = query.trim().to_ascii_lowercase();
    if q.is_empty() {
        return Some(1000);
    }
    let hay = format!("{} {}", label, keywords).to_ascii_lowercase();
    if let Some(pos) = label.to_ascii_lowercase().find(&q) {
        return Some(800 - (pos as i32) * 2);
    }
    if let Some(pos) = hay.find(&q) {
        return Some(500 - (pos as i32));
    }
    let mut score = 200_i32;
    let mut hay_iter = hay.chars();
    let mut consecutive = 0_i32;
    for qch in q.chars() {
        if qch.is_whitespace() {
            continue;
        }
        let mut matched = false;
        for hch in hay_iter.by_ref() {
            if hch == qch {
                score += 12 + consecutive * 6;
                consecutive += 1;
                matched = true;
                break;
            }
            consecutive = 0;
        }
        if !matched {
            return None;
        }
    }
    Some(score)
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
    gene_dir_filter: GeneDirectionFilter,
    gene_sort: Option<TableSort<GeneSortCol>>,
    pathway_search: String,
    regulator_search: String,
    ura_state_filter: RegStateFilter,
    ura_sort: Option<TableSort<UraSortCol>>,
    cna_search: String,
    cna_state_filter: RegStateFilter,
    cna_sort: Option<TableSort<CnaSortCol>>,
    selected_pathway_idx: usize,
    override_pathway_id: Option<String>,
    selected_pathway_node: Option<String>,
    map_user_clamps: BTreeMap<String, f64>,
    active_drug_rescue_name: Option<String>,
    selected_perturb_idx: usize,
    about_open: bool,
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
    diff_cohort_a_idx: usize,
    diff_cohort_b_idx: usize,
    show_command_palette: bool,
    command_palette_query: String,
    command_palette_selected: usize,
    focus_command_palette: bool,
    pathway_fullscreen: bool,
    sidebar_hidden: bool,
    pathway_minimap_visible: bool,
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

        let mut db_fallback_warn: Option<String> = None;
        let mut db = match SqliteKnowledgeStore::open(&db_path) {
            Ok(store) => store,
            Err(err) => {
                let warn_msg = format!(
                    "[WARN] Failed opening SQLite store at {} ({err}); falling back to in-memory store",
                    db_path.display()
                );
                eprintln!("{warn_msg}");
                db_fallback_warn = Some(warn_msg);
                SqliteKnowledgeStore::open_in_memory()
                    .expect("SQLite store initialization must succeed")
            }
        };
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
        let default_diff_b = usize::from(reports.len() > 1);
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
            gene_dir_filter: GeneDirectionFilter::All,
            gene_sort: None,
            pathway_search: String::new(),
            regulator_search: String::new(),
            ura_state_filter: RegStateFilter::All,
            ura_sort: None,
            cna_search: String::new(),
            cna_state_filter: RegStateFilter::All,
            cna_sort: None,
            selected_pathway_idx: 0,
            override_pathway_id: None,
            selected_pathway_node: None,
            map_user_clamps: BTreeMap::new(),
            active_drug_rescue_name: None,
            selected_perturb_idx: 0,
            about_open: false,
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
            diff_cohort_a_idx: 0,
            diff_cohort_b_idx: default_diff_b,
            show_command_palette: false,
            command_palette_query: String::new(),
            command_palette_selected: 0,
            focus_command_palette: false,
            pathway_fullscreen: false,
            sidebar_hidden: false,
            pathway_minimap_visible: true,
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
            sync_logs: {
                let mut logs = vec![format!(
                    "Opened SQLite WAL store ({} pathways, {} signed edges, {} cached reports)",
                    db_stats.cached_pathways, db_stats.cached_edges, db_stats.cached_analysis_reports
                )];
                if let Some(ref warn) = db_fallback_warn {
                    logs.push(warn.clone());
                }
                logs
            },
            cancel_flag: Arc::new(AtomicBool::new(false)),
            status_banner: db_fallback_warn.unwrap_or_else(|| {
                format!(
                    "SQLite WAL active ({}) — {} cache hits, {} computed",
                    db_stats.db_path, hits, misses
                )
            }),
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
                    let msg = if sum.errors.is_empty() {
                        format!(
                            "Batch sync finished: {} checked, {} cache hits, {} downloaded (+{} pathways, +{} edges)",
                            sum.streams_checked,
                            sum.cache_hits,
                            sum.streams_downloaded,
                            sum.pathways_ingested,
                            sum.causal_edges_ingested
                        )
                    } else {
                        for err in &sum.errors {
                            self.sync_logs.push(format!("[ERROR] {err}"));
                        }
                        format!(
                            "Batch sync finished with {} error(s): {} ({} checked, {} cache hits, {} downloaded, +{} pathways, +{} edges)",
                            sum.errors.len(),
                            sum.errors.join("; "),
                            sum.streams_checked,
                            sum.cache_hits,
                            sum.streams_downloaded,
                            sum.pathways_ingested,
                            sum.causal_edges_ingested
                        )
                    };
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
        self.selected_perturb_idx = 0;
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
                    self.selected_perturb_idx = 0;
                    self.map_user_clamps.clear();
                    self.active_drug_rescue_name = None;
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
                                    self.selected_pathway_idx = 0;
                                    self.override_pathway_id = None;
                                    self.selected_ura_idx = 0;
                                    self.selected_cna_idx = 0;
                                    self.selected_reg_effect_idx = 0;
                                    self.selected_perturb_idx = 0;
                                    self.map_user_clamps.clear();
                                    self.active_drug_rescue_name = None;
                                    self.mechanistic_node_offsets.clear();
                                    self.re_node_offsets.clear();
                                    self.pathway_fit_pending = true;
                                    self.ura_fit_pending = true;
                                    self.re_fit_pending = true;
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
                    .add_filter("JSON Report", &["json"])
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
                            .add_filter("SVG Vector Diagram", &["svg"])
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
                    .add_filter("SVG Vector Diagram", &["svg"])
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
            GuiAction::ExportGenesCsvDialog => {
                let csv_text = self.export_genes_csv_string();
                let ds_name = sanitize_filename_stem(&self.cohorts[self.active_cohort_idx].name);
                self.save_csv_dialog(&format!("{ds_name}_genes.csv"), &csv_text, "Measured Genes");
            }
            GuiAction::ExportUraCsvDialog => {
                let csv_text = self.export_ura_csv_string();
                let ds_name = sanitize_filename_stem(&self.cohorts[self.active_cohort_idx].name);
                self.save_csv_dialog(
                    &format!("{ds_name}_upstream_regulators.csv"),
                    &csv_text,
                    "Upstream Regulators",
                );
            }
            GuiAction::ExportCnaCsvDialog => {
                let csv_text = self.export_cna_csv_string();
                let ds_name = sanitize_filename_stem(&self.cohorts[self.active_cohort_idx].name);
                self.save_csv_dialog(
                    &format!("{ds_name}_causal_networks.csv"),
                    &csv_text,
                    "Causal Networks",
                );
            }
            GuiAction::SetDiffCohortA(idx) => {
                if idx < self.reports.len() {
                    self.diff_cohort_a_idx = idx;
                }
            }
            GuiAction::SetDiffCohortB(idx) => {
                if idx < self.reports.len() {
                    self.diff_cohort_b_idx = idx;
                }
            }
            GuiAction::ExportRankDiffJsonDialog => {
                if let Some((rep_a, rep_b)) = self.diff_report_pair() {
                    let digest = compute_regulator_rank_diff(rep_a, rep_b);
                    let default_name = format!(
                        "{}_vs_{}_rank_diff.json",
                        sanitize_filename_stem(&digest.cohort_a_name),
                        sanitize_filename_stem(&digest.cohort_b_name)
                    );
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("JSON Rank Diff Digest", &["json"])
                        .set_file_name(&default_name)
                        .save_file()
                    {
                        match format_regulator_rank_diff_json(&digest)
                            .and_then(|txt| std::fs::write(&path, txt).map_err(|e| e.to_string()))
                        {
                            Ok(()) => {
                                self.status_banner = format!(
                                    "Exported Regulator Rank Diff JSON ({} shared, {} flips) to {}",
                                    digest.shared_count,
                                    digest.direction_flip_count,
                                    path.display()
                                );
                            }
                            Err(e) => self.status_banner = format!("Rank diff JSON export error: {e}"),
                        }
                    }
                }
            }
            GuiAction::ExportRankDiffCsvDialog => {
                if let Some((rep_a, rep_b)) = self.diff_report_pair() {
                    let digest = compute_regulator_rank_diff(rep_a, rep_b);
                    let csv_text = format_regulator_rank_diff_csv(&digest);
                    let default_name = format!(
                        "{}_vs_{}_rank_diff.csv",
                        sanitize_filename_stem(&digest.cohort_a_name),
                        sanitize_filename_stem(&digest.cohort_b_name)
                    );
                    self.save_csv_dialog(&default_name, &csv_text, "Regulator Rank Diff");
                }
            }
            GuiAction::SaveSnapshotDialog => {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Compressed Knowledge Bundle", &["gz", "json"])
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
                        Ok(loaded) => match self.db.upsert_knowledge_graph(&loaded) {
                            Ok(()) => {
                                self.refresh_from_sqlite();
                                self.status_banner =
                                    format!("Merged snapshot into SQLite from {}", path.display());
                            }
                            Err(e) => {
                                self.status_banner =
                                    format!("Snapshot SQLite merge error: {e}");
                            }
                        },
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
                    if self.map_user_clamps.is_empty() {
                        self.active_drug_rescue_name = None;
                    }
                }
            }
            GuiAction::ClearAllClamps => {
                let had_clamps = !self.map_user_clamps.is_empty();
                self.map_user_clamps.clear();
                self.active_drug_rescue_name = None;
                if had_clamps {
                    self.status_banner = "Cleared all MAP node clamps".to_string();
                }
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
            GuiAction::TogglePathwayMinimap => {
                self.pathway_minimap_visible = !self.pathway_minimap_visible;
            }
            GuiAction::FitPathwayViewport => {
                self.pathway_fit_pending = true;
                self.camera_anim_progress = 1.0;
                self.camera_pan_pending_node = None;
                ctx.request_repaint();
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
                        "Recomputed 120-iteration force-directed compartment layout for {pw_id}"
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
            GuiAction::ToggleAbout => {
                self.about_open = !self.about_open;
            }
            GuiAction::SelectPerturbIndex(idx) => {
                self.selected_perturb_idx = idx;
            }
            GuiAction::ApplyDrugRescue {
                hit_idx,
                switch_tab,
            } => {
                let Some(hit) = self.reports[self.active_cohort_idx]
                    .perturbation_matches
                    .get(hit_idx)
                    .cloned()
                else {
                    return;
                };
                self.selected_perturb_idx = hit_idx;

                let active = self.cohorts[self.active_cohort_idx].filter_active(
                    self.config.abs_log2fc_cutoff,
                    self.config.fdr_cutoff,
                    &self.kg,
                );
                let mut new_clamps: BTreeMap<String, f64> = BTreeMap::new();
                for (g, sgn) in &hit.target_clamps {
                    new_clamps.insert(g.clone(), f64::from((*sgn).clamp(-1, 1)));
                }
                let is_reverser = hit.connectivity_score <= -25.0;
                for g in &hit.matched_up {
                    let fc = active.degs.get(g).map(|m| m.log2_fc).unwrap_or(0.0);
                    if is_reverser && fc < 0.0 {
                        new_clamps.entry(g.clone()).or_insert(1.0);
                    }
                }
                for g in &hit.matched_down {
                    let fc = active.degs.get(g).map(|m| m.log2_fc).unwrap_or(0.0);
                    if is_reverser && fc > 0.0 {
                        new_clamps.entry(g.clone()).or_insert(-1.0);
                    }
                }
                if new_clamps.is_empty() {
                    for g in &hit.matched_up {
                        new_clamps.entry(g.clone()).or_insert(1.0);
                    }
                    for g in &hit.matched_down {
                        new_clamps.entry(g.clone()).or_insert(-1.0);
                    }
                }

                if switch_tab {
                    let rep = &self.reports[self.active_cohort_idx];
                    let best_pw_id = self
                        .kg
                        .pathways
                        .iter()
                        .max_by_key(|pw| {
                            let has_node = |g: &str| {
                                pw.members.iter().any(|m| m == g)
                                    || pw.topology_edges.iter().any(|e| e.from == g || e.to == g)
                            };
                            let target_hits = hit
                                .target_clamps
                                .iter()
                                .filter(|(g, _)| has_node(g))
                                .count();
                            let deg_hits = hit
                                .matched_up
                                .iter()
                                .chain(hit.matched_down.iter())
                                .filter(|g| has_node(g))
                                .count();
                            let rep_bonus = rep
                                .canonical_pathways
                                .iter()
                                .position(|h| h.pathway_id == pw.id)
                                .map(|pos| 100usize.saturating_sub(pos))
                                .unwrap_or(0);
                            (target_hits * 1000 + deg_hits * 100 + rep_bonus, pw.members.len())
                        })
                        .map(|pw| pw.id.clone());

                    if let Some(pw_id) = best_pw_id {
                        if let Some(idx) = rep
                            .canonical_pathways
                            .iter()
                            .position(|h| h.pathway_id == pw_id)
                        {
                            self.selected_pathway_idx = idx;
                            self.override_pathway_id = None;
                        } else {
                            self.override_pathway_id = Some(pw_id.clone());
                        }
                        self.active_tab = WorkbenchTab::PathwayAtlas;
                        self.pathway_fit_pending = true;
                        self.camera_anim_progress = 1.0;
                        if let Some(pw_def) = self.kg.pathways.iter().find(|p| p.id == pw_id) {
                            let focus = hit
                                .target_clamps
                                .iter()
                                .map(|(g, _)| g)
                                .chain(hit.matched_down.iter())
                                .chain(hit.matched_up.iter())
                                .find(|g| pw_def.members.iter().any(|m| m == *g))
                                .cloned();
                            self.selected_pathway_node = focus;
                        }
                    }
                }

                let count = new_clamps.len();
                self.map_user_clamps = new_clamps;
                self.active_drug_rescue_name = Some(hit.perturbagen.clone());
                self.status_banner = format!(
                    "Applied Drug Rescue '{}' ({} MAP node clamps active)",
                    hit.perturbagen, count
                );
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
    #[allow(dead_code)]
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

        let pal = ScientificPalette::for_mode(self.palette_mode);
        let mut actions: Vec<GuiAction> = Vec::new();
        let wants_kb = ctx.wants_keyboard_input();

        self.process_keyboard_shortcuts(ctx, wants_kb, &mut actions);

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
                if ui.selectable_label(self.about_open, "About").clicked() {
                    actions.push(GuiAction::ToggleAbout);
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
                        RichText::new("Copyright CRV Science Network")
                            .size(11.5)
                            .color(pal.text_muted),
                    );
                    ui.separator();
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
                    ui.horizontal_wrapped(|ui| {
                        for &(name, fc_p, fdr_p, tooltip) in &CUTOFF_PRESETS {
                            let is_active = (self.config.abs_log2fc_cutoff - fc_p).abs() < 1e-4
                                && (self.config.fdr_cutoff - fdr_p).abs() < 1e-6;
                            if ui
                                .selectable_label(is_active, name)
                                .on_hover_text(tooltip)
                                .clicked()
                            {
                                self.config.abs_log2fc_cutoff = fc_p;
                                self.config.fdr_cutoff = fdr_p;
                                cutoff_changed = true;
                            }
                        }
                    });
                    cutoff_changed |= ui
                        .add(
                            egui::Slider::new(&mut self.config.abs_log2fc_cutoff, 0.25..=3.0)
                                .text("|log2FC|"),
                        )
                        .on_hover_text("Minimum absolute log2 fold-change threshold for including differentially expressed genes in causal scoring")
                        .changed();
                    cutoff_changed |= ui
                        .add(
                            egui::Slider::new(&mut self.config.fdr_cutoff, 0.0001..=0.25)
                                .logarithmic(true)
                                .text("Max FDR"),
                        )
                        .on_hover_text("Maximum Benjamini-Hochberg adjusted p-value (FDR) cutoff for significant dataset molecules")
                        .changed();
                    cutoff_changed |= ui
                        .add(
                            egui::Slider::new(&mut self.config.cna_max_depth, 1..=3)
                                .text("CNA Hop K"),
                        )
                        .on_hover_text("Maximum directed causal path length (1 to 3 hops) explored by Multi-Hop Causal Network Analysis")
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

                            ui.label("Sign Skew (mu_2):")
                                .on_hover_text("Mean signed fold-change polarity balance across active DEGs (used for Kramer bias-corrected Z-score calibration)");
                            ui.label(format!("{:+.2}", sum.mean_sign_skew))
                                .on_hover_text("Mean signed fold-change polarity balance across active DEGs (used for Kramer bias-corrected Z-score calibration)");
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
                self.render_dataset_volcano_tab(ui, &pal, &mut actions);
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

        if self.show_command_palette {
            self.render_command_palette_modal(ctx, &mut actions);
        }

        if self.about_open {
            let mut open = self.about_open;
            egui::Window::new("About Causal Atlas")
                .open(&mut open)
                .collapsible(false)
                .resizable(false)
                .anchor(Align2::CENTER_CENTER, Vec2::ZERO)
                .show(ctx, |ui| {
                    ui.set_min_width(420.0);
                    ui.heading(format!("Causal Atlas v{}", env!("CARGO_PKG_VERSION")));
                    ui.label(
                        RichText::new("Copyright CRV Science Network")
                            .strong()
                            .color(pal.activated_strong),
                    );
                    ui.add_space(4.0);
                    ui.label(
                        "Independent clean-room native Rust causal pathway, upstream regulator, mechanistic network, and transcriptomic perturbation reversal workbench backed by SQLite 3 WAL caching.",
                    );
                    ui.add_space(6.0);
                    ui.label(
                        RichText::new(
                            "Licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).",
                        )
                        .size(11.5)
                        .color(pal.text_muted),
                    );
                    ui.separator();
                    ui.hyperlink_to(
                        "CRV Science Network (https://crvscience.com)",
                        "https://crvscience.com",
                    );
                    ui.hyperlink_to(
                        "CRV Science Apps Hub (https://crvscience.com/apps)",
                        "https://crvscience.com/apps",
                    );
                    ui.hyperlink_to(
                        "Preprint APP-PREPRINT-04 (Causal-Atlas)",
                        "https://crvscience.com/apps/preprints/causal-atlas-pathway-perturbation",
                    );
                    ui.add_space(8.0);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                        if ui.button("Close").clicked() {
                            actions.push(GuiAction::ToggleAbout);
                        }
                    });
                });
            self.about_open = open;
        }

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
    fn render_dataset_volcano_tab(
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let fc_cut = self.config.abs_log2fc_cutoff;
        let q_cut = self.config.fdr_cutoff;
        let selected_gene = self.omnibox_selected_gene.clone();
        let filtered_genes: Vec<ca_core::GeneMeasurement> = self
            .filtered_sorted_genes()
            .into_iter()
            .cloned()
            .collect();

        ui.columns(2, |cols| {
            let ds = &self.cohorts[self.active_cohort_idx];
            cols[0].heading("Interactive Volcano Plot (log2FC vs. -log10 FDR)");
            let avail = cols[0].available_size();
            let plot_size = Vec2::new(avail.x.max(280.0), (avail.y - 32.0).max(280.0));
            let (resp, painter) = cols[0].allocate_painter(plot_size, Sense::click());
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, pal.bg_canvas);
            painter.rect_stroke(rect, 6.0, Stroke::new(1.0_f32, pal.grid_line));

            let max_x = ds
                .rows
                .iter()
                .map(|r| r.log2_fc.abs())
                .fold(3.5f64, f64::max)
                * 1.12;
            let max_y = ds
                .rows
                .iter()
                .map(|r| -r.fdr.max(1e-25).log10())
                .fold(6.0f64, f64::max)
                * 1.10;

            let pad_left = 54.0f32;
            let pad_right = 20.0f32;
            let pad_top = 24.0f32;
            let pad_bottom = 44.0f32;
            let inner = Rect::from_min_max(
                Pos2::new(rect.left() + pad_left, rect.top() + pad_top),
                Pos2::new(rect.right() - pad_right, rect.bottom() - pad_bottom),
            );

            let map_pt = |fc: f64, neg_log_q: f64| -> Pos2 {
                let nx = ((fc + max_x) / (2.0 * max_x)).clamp(0.0, 1.0) as f32;
                let ny = (neg_log_q / max_y).clamp(0.0, 1.0) as f32;
                Pos2::new(
                    inner.left() + nx * inner.width(),
                    inner.bottom() - ny * inner.height(),
                )
            };

            // Numeric X/Y ticks & gridlines
            let grid_stroke = Stroke::new(0.8_f32, pal.grid_line);
            let axis_stroke = Stroke::new(1.3_f32, pal.border_default);

            let x_ticks = compute_nice_ticks(-max_x, max_x, 7);
            for &xt in &x_ticks {
                let sx = map_pt(xt, 0.0).x;
                if sx >= inner.left() - 1.0 && sx <= inner.right() + 1.0 {
                    painter.line_segment(
                        [Pos2::new(sx, inner.top()), Pos2::new(sx, inner.bottom())],
                        if xt.abs() < 1e-6 {
                            Stroke::new(1.0_f32, pal.border_default)
                        } else {
                            grid_stroke
                        },
                    );
                    painter.line_segment(
                        [
                            Pos2::new(sx, inner.bottom()),
                            Pos2::new(sx, inner.bottom() + 4.0),
                        ],
                        axis_stroke,
                    );
                    let lbl = if xt.abs() < 1e-6 {
                        "0".to_string()
                    } else if (xt - xt.round()).abs() < 1e-6 {
                        format!("{xt:+.0}")
                    } else {
                        format!("{xt:+.1}")
                    };
                    painter.text(
                        Pos2::new(sx, inner.bottom() + 6.0),
                        Align2::CENTER_TOP,
                        lbl,
                        FontId::proportional(10.0),
                        pal.text_muted,
                    );
                }
            }

            let y_ticks = compute_nice_ticks(0.0, max_y, 6);
            for &yt in &y_ticks {
                if yt < -1e-6 {
                    continue;
                }
                let sy = map_pt(0.0, yt).y;
                if sy >= inner.top() - 1.0 && sy <= inner.bottom() + 1.0 {
                    painter.line_segment(
                        [Pos2::new(inner.left(), sy), Pos2::new(inner.right(), sy)],
                        grid_stroke,
                    );
                    painter.line_segment(
                        [Pos2::new(inner.left() - 4.0, sy), Pos2::new(inner.left(), sy)],
                        axis_stroke,
                    );
                    let lbl = if (yt - yt.round()).abs() < 1e-6 {
                        format!("{yt:.0}")
                    } else {
                        format!("{yt:.1}")
                    };
                    painter.text(
                        Pos2::new(inner.left() - 7.0, sy),
                        Align2::RIGHT_CENTER,
                        lbl,
                        FontId::proportional(10.0),
                        pal.text_muted,
                    );
                }
            }

            // Axis spines & titles
            painter.line_segment(
                [inner.left_top(), inner.left_bottom()],
                axis_stroke,
            );
            painter.line_segment(
                [inner.left_bottom(), inner.right_bottom()],
                axis_stroke,
            );
            painter.text(
                Pos2::new(inner.center().x, rect.bottom() - 6.0),
                Align2::CENTER_BOTTOM,
                "log2 Fold Change (log2FC)",
                FontId::proportional(11.0),
                pal.text_primary,
            );
            painter.text(
                Pos2::new(rect.left() + 10.0, rect.top() + 6.0),
                Align2::LEFT_TOP,
                "-log10(BH-FDR)",
                FontId::proportional(11.0),
                pal.text_primary,
            );

            // Cutoff guide lines & value labels
            let x_pos = map_pt(fc_cut, 0.0).x;
            let x_neg = map_pt(-fc_cut, 0.0).x;
            let y_q = map_pt(0.0, -q_cut.max(1e-25).log10()).y;
            let cut_stroke = Stroke::new(1.2_f32, pal.inconsistent_amber);
            painter.line_segment(
                [Pos2::new(x_pos, inner.top()), Pos2::new(x_pos, inner.bottom())],
                cut_stroke,
            );
            painter.line_segment(
                [Pos2::new(x_neg, inner.top()), Pos2::new(x_neg, inner.bottom())],
                cut_stroke,
            );
            painter.line_segment(
                [Pos2::new(inner.left(), y_q), Pos2::new(inner.right(), y_q)],
                cut_stroke,
            );
            painter.text(
                Pos2::new(x_pos + 3.0, inner.top() + 2.0),
                Align2::LEFT_TOP,
                format!("+{fc_cut:.2}"),
                FontId::proportional(9.5),
                pal.inconsistent_amber,
            );
            painter.text(
                Pos2::new(x_neg - 3.0, inner.top() + 2.0),
                Align2::RIGHT_TOP,
                format!("-{fc_cut:.2}"),
                FontId::proportional(9.5),
                pal.inconsistent_amber,
            );
            let fdr_lbl = if q_cut < 0.001 {
                format!("FDR ≤ {q_cut:.1e}")
            } else {
                format!("FDR ≤ {q_cut:.3}")
            };
            painter.text(
                Pos2::new(inner.right() - 4.0, (y_q - 3.0).max(inner.top() + 12.0)),
                Align2::RIGHT_BOTTOM,
                fdr_lbl,
                FontId::proportional(9.5),
                pal.inconsistent_amber,
            );

            let hover_pos = resp.hover_pos();
            let click_pos = if resp.clicked() {
                resp.interact_pointer_pos()
            } else {
                None
            };
            let mut nearest_hover: Option<(f32, &str, f64, f64)> = None;
            let mut nearest_click: Option<(f32, &str)> = None;
            let mut label_candidates: Vec<(f64, Pos2, &str)> = Vec::new();

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
                let is_sel = selected_gene.as_deref() == Some(row.symbol.as_str());
                let rad = if is_sel {
                    6.0
                } else if active {
                    4.8
                } else {
                    3.0
                };
                if is_sel {
                    painter.circle_stroke(
                        pt,
                        rad + 3.0,
                        Stroke::new(2.0_f32, Color32::from_rgb(250, 204, 21)),
                    );
                }
                painter.circle_filled(pt, rad, col);
                if active && (yv > max_y * 0.28 || is_sel) {
                    label_candidates.push((yv + row.log2_fc.abs() * 0.1, pt, row.symbol.as_str()));
                }
                if let Some(hp) = hover_pos {
                    let d = hp.distance(pt);
                    if d <= 10.0 && nearest_hover.is_none_or(|(bd, ..)| d < bd) {
                        nearest_hover = Some((d, row.symbol.as_str(), row.log2_fc, row.fdr));
                    }
                }
                if let Some(cp) = click_pos {
                    let d = cp.distance(pt);
                    if d <= 10.0 && nearest_click.is_none_or(|(bd, _)| d < bd) {
                        nearest_click = Some((d, row.symbol.as_str()));
                    }
                }
            }

            // Collision-avoidant DEG label placement sorted by significance priority
            label_candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
            let mut placed_rects: Vec<Rect> = Vec::new();
            let candidate_offsets = [
                (Vec2::new(0.0, -9.0), Align2::CENTER_BOTTOM),
                (Vec2::new(7.0, -7.0), Align2::LEFT_BOTTOM),
                (Vec2::new(-7.0, -7.0), Align2::RIGHT_BOTTOM),
                (Vec2::new(0.0, 9.0), Align2::CENTER_TOP),
                (Vec2::new(8.0, 0.0), Align2::LEFT_CENTER),
                (Vec2::new(-8.0, 0.0), Align2::RIGHT_CENTER),
            ];
            for (_score, pt, sym) in label_candidates {
                let est_w = (sym.chars().count() as f32) * 6.4 + 6.0;
                let est_h = 13.0f32;
                for &(off, align) in &candidate_offsets {
                    let anchor = pt + off;
                    let r = match align {
                        Align2::CENTER_BOTTOM => Rect::from_min_size(
                            Pos2::new(anchor.x - est_w * 0.5, anchor.y - est_h),
                            Vec2::new(est_w, est_h),
                        ),
                        Align2::LEFT_BOTTOM => Rect::from_min_size(
                            Pos2::new(anchor.x, anchor.y - est_h),
                            Vec2::new(est_w, est_h),
                        ),
                        Align2::RIGHT_BOTTOM => Rect::from_min_size(
                            Pos2::new(anchor.x - est_w, anchor.y - est_h),
                            Vec2::new(est_w, est_h),
                        ),
                        Align2::CENTER_TOP => Rect::from_min_size(
                            Pos2::new(anchor.x - est_w * 0.5, anchor.y),
                            Vec2::new(est_w, est_h),
                        ),
                        Align2::LEFT_CENTER => Rect::from_min_size(
                            Pos2::new(anchor.x, anchor.y - est_h * 0.5),
                            Vec2::new(est_w, est_h),
                        ),
                        _ => Rect::from_min_size(
                            Pos2::new(anchor.x - est_w, anchor.y - est_h * 0.5),
                            Vec2::new(est_w, est_h),
                        ),
                    };
                    let in_bounds = r.left() >= inner.left() + 2.0
                        && r.right() <= inner.right() - 2.0
                        && r.top() >= inner.top() + 2.0
                        && r.bottom() <= inner.bottom() - 2.0;
                    if in_bounds && !placed_rects.iter().any(|prev| prev.intersects(r)) {
                        painter.text(
                            anchor,
                            align,
                            sym,
                            FontId::proportional(10.0),
                            pal.text_primary,
                        );
                        placed_rects.push(r);
                        break;
                    }
                }
            }

            if let Some((_, sym)) = nearest_click {
                actions.push(GuiAction::JumpToGene {
                    symbol: sym.to_string(),
                    pathway_id: None,
                });
            }
            if let Some((_, sym, fc, q)) = nearest_hover {
                resp.on_hover_text(format!(
                    "{sym}\nlog2FC: {fc:+.3}\nFDR: {q:.2e}\nClick point to focus in Omnibox & Pathway Atlas"
                ));
            }

            cols[1].horizontal_wrapped(|ui| {
                ui.heading("Measured Gene Table");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Export CSV").clicked() {
                        actions.push(GuiAction::ExportGenesCsvDialog);
                    }
                });
            });
            cols[1].horizontal_wrapped(|ui| {
                ui.label("Filter:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.gene_search)
                        .hint_text("Symbol...")
                        .desired_width(120.0),
                );
                if !self.gene_search.is_empty() && ui.small_button("✕").clicked() {
                    self.gene_search.clear();
                }
                for (filt, lbl) in [
                    (GeneDirectionFilter::All, "All"),
                    (GeneDirectionFilter::DegsOnly, "DEGs"),
                    (GeneDirectionFilter::UpOnly, "Up (+)"),
                    (GeneDirectionFilter::DownOnly, "Down (-)"),
                ] {
                    if ui
                        .selectable_label(self.gene_dir_filter == filt, lbl)
                        .clicked()
                    {
                        self.gene_dir_filter = filt;
                    }
                }
            });
            cols[1].push_id("dataset_gene_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("gene_table_grid")
                        .striped(true)
                        .num_columns(4)
                        .spacing([14.0, 5.0])
                        .show(ui, |ui| {
                            render_sort_header(
                                ui,
                                "Symbol",
                                GeneSortCol::Symbol,
                                SortDir::Asc,
                                &mut self.gene_sort,
                            );
                            render_sort_header(
                                ui,
                                "log2FC",
                                GeneSortCol::Log2Fc,
                                SortDir::Desc,
                                &mut self.gene_sort,
                            );
                            render_sort_header(
                                ui,
                                "P-value",
                                GeneSortCol::PValue,
                                SortDir::Asc,
                                &mut self.gene_sort,
                            );
                            render_sort_header(
                                ui,
                                "BH-FDR",
                                GeneSortCol::Fdr,
                                SortDir::Asc,
                                &mut self.gene_sort,
                            );
                            ui.end_row();

                            for r in &filtered_genes {
                                let pass = r.log2_fc.abs() >= fc_cut && r.fdr <= q_cut;
                                let c = if !pass {
                                    pal.text_muted
                                } else if r.log2_fc > 0.0 {
                                    pal.expr_up
                                } else {
                                    pal.expr_down
                                };
                                let is_sel = selected_gene.as_deref() == Some(r.symbol.as_str());
                                if ui
                                    .selectable_label(is_sel, RichText::new(&r.symbol).color(c))
                                    .on_hover_text("Focus gene in Omnibox & Pathway Atlas")
                                    .clicked()
                                {
                                    actions.push(GuiAction::JumpToGene {
                                        symbol: r.symbol.clone(),
                                        pathway_id: None,
                                    });
                                }
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
                                                .map(|z| format!("Z={}", ca_core::format_z_2dp(z)))
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
                    if ui
                        .button("Fit All")
                        .on_hover_text("Fit entire pathway diagram to viewport (Shortcut: F or 0)")
                        .clicked()
                    {
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
                    let mm_lbl = if self.pathway_minimap_visible {
                        "Hide Mini-Map"
                    } else {
                        "Show Mini-Map"
                    };
                    if ui.button(mm_lbl).clicked() {
                        actions.push(GuiAction::TogglePathwayMinimap);
                    }
                    let fs_lbl = if is_full {
                        "Exit Full Screen"
                    } else {
                        "Full Screen"
                    };
                    if ui.button(fs_lbl).clicked() {
                        actions.push(GuiAction::TogglePathwayFullscreen);
                    }
                    if ui
                        .button("Export SVG")
                        .on_hover_text("Export vector MAP diagram to SVG (Shortcut: Ctrl+E)")
                        .clicked()
                    {
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
                    if !self.map_user_clamps.is_empty()
                        && ui
                            .button("Reset MAP Clamps")
                            .on_hover_text("Clear all user & drug-rescue MAP clamps (Shortcut: R)")
                            .clicked()
                    {
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

                // Selected node MAP clamp bar + quick Drug Rescue
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
                    let pw_rescue_hits: Vec<(usize, &ca_core::PerturbationHit)> = rep
                        .perturbation_matches
                        .iter()
                        .enumerate()
                        .filter(|(_, h)| h.connectivity_score <= -25.0)
                        .filter(|(_, h)| {
                            h.target_clamps
                                .iter()
                                .any(|(g, _)| world_positions.contains_key(g))
                                || h.matched_up
                                    .iter()
                                    .chain(h.matched_down.iter())
                                    .any(|g| world_positions.contains_key(g))
                        })
                        .take(3)
                        .collect();
                    if !pw_rescue_hits.is_empty() {
                        ui.separator();
                        ui.label(
                            RichText::new("Drug Rescue:")
                                .size(11.5)
                                .strong()
                                .color(pal.inhibited_strong),
                        );
                        for (h_idx, hit) in pw_rescue_hits {
                            let short_name = short_perturbagen_label(&hit.perturbagen);
                            if ui
                                .small_button(format!("{short_name} ({:+.1})", hit.connectivity_score))
                                .on_hover_text(format!(
                                    "1-click MAP Drug Rescue: {}\n{}",
                                    hit.perturbagen, hit.mechanism_of_action
                                ))
                                .clicked()
                            {
                                actions.push(GuiAction::ApplyDrugRescue {
                                    hit_idx: h_idx,
                                    switch_tab: false,
                                });
                            }
                        }
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

                // Active Clamps chip bar
                if !self.map_user_clamps.is_empty() {
                    pane_ui.horizontal_wrapped(|ui| {
                        ui.label(
                            RichText::new(format!(
                                "Active Clamps ({}):",
                                self.map_user_clamps.len()
                            ))
                            .size(11.5)
                            .strong()
                            .color(Color32::from_rgb(250, 204, 21)),
                        );
                        if let Some(ref drug) = self.active_drug_rescue_name {
                            ui.label(
                                RichText::new(format!("[{drug}]"))
                                    .size(11.0)
                                    .color(pal.inhibited_strong),
                            );
                        }
                        for (nid, &val) in &self.map_user_clamps {
                            let is_act = val >= 0.0;
                            let accent = if is_act {
                                pal.activated_strong
                            } else {
                                pal.inhibited_strong
                            };
                            let fill = lerp_color(pal.bg_panel, accent, 0.22);
                            egui::Frame::none()
                                .fill(fill)
                                .stroke(Stroke::new(1.0_f32, accent))
                                .rounding(4.0)
                                .inner_margin(egui::Margin::symmetric(5.0, 2.0))
                                .show(ui, |ui| {
                                    ui.spacing_mut().item_spacing.x = 4.0;
                                    if ui
                                        .add(
                                            egui::Button::new(
                                                RichText::new(nid).size(11.0).strong(),
                                            )
                                            .frame(false),
                                        )
                                        .on_hover_text("Focus clamped node")
                                        .clicked()
                                    {
                                        if world_positions.contains_key(nid) {
                                            actions
                                                .push(GuiAction::SelectPathwayNode(nid.clone()));
                                            actions.push(GuiAction::CenterSelectedPathwayNode);
                                        } else {
                                            actions.push(GuiAction::JumpToGene {
                                                symbol: nid.clone(),
                                                pathway_id: None,
                                            });
                                        }
                                    }
                                    let pol_lbl = if is_act { "+1" } else { "-1" };
                                    if ui
                                        .add(
                                            egui::Button::new(
                                                RichText::new(pol_lbl)
                                                    .size(10.5)
                                                    .strong()
                                                    .color(accent),
                                            )
                                            .frame(false),
                                        )
                                        .on_hover_text("Toggle clamp polarity (+1 / -1)")
                                        .clicked()
                                    {
                                        let next = if is_act { -1.0 } else { 1.0 };
                                        actions.push(GuiAction::SetNodeClamp(
                                            nid.clone(),
                                            Some(next),
                                        ));
                                    }
                                    if ui
                                        .add(
                                            egui::Button::new(
                                                RichText::new("×")
                                                    .size(11.0)
                                                    .color(pal.text_muted),
                                            )
                                            .frame(false),
                                        )
                                        .on_hover_text("Remove this clamp")
                                        .clicked()
                                    {
                                        actions.push(GuiAction::SetNodeClamp(nid.clone(), None));
                                    }
                                });
                        }
                        if ui.small_button("Clear All").clicked() {
                            actions.push(GuiAction::ClearAllClamps);
                        }
                    });
                }

                let canvas_sz = Vec2::new(
                    pane_ui.available_width().max(320.0),
                    (pane_ui.available_height() - 6.0).max(340.0),
                );
                let (resp, painter) =
                    pane_ui.allocate_painter(canvas_sz, Sense::click_and_drag());
                let rect = resp.rect;
                let painter = painter.with_clip_rect(rect);
                painter.rect_filled(rect, 6.0, pal.bg_canvas);
                painter.rect_stroke(rect, 6.0, Stroke::new(1.0_f32, pal.grid_line));

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

                let minimap_header_h = 18.0;
                let minimap_rect = Rect::from_min_size(
                    Pos2::new(rect.right() - 168.0, rect.bottom() - 140.0),
                    Vec2::new(156.0, 128.0),
                );
                let minimap_close_rect = Rect::from_min_size(
                    Pos2::new(minimap_rect.right() - 18.0, minimap_rect.top() + 1.0),
                    Vec2::new(16.0, minimap_header_h - 2.0),
                );
                let minimap_restore_rect = Rect::from_min_size(
                    Pos2::new(rect.right() - 104.0, rect.bottom() - 32.0),
                    Vec2::new(92.0, 22.0),
                );
                let minimap_inner = Rect::from_min_max(
                    Pos2::new(minimap_rect.left() + 6.0, minimap_rect.top() + minimap_header_h + 4.0),
                    Pos2::new(minimap_rect.right() - 6.0, minimap_rect.bottom() - 6.0),
                );

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

                let pointer_in_minimap_close =
                    self.pathway_minimap_visible && pointer_pos.is_some_and(|pp| minimap_close_rect.contains(pp));
                let pointer_in_minimap_restore =
                    !self.pathway_minimap_visible && pointer_pos.is_some_and(|pp| minimap_restore_rect.contains(pp));
                let pointer_in_minimap = (self.pathway_minimap_visible
                    && pointer_pos.is_some_and(|pp| minimap_rect.contains(pp)))
                    || pointer_in_minimap_restore;

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
                    if pointer_in_minimap_close || pointer_in_minimap_restore {
                        actions.push(GuiAction::TogglePathwayMinimap);
                    } else if let Some(pp) = resp.interact_pointer_pos() {
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

                // 1. Subtle 5-Compartment Subcellular Architectural Band Fills & Dividers
                let pm_l = to_screen(24.0, wl.bands.plasma_membrane_y);
                let pm_r = to_screen(wl.world_w - 24.0, wl.bands.plasma_membrane_y);
                let org_l = to_screen(48.0, wl.bands.organelle_divider_y);
                let org_r = to_screen(wl.world_w - 48.0, wl.bands.organelle_divider_y);
                let nm_l = to_screen(52.0, wl.bands.nuclear_envelope_y);
                let nm_r = to_screen(wl.world_w - 52.0, wl.bands.nuclear_envelope_y);

                draw_pathway_compartment_bands(
                    &painter,
                    rect,
                    pm_l.y,
                    org_l.y,
                    nm_l.y,
                    zoom,
                    self.palette_mode,
                );

                painter.line_segment(
                    [pm_l, pm_r],
                    Stroke::new((6.0 * zoom).clamp(2.0, 8.0), pal.border_default),
                );
                painter.line_segment([org_l, org_r], Stroke::new(1.2_f32, pal.grid_line));
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
                    Pos2::new(pm_l.x.max(rect.left() + 12.0), pm_l.y - 6.0 * zoom - 2.0),
                    Align2::LEFT_BOTTOM,
                    SubcellularLocation::PlasmaMembrane.label(),
                    FontId::proportional((label_sz * 0.88).max(8.5)),
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

                // 3. Draw cubic Bezier curved causal edges deflected around intermediate nodes
                let base_half = Vec2::new((32.0 * zoom).max(16.0), (16.0 * zoom).max(9.0));
                let obs_pad = (7.0 * zoom).clamp(4.0, 10.0);
                let mut edge_obstacles: Vec<EdgeObstacle> = Vec::with_capacity(world_positions.len());
                for (nid, &wpos) in &world_positions {
                    if let Some(ref ego) = ego_visible {
                        if !ego.contains(nid) {
                            continue;
                        }
                    }
                    let kind = self
                        .kg
                        .nodes
                        .get(nid)
                        .map(|n| n.kind)
                        .unwrap_or(EntityKind::GeneOther);
                    let scale = match kind {
                        EntityKind::TransmembraneReceptor | EntityKind::IonChannel => {
                            Vec2::new(1.0, 1.42)
                        }
                        EntityKind::Kinase | EntityKind::Phosphatase | EntityKind::Enzyme => {
                            Vec2::new(1.05, 1.28)
                        }
                        _ => Vec2::new(1.0, 1.05),
                    };
                    edge_obstacles.push(EdgeObstacle {
                        id: nid.clone(),
                        center: to_screen(wpos.x, wpos.y),
                        half_extents: Vec2::new(
                            base_half.x * scale.x + obs_pad,
                            base_half.y * scale.y + obs_pad,
                        ),
                    });
                }

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

                        draw_deflected_bezier_causal_edge(
                            &painter,
                            p1,
                            p2,
                            (20.0 * zoom).clamp(8.0, 36.0),
                            (24.0 * zoom).clamp(10.0, 42.0),
                            bow,
                            edge.sign,
                            Stroke::new((1.65 * zoom).clamp(0.9, 2.2), edge_col),
                            &edge_obstacles,
                            &edge.from,
                            &edge.to,
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
                            Stroke::new(2.4_f32, Color32::from_rgb(250, 204, 21))
                        } else if st.observed_log2fc.is_some() {
                            Stroke::new(1.8_f32, pal.text_primary)
                        } else {
                            Stroke::new(1.2_f32, pal.border_default)
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

                // 5. Bottom-left Edge-Sign Legend
                draw_pathway_edge_legend(&painter, rect, pal, self.palette_mode);

                // 6. Toggleable Bottom-right Mini-Map with header hide button & live viewport box
                if self.pathway_minimap_visible {
                    painter.rect_filled(
                        minimap_rect,
                        5.0,
                        Color32::from_rgba_unmultiplied(12, 16, 24, 224),
                    );
                    let mm_header_rect = Rect::from_min_max(
                        minimap_rect.min,
                        Pos2::new(minimap_rect.right(), minimap_rect.top() + minimap_header_h),
                    );
                    painter.rect_filled(
                        mm_header_rect,
                        egui::Rounding {
                            nw: 5.0,
                            ne: 5.0,
                            sw: 0.0,
                            se: 0.0,
                        },
                        Color32::from_rgba_unmultiplied(24, 31, 44, 238),
                    );
                    painter.line_segment(
                        [
                            Pos2::new(minimap_rect.left(), minimap_rect.top() + minimap_header_h),
                            Pos2::new(minimap_rect.right(), minimap_rect.top() + minimap_header_h),
                        ],
                        Stroke::new(1.0_f32, pal.grid_line),
                    );
                    painter.rect_stroke(minimap_rect, 5.0, Stroke::new(1.2_f32, pal.border_default));

                    let map_to_mini = |wx: f32, wy: f32| -> Pos2 {
                        Pos2::new(
                            minimap_inner.left()
                                + (wx / wl.world_w.max(1.0)).clamp(0.0, 1.0)
                                    * minimap_inner.width(),
                            minimap_inner.top()
                                + (wy / wl.world_h.max(1.0)).clamp(0.0, 1.0)
                                    * minimap_inner.height(),
                        )
                    };

                    let mm_pm_y = map_to_mini(0.0, wl.bands.plasma_membrane_y).y;
                    let mm_nm_y = map_to_mini(0.0, wl.bands.nuclear_envelope_y).y;
                    painter.line_segment(
                        [
                            Pos2::new(minimap_inner.left(), mm_pm_y),
                            Pos2::new(minimap_inner.right(), mm_pm_y),
                        ],
                        Stroke::new(0.8_f32, pal.grid_line),
                    );
                    painter.line_segment(
                        [
                            Pos2::new(minimap_inner.left(), mm_nm_y),
                            Pos2::new(minimap_inner.right(), mm_nm_y),
                        ],
                        Stroke::new(0.8_f32, pal.grid_line),
                    );

                    for edge in &wl.edges {
                        if let (Some(&wp1), Some(&wp2)) =
                            (world_positions.get(&edge.from), world_positions.get(&edge.to))
                        {
                            painter.line_segment(
                                [map_to_mini(wp1.x, wp1.y), map_to_mini(wp2.x, wp2.y)],
                                Stroke::new(
                                    0.7_f32,
                                    Color32::from_rgba_unmultiplied(120, 136, 160, 110),
                                ),
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
                        Stroke::new(1.4_f32, Color32::from_rgb(56, 189, 248)),
                    );
                    painter.text(
                        Pos2::new(
                            minimap_rect.left() + 7.0,
                            minimap_rect.top() + minimap_header_h * 0.5,
                        ),
                        Align2::LEFT_CENTER,
                        "MINI-MAP",
                        FontId::proportional(8.5),
                        pal.text_muted,
                    );
                    painter.text(
                        Pos2::new(
                            minimap_rect.right() - 22.0,
                            minimap_rect.top() + minimap_header_h * 0.5,
                        ),
                        Align2::RIGHT_CENTER,
                        format!("{:.0}%", zoom * 100.0),
                        FontId::proportional(8.5),
                        pal.text_muted,
                    );
                    painter.text(
                        minimap_close_rect.center(),
                        Align2::CENTER_CENTER,
                        "×",
                        FontId::proportional(11.0),
                        pal.text_primary,
                    );
                } else {
                    painter.rect_filled(
                        minimap_restore_rect,
                        4.0,
                        Color32::from_rgba_unmultiplied(24, 31, 44, 215),
                    );
                    painter.rect_stroke(
                        minimap_restore_rect,
                        4.0,
                        Stroke::new(1.0_f32, pal.border_default),
                    );
                    painter.text(
                        minimap_restore_rect.center(),
                        Align2::CENTER_CENTER,
                        "Show Mini-Map",
                        FontId::proportional(10.0),
                        Color32::from_rgb(232, 237, 244),
                    );
                }
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
        let ura_rows = self.filtered_sorted_ura(&rep);

        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(
                Vec2::new(560.0, total_avail_h),
                egui::Layout::top_down(egui::Align::Min),
                |left_ui| {
                    left_ui.horizontal_wrapped(|ui| {
                        ui.heading("Upstream Regulators (URA)");
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.button("Export CSV").clicked() {
                                actions.push(GuiAction::ExportUraCsvDialog);
                            }
                        });
                    });
                    left_ui.horizontal_wrapped(|ui| {
                        ui.add(
                            egui::TextEdit::singleline(&mut self.regulator_search)
                                .hint_text("Filter regulator or type...")
                                .desired_width(155.0),
                        );
                        if !self.regulator_search.is_empty() && ui.small_button("✕").clicked() {
                            self.regulator_search.clear();
                        }
                        for (st_f, lbl) in [
                            (RegStateFilter::All, "All"),
                            (RegStateFilter::Activated, "Activated"),
                            (RegStateFilter::Inhibited, "Inhibited"),
                        ] {
                            if ui
                                .selectable_label(self.ura_state_filter == st_f, lbl)
                                .clicked()
                            {
                                self.ura_state_filter = st_f;
                            }
                        }
                    });

                    left_ui.push_id("ura_table_scroll", |ui| {
                        egui::ScrollArea::vertical()
                            .max_height((total_avail_h * 0.44).max(210.0))
                            .show(ui, |ui| {
                                egui::Grid::new("ura_grid")
                                    .striped(true)
                                    .num_columns(7)
                                    .spacing([10.0, 5.0])
                                    .show(ui, |ui| {
                                        render_sort_header(
                                            ui,
                                            "Regulator",
                                            UraSortCol::Regulator,
                                            SortDir::Asc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "Type",
                                            UraSortCol::Type,
                                            SortDir::Asc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "State",
                                            UraSortCol::State,
                                            SortDir::Asc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "Act Z",
                                            UraSortCol::ActZ,
                                            SortDir::Desc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "Bias Z",
                                            UraSortCol::BiasZ,
                                            SortDir::Desc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "Overlap P",
                                            UraSortCol::OverlapP,
                                            SortDir::Asc,
                                            &mut self.ura_sort,
                                        );
                                        render_sort_header(
                                            ui,
                                            "Targets",
                                            UraSortCol::Targets,
                                            SortDir::Desc,
                                            &mut self.ura_sort,
                                        );
                                        ui.end_row();

                                        for (orig_idx, h) in &ura_rows {
                                            let idx = *orig_idx;
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
                                                format!("{}{skew_mark}", ca_core::format_z_2dp(h.activation_z)),
                                            );
                                            ui.label(ca_core::format_z_2dp(h.bias_corrected_z));
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
                        "Force-directed layered causal cascade connecting upstream co-regulators & targets (scroll to zoom, drag background to pan, drag nodes to pin).",
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

                let shown_targets: Vec<TargetGeneRef<'_>> =
                    select_stratified_mechanistic_targets(&sel_u.target_evidences, 14)
                        .into_iter()
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
                painter.rect_stroke(rect, 6.0, Stroke::new(1.0_f32, pal.grid_line));

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
                            let base_col = match tg.sign {
                                EdgeSign::Activating => pal.activated_strong,
                                EdgeSign::Inhibiting => pal.inhibited_strong,
                                EdgeSign::Ambiguous => pal.border_default,
                            };
                            let edge_col = lerp_color(pal.bg_canvas, base_col, 0.84);
                            draw_bezier_causal_edge(
                                &painter,
                                rp,
                                tp,
                                (18.0 * zoom).clamp(8.0, 28.0),
                                (20.0 * zoom).clamp(10.0, 30.0),
                                bow,
                                tg.sign,
                                Stroke::new((1.55 * zoom).clamp(1.0, 2.1), edge_col),
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
                            Stroke::new(2.0_f32, pal.text_primary),
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
                            Stroke::new(1.3_f32, pal.border_default),
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
        &mut self,
        ui: &mut Ui,
        pal: &ScientificPalette,
        actions: &mut Vec<GuiAction>,
    ) {
        let rep = self.reports[self.active_cohort_idx].clone();
        let cna_rows = self.filtered_sorted_cna(&rep);
        ui.columns(2, |cols| {
            cols[0].horizontal_wrapped(|ui| {
                ui.heading("Multi-Hop Causal Networks (CNA, Depth K=1..3)");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Export CSV").clicked() {
                        actions.push(GuiAction::ExportCnaCsvDialog);
                    }
                });
            });
            cols[0].label(
                RichText::new(
                    "Identifies root master regulators (even when transcriptionally unchanged) via multi-step signaling cascades and degree-binned Monte Carlo p-values.",
                )
                .size(11.5)
                .color(pal.text_muted),
            );
            cols[0].horizontal_wrapped(|ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut self.cna_search)
                        .hint_text("Filter root or relay...")
                        .desired_width(170.0),
                );
                if !self.cna_search.is_empty() && ui.small_button("✕").clicked() {
                    self.cna_search.clear();
                }
                for (st_f, lbl) in [
                    (RegStateFilter::All, "All"),
                    (RegStateFilter::Activated, "Activated"),
                    (RegStateFilter::Inhibited, "Inhibited"),
                ] {
                    if ui
                        .selectable_label(self.cna_state_filter == st_f, lbl)
                        .clicked()
                    {
                        self.cna_state_filter = st_f;
                    }
                }
            });

            cols[0].push_id("cna_scroll", |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("cna_grid")
                        .striped(true)
                        .num_columns(6)
                        .spacing([10.0, 5.0])
                        .show(ui, |ui| {
                            render_sort_header(
                                ui,
                                "Master Root",
                                CnaSortCol::Root,
                                SortDir::Asc,
                                &mut self.cna_sort,
                            );
                            render_sort_header(
                                ui,
                                "Depth",
                                CnaSortCol::Depth,
                                SortDir::Asc,
                                &mut self.cna_sort,
                            );
                            render_sort_header(
                                ui,
                                "State",
                                CnaSortCol::State,
                                SortDir::Asc,
                                &mut self.cna_sort,
                            );
                            render_sort_header(
                                ui,
                                "Act Z",
                                CnaSortCol::ActZ,
                                SortDir::Desc,
                                &mut self.cna_sort,
                            );
                            render_sort_header(
                                ui,
                                "Overlap P",
                                CnaSortCol::OverlapP,
                                SortDir::Asc,
                                &mut self.cna_sort,
                            );
                            render_sort_header(
                                ui,
                                "Net-Bias P",
                                CnaSortCol::NetBiasP,
                                SortDir::Asc,
                                &mut self.cna_sort,
                            );
                            ui.end_row();

                            for (orig_idx, h) in &cna_rows {
                                let idx = *orig_idx;
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
                                ui.colored_label(c, ca_core::format_z_2dp(h.activation_z));
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
                    Stroke::new(if is_sel { 2.4 } else { 1.0_f32 }, pal.text_primary),
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
                            format!("Z = {} ({} genes)", ca_core::format_z_2dp(r.activation_z), r.overlap_count),
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
                            "{}\nCategory: {}\nActivation Z: {}\n-log10(P): {:.2}",
                            r.term_name,
                            r.category,
                            ca_core::format_z_2dp(r.activation_z),
                            r.neg_log10_p
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
                                ui.colored_label(c, ca_core::format_z_2dp(d.activation_z));
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
            painter.rect_stroke(rect, 6.0, Stroke::new(1.0_f32, pal.grid_line));

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
                        Stroke::new(2.0_f32, pal.text_primary),
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
                        Stroke::new(1.5_f32, pal.text_primary),
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
                        Stroke::new(2.0_f32, pal.text_primary),
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
        let total_avail_w = ui.available_width().max(760.0);
        let total_avail_h = ui.available_height().max(380.0);
        let gap = 12.0;
        let left_w = (total_avail_w * 0.40)
            .clamp(470.0, 525.0)
            .min((total_avail_w - 420.0).max(320.0));
        let right_w = (total_avail_w - left_w - gap - 8.0).max(380.0);

        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(
                Vec2::new(left_w, total_avail_h),
                egui::Layout::top_down(egui::Align::Min),
                |left_ui| {
                    left_ui.set_max_width(left_w);
                    left_ui.heading("Multi-Cohort Comparison Heatmap");
                    left_ui.horizontal_wrapped(|ui| {
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
                    if let Some((rep_a, rep_b)) = self.diff_report_pair() {
                        let digest = compute_regulator_rank_diff(rep_a, rep_b);
                        left_ui.add_space(4.0);
                        left_ui.horizontal_wrapped(|ui| {
                            ui.label(RichText::new("Rank Diff:").strong().size(11.5));
                            egui::ComboBox::from_id_salt("diff_cohort_a")
                                .selected_text(rep_a.dataset_name.chars().take(16).collect::<String>())
                                .width(120.0)
                                .show_ui(ui, |ui| {
                                    for (idx, r) in self.reports.iter().enumerate() {
                                        if ui
                                            .selectable_label(self.diff_cohort_a_idx == idx, &r.dataset_name)
                                            .clicked()
                                        {
                                            actions.push(GuiAction::SetDiffCohortA(idx));
                                        }
                                    }
                                });
                            ui.label("vs");
                            egui::ComboBox::from_id_salt("diff_cohort_b")
                                .selected_text(rep_b.dataset_name.chars().take(16).collect::<String>())
                                .width(120.0)
                                .show_ui(ui, |ui| {
                                    for (idx, r) in self.reports.iter().enumerate() {
                                        if ui
                                            .selectable_label(self.diff_cohort_b_idx == idx, &r.dataset_name)
                                            .clicked()
                                        {
                                            actions.push(GuiAction::SetDiffCohortB(idx));
                                        }
                                    }
                                });
                            if ui.small_button("Diff JSON").clicked() {
                                actions.push(GuiAction::ExportRankDiffJsonDialog);
                            }
                            if ui.small_button("Diff CSV").clicked() {
                                actions.push(GuiAction::ExportRankDiffCsvDialog);
                            }
                        });
                        left_ui.label(
                            RichText::new(format!(
                                "Regulator Rank Diff: {} shared | {} flips | Spearman rho = {:.3}",
                                digest.shared_count,
                                digest.direction_flip_count,
                                digest.spearman_rank_rho
                            ))
                            .size(11.0)
                            .color(pal.text_muted),
                        );
                        left_ui.add_space(4.0);
                    }
                    let entity_max_w = (left_w - (comp.cohort_names.len() as f32 * 72.0) - 28.0).max(150.0);
                    left_ui.push_id("comp_heatmap_scroll", |ui| {
                        egui::ScrollArea::both().show(ui, |ui| {
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
                                        ui.vertical(|ui| {
                                            ui.set_max_width(entity_max_w);
                                            ui.add(egui::Label::new(label).wrap());
                                        });
                                        for &z in &comp.values[r_idx] {
                                            let (cell_r, cell_p) = ui.allocate_painter(
                                                Vec2::new(66.0, 20.0),
                                                Sense::hover(),
                                            );
                                            cell_p.rect_filled(cell_r.rect, 3.0, pal.color_for_z(z));
                                            cell_p.text(
                                                cell_r.rect.center(),
                                                Align2::CENTER_CENTER,
                                                ca_core::format_z_2dp(z),
                                                FontId::proportional(10.5),
                                                pal.text_primary,
                                            );
                                        }
                                        ui.end_row();
                                    }
                                });
                        });
                    });
                },
            );

            ui.add_space(gap);

            ui.allocate_ui_with_layout(
                Vec2::new(right_w, total_avail_h),
                egui::Layout::top_down(egui::Align::Min),
                |right_ui| {
                    right_ui.set_max_width(right_w);
                    let rep = &self.reports[self.active_cohort_idx];
                    let active = self.cohorts[self.active_cohort_idx].filter_active(
                        self.config.abs_log2fc_cutoff,
                        self.config.fdr_cutoff,
                        &self.kg,
                    );
                    let deg_fc: HashMap<&str, f64> = active
                        .degs
                        .iter()
                        .map(|(g, m)| (g.as_str(), m.log2_fc))
                        .collect();

                    right_ui.heading("Perturbation Connectivity (SigCom LINCS Open Signatures)");
                    right_ui.add(
                        egui::Label::new(
                            RichText::new(
                                "Negative scores (blue) identify candidate Reversers opposing the disease/stimulus state; positive scores (orange) identify Mimickers. Click any row to inspect matched genes or Rescue in Pathway Atlas.",
                            )
                            .size(11.5)
                            .color(pal.text_muted),
                        )
                        .wrap(),
                    );

                    // Allocate explicit responsive budgets across all 6 columns so none clip at 1600x1000.
                    let grid_avail = (right_w - 28.0 - 30.0).max(520.0);
                    let scale = (grid_avail / 720.0).clamp(0.75, 1.2);
                    let w_pert = 135.0 * scale;
                    let w_score = 44.0 * scale;
                    let w_role = 115.0 * scale;
                    let w_genes = 200.0 * scale;
                    let w_moa = 170.0 * scale;
                    let w_action = 56.0 * scale;

                    right_ui.push_id("perturb_scroll", |ui| {
                        egui::ScrollArea::both()
                            .max_height((total_avail_h * 0.58).max(220.0))
                            .show(ui, |ui| {
                                egui::Grid::new("perturb_grid")
                                    .striped(true)
                                    .num_columns(6)
                                    .spacing([6.0, 5.0])
                                    .show(ui, |ui| {
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_pert);
                                            ui.strong("Perturbagen");
                                        });
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_score);
                                            ui.strong("Score");
                                        });
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_role);
                                            ui.strong("Role");
                                        });
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_genes);
                                            ui.strong("Reversed / Matched Genes");
                                        });
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_moa);
                                            ui.strong("Mechanism of Action");
                                        });
                                        ui.vertical(|ui| {
                                            ui.set_max_width(w_action);
                                            ui.strong("Action");
                                        });
                                        ui.end_row();

                                        for (idx, p) in rep.perturbation_matches.iter().enumerate() {
                                            let is_sel = self.selected_perturb_idx == idx;
                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_pert);
                                                ui.style_mut().wrap_mode =
                                                    Some(egui::TextWrapMode::Wrap);
                                                if ui
                                                    .selectable_label(is_sel, &p.perturbagen)
                                                    .clicked()
                                                {
                                                    actions.push(GuiAction::SelectPerturbIndex(idx));
                                                }
                                            });
                                            let c = pal.color_for_z(p.connectivity_score / 25.0);
                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_score);
                                                ui.colored_label(
                                                    c,
                                                    format!("{:+.1}", p.connectivity_score),
                                                );
                                            });
                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_role);
                                                ui.add(
                                                    egui::Label::new(
                                                        RichText::new(format!(
                                                            "{} ({})",
                                                            if p.connectivity_score <= -25.0 {
                                                                "Reverser"
                                                            } else if p.connectivity_score >= 25.0 {
                                                                "Mimicker"
                                                            } else {
                                                                "Orthogonal"
                                                            },
                                                            p.kind.label()
                                                        ))
                                                        .color(c),
                                                    )
                                                    .wrap(),
                                                );
                                            });

                                            let is_rev = p.connectivity_score <= -25.0;
                                            let rev_up: Vec<&str> = p
                                                .matched_up
                                                .iter()
                                                .filter(|g| {
                                                    deg_fc.get(g.as_str()).copied().unwrap_or(0.0)
                                                        < 0.0
                                                })
                                                .map(String::as_str)
                                                .collect();
                                            let rev_down: Vec<&str> = p
                                                .matched_down
                                                .iter()
                                                .filter(|g| {
                                                    deg_fc.get(g.as_str()).copied().unwrap_or(0.0)
                                                        > 0.0
                                                })
                                                .map(String::as_str)
                                                .collect();

                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_genes);
                                                let mut shown_any = false;
                                                if !p.matched_up.is_empty() {
                                                    shown_any = true;
                                                    let preview: Vec<&str> =
                                                        if is_rev && !rev_up.is_empty() {
                                                            rev_up.iter().take(4).copied().collect()
                                                        } else {
                                                            p.matched_up
                                                                .iter()
                                                                .take(4)
                                                                .map(String::as_str)
                                                                .collect()
                                                        };
                                                    let lbl = if is_rev {
                                                        format!(
                                                            "Up(+{}): {}",
                                                            p.matched_up.len(),
                                                            preview.join(", ")
                                                        )
                                                    } else {
                                                        format!(
                                                            "Up({}): {}",
                                                            p.matched_up.len(),
                                                            preview.join(", ")
                                                        )
                                                    };
                                                    ui.add(
                                                        egui::Label::new(
                                                            RichText::new(lbl)
                                                                .size(11.0)
                                                                .color(pal.expr_up),
                                                        )
                                                        .wrap(),
                                                    );
                                                }
                                                if !p.matched_down.is_empty() {
                                                    shown_any = true;
                                                    let preview: Vec<&str> =
                                                        if is_rev && !rev_down.is_empty() {
                                                            rev_down
                                                                .iter()
                                                                .take(4)
                                                                .copied()
                                                                .collect()
                                                        } else {
                                                            p.matched_down
                                                                .iter()
                                                                .take(4)
                                                                .map(String::as_str)
                                                                .collect()
                                                        };
                                                    let lbl = if is_rev {
                                                        format!(
                                                            "Down(-{}): {}",
                                                            p.matched_down.len(),
                                                            preview.join(", ")
                                                        )
                                                    } else {
                                                        format!(
                                                            "Down({}): {}",
                                                            p.matched_down.len(),
                                                            preview.join(", ")
                                                        )
                                                    };
                                                    ui.add(
                                                        egui::Label::new(
                                                            RichText::new(lbl)
                                                                .size(11.0)
                                                                .color(pal.inhibited_strong),
                                                        )
                                                        .wrap(),
                                                    );
                                                }
                                                if !shown_any {
                                                    ui.label(
                                                        RichText::new("-").color(pal.text_muted),
                                                    );
                                                }
                                            });

                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_moa);
                                                ui.add(
                                                    egui::Label::new(
                                                        RichText::new(&p.mechanism_of_action)
                                                            .size(11.0),
                                                    )
                                                    .wrap(),
                                                );
                                            });
                                            ui.vertical(|ui| {
                                                ui.set_max_width(w_action);
                                                if ui
                                                    .small_button("Rescue")
                                                    .on_hover_text(
                                                        "Clamp targets & reversed genes in Pathway Atlas",
                                                    )
                                                    .clicked()
                                                {
                                                    actions.push(GuiAction::ApplyDrugRescue {
                                                        hit_idx: idx,
                                                        switch_tab: true,
                                                    });
                                                }
                                            });
                                            ui.end_row();
                                        }
                                    });
                            });
                    });

                    if let Some(sel_p) = rep.perturbation_matches.get(self.selected_perturb_idx) {
                        right_ui.add_space(6.0);
                        right_ui.group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(
                                    RichText::new(format!(
                                        "Selected Perturbagen: {} ({:+.1} — {})",
                                        sel_p.perturbagen, sel_p.connectivity_score, sel_p.mode
                                    ))
                                    .strong(),
                                );
                                if ui.small_button("Clamp Drug Rescue in Pathway Atlas").clicked() {
                                    actions.push(GuiAction::ApplyDrugRescue {
                                        hit_idx: self.selected_perturb_idx,
                                        switch_tab: true,
                                    });
                                }
                            });
                            ui.label(
                                RichText::new(format!(
                                    "Context: {} | Source: {}",
                                    sel_p.cell_context, sel_p.source_db
                                ))
                                .size(11.0)
                                .color(pal.text_muted),
                            );
                            if !sel_p.matched_up.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(
                                        RichText::new(format!(
                                            "Matched Up ({}):",
                                            sel_p.matched_up.len()
                                        ))
                                        .size(11.5)
                                        .strong()
                                        .color(pal.expr_up),
                                    );
                                    for g in &sel_p.matched_up {
                                        let fc = deg_fc.get(g.as_str()).copied().unwrap_or(0.0);
                                        let rev_tag = if fc < 0.0 { " [REV]" } else { "" };
                                        if ui
                                            .small_button(format!("{g} ({fc:+.2}{rev_tag})"))
                                            .on_hover_text("Jump to gene in Pathway Atlas")
                                            .clicked()
                                        {
                                            actions.push(GuiAction::JumpToGene {
                                                symbol: g.clone(),
                                                pathway_id: None,
                                            });
                                        }
                                    }
                                });
                            }
                            if !sel_p.matched_down.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(
                                        RichText::new(format!(
                                            "Matched Down ({}):",
                                            sel_p.matched_down.len()
                                        ))
                                        .size(11.5)
                                        .strong()
                                        .color(pal.inhibited_strong),
                                    );
                                    for g in &sel_p.matched_down {
                                        let fc = deg_fc.get(g.as_str()).copied().unwrap_or(0.0);
                                        let rev_tag = if fc > 0.0 { " [REV]" } else { "" };
                                        if ui
                                            .small_button(format!("{g} ({fc:+.2}{rev_tag})"))
                                            .on_hover_text("Jump to gene in Pathway Atlas")
                                            .clicked()
                                        {
                                            actions.push(GuiAction::JumpToGene {
                                                symbol: g.clone(),
                                                pathway_id: None,
                                            });
                                        }
                                    }
                                });
                            }
                        });
                    }
                },
            );
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
                    if ui.button("Import .json.gz").clicked() {
                        actions.push(GuiAction::LoadSnapshotDialog);
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

    fn process_keyboard_shortcuts(
        &mut self,
        ctx: &egui::Context,
        wants_kb: bool,
        actions: &mut Vec<GuiAction>,
    ) {
        // Always active: Ctrl/Cmd+K (Command Palette), Ctrl/Cmd+F (Omnibox), F1 (About), layered Escape
        if ctx.input(|i| (i.modifiers.ctrl || i.modifiers.command) && i.key_pressed(egui::Key::K)) {
            self.show_command_palette = !self.show_command_palette;
            if self.show_command_palette {
                self.command_palette_query.clear();
                self.command_palette_selected = 0;
                self.focus_command_palette = true;
            }
            return;
        }
        if self.show_command_palette {
            let (esc, up, down, enter) = ctx.input(|i| {
                (
                    i.key_pressed(egui::Key::Escape),
                    i.key_pressed(egui::Key::ArrowUp),
                    i.key_pressed(egui::Key::ArrowDown),
                    i.key_pressed(egui::Key::Enter),
                )
            });
            let matches = self.filtered_palette_commands();
            if esc {
                self.show_command_palette = false;
            } else if !matches.is_empty() {
                if self.command_palette_selected >= matches.len() {
                    self.command_palette_selected = 0;
                }
                if up {
                    self.command_palette_selected = self
                        .command_palette_selected
                        .checked_sub(1)
                        .unwrap_or(matches.len() - 1);
                } else if down {
                    self.command_palette_selected =
                        (self.command_palette_selected + 1) % matches.len();
                } else if enter {
                    let chosen = matches[self.command_palette_selected].action;
                    self.show_command_palette = false;
                    self.execute_palette_command(chosen, actions);
                }
            }
            return;
        }
        if ctx.input(|i| (i.modifiers.ctrl || i.modifiers.command) && i.key_pressed(egui::Key::F)) {
            self.omnibox_focus_requested = true;
            self.omnibox_dropdown_open = true;
        }
        if ctx.input(|i| i.key_pressed(egui::Key::F1)) {
            self.about_open = !self.about_open;
        }
        if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            if self.omnibox_dropdown_open {
                self.omnibox_dropdown_open = false;
            } else if self.about_open {
                self.about_open = false;
            } else if self.pathway_fullscreen {
                self.pathway_fullscreen = false;
                self.pathway_fit_pending = true;
            }
        }

        if wants_kb {
            return;
        }

        let (cmd_o, cmd_s, cmd_e, plain_key) = ctx.input(|i| {
            let cmd = i.modifiers.ctrl || i.modifiers.command;
            let no_mods = !cmd && !i.modifiers.alt;
            let mut plain = None;
            if no_mods {
                for (k, tag) in [
                    (egui::Key::Num1, 1u8),
                    (egui::Key::Num2, 2),
                    (egui::Key::Num3, 3),
                    (egui::Key::Num4, 4),
                    (egui::Key::Num5, 5),
                    (egui::Key::F, 6),
                    (egui::Key::Num0, 6),
                    (egui::Key::R, 7),
                ] {
                    if i.key_pressed(k) {
                        plain = Some(tag);
                        break;
                    }
                }
            }
            (
                cmd && i.key_pressed(egui::Key::O),
                cmd && i.key_pressed(egui::Key::S),
                cmd && i.key_pressed(egui::Key::E),
                plain,
            )
        });

        if cmd_o {
            actions.push(GuiAction::ImportDatasetDialog);
        }
        if cmd_s {
            actions.push(GuiAction::ExportJsonReportDialog);
        }
        if cmd_e {
            if self.active_tab == WorkbenchTab::PhenotypeDrugMatch
                && self.phenotype_subtab == PhenotypeDrugSubTab::DownstreamTreeMap
            {
                actions.push(GuiAction::ExportTreeMapSvgDialog);
            } else {
                actions.push(GuiAction::ExportPathwaySvgDialog);
            }
        }
        match plain_key {
            Some(1) => self.active_tab = WorkbenchTab::CohortVolcano,
            Some(2) => self.active_tab = WorkbenchTab::PathwayAtlas,
            Some(3) => self.active_tab = WorkbenchTab::CausalRegulators,
            Some(4) => self.active_tab = WorkbenchTab::PhenotypeDrugMatch,
            Some(5) => self.active_tab = WorkbenchTab::KnowledgeDatabases,
            Some(6) => actions.push(GuiAction::FitPathwayViewport),
            Some(7) => actions.push(GuiAction::ClearAllClamps),
            _ => {}
        }
    }

    fn diff_report_pair(&self) -> Option<(&CoreAnalysisReport, &CoreAnalysisReport)> {
        let a = self.reports.get(self.diff_cohort_a_idx)?;
        let b = self
            .reports
            .get(self.diff_cohort_b_idx)
            .or_else(|| self.reports.first())?;
        Some((a, b))
    }

    fn filtered_palette_commands(&self) -> Vec<&'static CaPaletteCommand> {
        let mut scored: Vec<(i32, usize, &'static CaPaletteCommand)> = CA_PALETTE_COMMANDS
            .iter()
            .enumerate()
            .filter_map(|(idx, cmd)| {
                ca_fuzzy_palette_score(&self.command_palette_query, cmd.label, cmd.keywords)
                    .map(|sc| (sc, idx, cmd))
            })
            .collect();
        scored.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
        scored.into_iter().map(|(_, _, cmd)| cmd).collect()
    }

    fn execute_palette_command(
        &mut self,
        action: CaPaletteAction,
        actions: &mut Vec<GuiAction>,
    ) {
        match action {
            CaPaletteAction::TabCohortVolcano => self.active_tab = WorkbenchTab::CohortVolcano,
            CaPaletteAction::TabPathwayAtlas => self.active_tab = WorkbenchTab::PathwayAtlas,
            CaPaletteAction::TabCausalRegulators => {
                self.active_tab = WorkbenchTab::CausalRegulators;
            }
            CaPaletteAction::TabPhenotypeComparison => {
                self.active_tab = WorkbenchTab::PhenotypeDrugMatch;
                self.phenotype_subtab = PhenotypeDrugSubTab::ComparisonAndReversers;
            }
            CaPaletteAction::TabKnowledgeDatabases => {
                self.active_tab = WorkbenchTab::KnowledgeDatabases;
            }
            CaPaletteAction::FocusOmnibox => {
                self.omnibox_focus_requested = true;
                self.omnibox_dropdown_open = true;
            }
            CaPaletteAction::RunCoreAnalysis => actions.push(GuiAction::RecomputeActive),
            CaPaletteAction::SyncAllCohorts => actions.push(GuiAction::RecomputeAllComparison),
            CaPaletteAction::FitPathwayCamera => actions.push(GuiAction::FitPathwayViewport),
            CaPaletteAction::ClearMapClamps => actions.push(GuiAction::ClearAllClamps),
            CaPaletteAction::ToggleTheme => actions.push(GuiAction::TogglePalette),
            CaPaletteAction::ExportReportJson => actions.push(GuiAction::ExportJsonReportDialog),
            CaPaletteAction::ExportPathwaySvg => actions.push(GuiAction::ExportPathwaySvgDialog),
            CaPaletteAction::ExportRankDiffJson => {
                actions.push(GuiAction::ExportRankDiffJsonDialog);
            }
            CaPaletteAction::ExportRankDiffCsv => {
                actions.push(GuiAction::ExportRankDiffCsvDialog);
            }
            CaPaletteAction::ToggleAbout => actions.push(GuiAction::ToggleAbout),
        }
    }

    fn render_command_palette_modal(
        &mut self,
        ctx: &egui::Context,
        actions: &mut Vec<GuiAction>,
    ) {
        let mut chosen_action: Option<CaPaletteAction> = None;
        let matches = self.filtered_palette_commands();
        if !matches.is_empty() && self.command_palette_selected >= matches.len() {
            self.command_palette_selected = 0;
        }
        egui::Window::new("Command Palette (Ctrl+K)")
            .collapsible(false)
            .resizable(false)
            .fixed_size(Vec2::new(540.0, 360.0))
            .anchor(Align2::CENTER_TOP, Vec2::new(0.0, 68.0))
            .show(ctx, |ui| {
                let edit = ui.add(
                    egui::TextEdit::singleline(&mut self.command_palette_query)
                        .hint_text("Type a tab, regulator diff export, MAP action, or shortcut...")
                        .desired_width(f32::INFINITY),
                );
                if self.focus_command_palette {
                    edit.request_focus();
                    self.focus_command_palette = false;
                }
                if edit.changed() {
                    self.command_palette_selected = 0;
                }
                ui.label(
                    RichText::new("Up/Down: navigate • Enter: run • Esc: close")
                        .monospace()
                        .size(10.5),
                );
                ui.separator();
                egui::ScrollArea::vertical()
                    .max_height(265.0)
                    .show(ui, |ui| {
                        if matches.is_empty() {
                            ui.label("No matching commands.");
                        } else {
                            for (idx, cmd) in matches.iter().enumerate() {
                                let selected = idx == self.command_palette_selected;
                                let row_text = if cmd.hotkey.is_empty() {
                                    format!("[{:<10}] {}", cmd.category, cmd.label)
                                } else {
                                    format!("[{:<10}] {}   ({})", cmd.category, cmd.label, cmd.hotkey)
                                };
                                let resp = ui.selectable_label(
                                    selected,
                                    RichText::new(row_text).monospace().size(11.5),
                                );
                                if resp.hovered() {
                                    self.command_palette_selected = idx;
                                }
                                if resp.clicked() {
                                    chosen_action = Some(cmd.action);
                                }
                            }
                        }
                    });
            });
        if let Some(act) = chosen_action {
            self.show_command_palette = false;
            self.execute_palette_command(act, actions);
        }
    }

    fn save_csv_dialog(&mut self, default_name: &str, csv_text: &str, label: &str) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("Comma-Separated Values (CSV)", &["csv"])
            .set_file_name(default_name)
            .save_file()
        {
            match std::fs::write(&path, csv_text) {
                Ok(()) => {
                    self.status_banner = format!("Exported {label} CSV to {}", path.display());
                }
                Err(e) => {
                    self.status_banner = format!("CSV export error: {e}");
                }
            }
        }
    }

    fn filtered_sorted_genes(&self) -> Vec<&ca_core::GeneMeasurement> {
        let ds = &self.cohorts[self.active_cohort_idx];
        let fc_cut = self.config.abs_log2fc_cutoff;
        let q_cut = self.config.fdr_cutoff;
        let needle = self.gene_search.trim().to_uppercase();

        let mut rows: Vec<&ca_core::GeneMeasurement> = ds
            .rows
            .iter()
            .filter(|r| {
                if !needle.is_empty() && !r.symbol.to_uppercase().contains(&needle) {
                    return false;
                }
                let pass = r.log2_fc.abs() >= fc_cut && r.fdr <= q_cut;
                match self.gene_dir_filter {
                    GeneDirectionFilter::All => true,
                    GeneDirectionFilter::DegsOnly => pass,
                    GeneDirectionFilter::UpOnly => pass && r.log2_fc > 0.0,
                    GeneDirectionFilter::DownOnly => pass && r.log2_fc < 0.0,
                }
            })
            .collect();

        if let Some(sort) = self.gene_sort {
            rows.sort_by(|a, b| {
                let ord = match sort.col {
                    GeneSortCol::Symbol => a.symbol.cmp(&b.symbol),
                    GeneSortCol::Log2Fc => a.log2_fc.total_cmp(&b.log2_fc),
                    GeneSortCol::PValue => a.p_value.total_cmp(&b.p_value),
                    GeneSortCol::Fdr => a.fdr.total_cmp(&b.fdr),
                };
                sort.apply(ord).then_with(|| a.symbol.cmp(&b.symbol))
            });
        }
        rows
    }

    fn filtered_sorted_ura<'a>(
        &self,
        rep: &'a CoreAnalysisReport,
    ) -> Vec<(usize, &'a ca_core::UpstreamRegulatorHit)> {
        let needle = self.regulator_search.trim().to_uppercase();
        let mut rows: Vec<(usize, &'a ca_core::UpstreamRegulatorHit)> = rep
            .upstream_regulators
            .iter()
            .enumerate()
            .filter(|(_, h)| {
                self.ura_state_filter.matches(h.predicted_state)
                    && (needle.is_empty()
                        || h.regulator_id.to_uppercase().contains(&needle)
                        || h.regulator_name.to_uppercase().contains(&needle)
                        || h.kind.label().to_uppercase().contains(&needle))
            })
            .collect();

        if let Some(sort) = self.ura_sort {
            rows.sort_by(|(_, a), (_, b)| {
                let ord = match sort.col {
                    UraSortCol::Regulator => a.regulator_id.cmp(&b.regulator_id),
                    UraSortCol::Type => a.kind.label().cmp(b.kind.label()),
                    UraSortCol::State => a.predicted_state.label().cmp(b.predicted_state.label()),
                    UraSortCol::ActZ => a.activation_z.total_cmp(&b.activation_z),
                    UraSortCol::BiasZ => a.bias_corrected_z.total_cmp(&b.bias_corrected_z),
                    UraSortCol::OverlapP => a.p_value.total_cmp(&b.p_value),
                    UraSortCol::Targets => a.overlap_count.cmp(&b.overlap_count),
                };
                sort.apply(ord)
                    .then_with(|| a.regulator_id.cmp(&b.regulator_id))
            });
        }
        rows
    }

    fn filtered_sorted_cna<'a>(
        &self,
        rep: &'a CoreAnalysisReport,
    ) -> Vec<(usize, &'a ca_core::CausalNetworkHit)> {
        let needle = self.cna_search.trim().to_uppercase();
        let mut rows: Vec<(usize, &'a ca_core::CausalNetworkHit)> = rep
            .causal_networks
            .iter()
            .enumerate()
            .filter(|(_, h)| {
                self.cna_state_filter.matches(h.predicted_state)
                    && (needle.is_empty()
                        || h.root_regulator.to_uppercase().contains(&needle)
                        || h.root_name.to_uppercase().contains(&needle)
                        || h.kind.label().to_uppercase().contains(&needle)
                        || h.intermediate_regulators
                            .iter()
                            .any(|mid| mid.to_uppercase().contains(&needle)))
            })
            .collect();

        if let Some(sort) = self.cna_sort {
            rows.sort_by(|(_, a), (_, b)| {
                let ord = match sort.col {
                    CnaSortCol::Root => a.root_regulator.cmp(&b.root_regulator),
                    CnaSortCol::Depth => a.max_depth_k.cmp(&b.max_depth_k),
                    CnaSortCol::State => a.predicted_state.label().cmp(b.predicted_state.label()),
                    CnaSortCol::ActZ => a.activation_z.total_cmp(&b.activation_z),
                    CnaSortCol::OverlapP => a.overlap_p_value.total_cmp(&b.overlap_p_value),
                    CnaSortCol::NetBiasP => a.network_bias_p_value.total_cmp(&b.network_bias_p_value),
                };
                sort.apply(ord)
                    .then_with(|| a.root_regulator.cmp(&b.root_regulator))
            });
        }
        rows
    }

    fn export_genes_csv_string(&self) -> String {
        use std::fmt::Write as _;
        let fc_cut = self.config.abs_log2fc_cutoff;
        let q_cut = self.config.fdr_cutoff;
        let mut out = String::from("symbol,log2_fc,p_value,bh_fdr,passes_cutoffs\n");
        for r in self.filtered_sorted_genes() {
            let pass = r.log2_fc.abs() >= fc_cut && r.fdr <= q_cut;
            let _ = writeln!(
                out,
                "{},{:.6},{:.6e},{:.6e},{}",
                csv_escape(&r.symbol),
                r.log2_fc,
                r.p_value,
                r.fdr,
                pass
            );
        }
        out
    }

    fn export_ura_csv_string(&self) -> String {
        use std::fmt::Write as _;
        let rep = self.current_report();
        let mut out = String::from(
            "rank,regulator_id,regulator_name,kind,predicted_state,activation_z,bias_corrected_z,skew_flagged,overlap_count,regulon_size,overlap_p_value,bh_fdr,observed_log2fc,mechanistic_summary\n",
        );
        for (display_rank, (_, h)) in self.filtered_sorted_ura(rep).into_iter().enumerate() {
            let obs = h
                .observed_log2fc
                .map(|v| format!("{v:.4}"))
                .unwrap_or_else(|| "NA".to_string());
            let _ = writeln!(
                out,
                "{},{},{},{},{},{:.6},{:.6},{},{},{},{:.6e},{:.6e},{},{}",
                display_rank + 1,
                csv_escape(&h.regulator_id),
                csv_escape(&h.regulator_name),
                csv_escape(h.kind.label()),
                csv_escape(h.predicted_state.label()),
                h.activation_z,
                h.bias_corrected_z,
                h.skew_flagged,
                h.overlap_count,
                h.regulon_size,
                h.p_value,
                h.fdr,
                obs,
                csv_escape(h.mechanistic_summary.as_deref().unwrap_or("")),
            );
        }
        out
    }

    fn export_cna_csv_string(&self) -> String {
        use std::fmt::Write as _;
        let rep = self.current_report();
        let mut out = String::from(
            "rank,root_regulator,root_name,kind,max_depth_k,predicted_state,activation_z,overlap_p_value,network_bias_p_value,target_count,intermediate_regulators\n",
        );
        for (display_rank, (_, h)) in self.filtered_sorted_cna(rep).into_iter().enumerate() {
            let relays = h.intermediate_regulators.join("; ");
            let _ = writeln!(
                out,
                "{},{},{},{},{},{},{:.6},{:.6e},{:.6e},{},{}",
                display_rank + 1,
                csv_escape(&h.root_regulator),
                csv_escape(&h.root_name),
                csv_escape(h.kind.label()),
                h.max_depth_k,
                csv_escape(h.predicted_state.label()),
                h.activation_z,
                h.overlap_p_value,
                h.network_bias_p_value,
                h.target_count,
                csv_escape(&relays),
            );
        }
        out
    }
}

fn csv_escape(field: &str) -> String {
    if field.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}

fn sanitize_filename_stem(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect();
    let collapsed: Vec<&str> = cleaned.split('_').filter(|s| !s.is_empty()).collect();
    if collapsed.is_empty() {
        "cohort".to_string()
    } else {
        collapsed.join("_")
    }
}

fn render_sort_header<C: Copy + Eq>(
    ui: &mut Ui,
    label: &str,
    col: C,
    default_dir: SortDir,
    sort: &mut Option<TableSort<C>>,
) {
    let is_active = sort.is_some_and(|s| s.col == col);
    let text = match *sort {
        Some(s) if s.col == col => format!("{label} {}", s.dir.glyph()),
        _ => label.to_string(),
    };
    if ui
        .selectable_label(is_active, RichText::new(text).strong())
        .on_hover_text(format!("Click to sort by {label}"))
        .clicked()
    {
        TableSort::toggle_on(sort, col, default_dir);
    }
}

fn compute_nice_ticks(min_v: f64, max_v: f64, max_ticks: usize) -> Vec<f64> {
    let span = (max_v - min_v).abs();
    if span <= 1e-9 || max_ticks < 2 {
        return vec![0.0];
    }
    let raw_step = span / (max_ticks as f64);
    let mag = 10.0_f64.powf(raw_step.log10().floor());
    let norm = raw_step / mag;
    let nice_norm = if norm <= 1.0 {
        1.0
    } else if norm <= 2.0 {
        2.0
    } else if norm <= 5.0 {
        5.0
    } else {
        10.0
    };
    let step = (nice_norm * mag).max(0.5);
    let start_k = (min_v / step).ceil() as i64;
    let end_k = (max_v / step).floor() as i64;
    (start_k..=end_k).map(|k| (k as f64) * step).collect()
}

fn draw_pathway_compartment_bands(
    painter: &egui::Painter,
    rect: Rect,
    pm_y: f32,
    org_y: f32,
    nm_y: f32,
    zoom: f32,
    mode: WorkbenchPaletteMode,
) {
    let (ext_col, pm_col, cyto_col, org_col, nuc_col) = match mode {
        WorkbenchPaletteMode::AcademicDark => (
            Color32::from_rgba_unmultiplied(56, 189, 248, 12),
            Color32::from_rgba_unmultiplied(148, 163, 184, 22),
            Color32::from_rgba_unmultiplied(30, 41, 59, 18),
            Color32::from_rgba_unmultiplied(99, 102, 241, 12),
            Color32::from_rgba_unmultiplied(139, 92, 246, 16),
        ),
        WorkbenchPaletteMode::PublicationLight => (
            Color32::from_rgba_unmultiplied(191, 219, 254, 48),
            Color32::from_rgba_unmultiplied(148, 163, 184, 34),
            Color32::from_rgba_unmultiplied(226, 232, 240, 44),
            Color32::from_rgba_unmultiplied(199, 210, 254, 42),
            Color32::from_rgba_unmultiplied(221, 214, 254, 46),
        ),
    };

    let fill_band = |y0: f32, y1: f32, col: Color32| {
        let top = y0.max(rect.top());
        let bot = y1.min(rect.bottom());
        if bot > top + 1.0 {
            let r = Rect::from_min_max(Pos2::new(rect.left(), top), Pos2::new(rect.right(), bot));
            painter.rect_filled(r, 0.0, col);
        }
    };

    let pm_half = (12.0 * zoom).clamp(5.0, 18.0);
    fill_band(rect.top(), pm_y - pm_half, ext_col);
    fill_band(pm_y - pm_half, pm_y + pm_half, pm_col);
    fill_band(pm_y + pm_half, org_y, cyto_col);
    fill_band(org_y, nm_y, org_col);
    fill_band(nm_y, rect.bottom(), nuc_col);
}

fn draw_pathway_edge_legend(
    painter: &egui::Painter,
    rect: Rect,
    pal: &ScientificPalette,
    mode: WorkbenchPaletteMode,
) {
    let leg_w = (rect.width() - 200.0).clamp(240.0, 356.0);
    if rect.width() < 280.0 || rect.height() < 120.0 {
        return;
    }
    let leg_rect = Rect::from_min_size(
        Pos2::new(rect.left() + 10.0, rect.bottom() - 32.0),
        Vec2::new(leg_w, 22.0),
    );
    let bg = match mode {
        WorkbenchPaletteMode::AcademicDark => Color32::from_rgba_unmultiplied(15, 20, 29, 218),
        WorkbenchPaletteMode::PublicationLight => Color32::from_rgba_unmultiplied(248, 250, 252, 232),
    };
    painter.rect_filled(leg_rect, 4.0, bg);
    painter.rect_stroke(leg_rect, 4.0, Stroke::new(1.0_f32, pal.grid_line));

    let cy = leg_rect.center().y;
    let mut x = leg_rect.left() + 8.0;

    // Activating arrow
    painter.line_segment(
        [Pos2::new(x, cy), Pos2::new(x + 14.0, cy)],
        Stroke::new(1.8_f32, pal.activated_strong),
    );
    painter.add(egui::Shape::convex_polygon(
        vec![
            Pos2::new(x + 17.0, cy),
            Pos2::new(x + 11.0, cy - 3.5),
            Pos2::new(x + 11.0, cy + 3.5),
        ],
        pal.activated_strong,
        Stroke::NONE,
    ));
    painter.text(
        Pos2::new(x + 21.0, cy),
        Align2::LEFT_CENTER,
        "Activating",
        FontId::proportional(9.5),
        pal.text_primary,
    );
    x += 78.0;

    // Inhibiting T-bar
    painter.line_segment(
        [Pos2::new(x, cy), Pos2::new(x + 14.0, cy)],
        Stroke::new(1.8_f32, pal.inhibited_strong),
    );
    painter.line_segment(
        [Pos2::new(x + 14.0, cy - 4.5), Pos2::new(x + 14.0, cy + 4.5)],
        Stroke::new(2.0_f32, pal.inhibited_strong),
    );
    painter.text(
        Pos2::new(x + 19.0, cy),
        Align2::LEFT_CENTER,
        "Inhibiting",
        FontId::proportional(9.5),
        pal.text_primary,
    );
    x += 75.0;

    // Ambiguous / Binding
    painter.line_segment(
        [Pos2::new(x, cy), Pos2::new(x + 12.0, cy)],
        Stroke::new(1.6_f32, pal.border_default),
    );
    painter.circle_filled(Pos2::new(x + 13.0, cy), 2.8, pal.border_default);
    painter.text(
        Pos2::new(x + 19.0, cy),
        Align2::LEFT_CENTER,
        "Binding",
        FontId::proportional(9.5),
        pal.text_primary,
    );
    x += 64.0;

    if x + 76.0 <= leg_rect.right() {
        painter.line_segment(
            [Pos2::new(x, cy), Pos2::new(x + 14.0, cy)],
            Stroke::new(1.8_f32, pal.inconsistent_amber),
        );
        painter.text(
            Pos2::new(x + 18.0, cy),
            Align2::LEFT_CENTER,
            "Inconsistent",
            FontId::proportional(9.5),
            pal.text_primary,
        );
    }
}

struct TargetGeneRef<'a> {
    sym: &'a str,
    fc: f64,
    sign: EdgeSign,
}

fn short_perturbagen_label(name: &str) -> String {
    let base = match name.find('(') {
        Some(idx) => {
            let prefix = name[..idx].trim();
            if prefix.is_empty() {
                name.trim()
            } else {
                prefix
            }
        }
        None => name.trim(),
    };
    const MAX_LEN: usize = 20;
    if base.chars().count() <= MAX_LEN {
        return base.to_string();
    }
    let mut out = String::new();
    for word in base.split_whitespace() {
        let candidate_len = if out.is_empty() {
            word.chars().count()
        } else {
            out.chars().count() + 1 + word.chars().count()
        };
        if candidate_len <= MAX_LEN {
            if !out.is_empty() {
                out.push(' ');
            }
            out.push_str(word);
        } else if out.is_empty() {
            out = word.chars().take(MAX_LEN).collect();
            break;
        } else {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drug_rescue_short_label_strips_parenthetical_cleanly() {
        assert_eq!(
            short_perturbagen_label("KLF2 Lentiviral Overexpression"),
            "KLF2 Lentiviral"
        );
        assert_eq!(
            short_perturbagen_label("BMS-345541 (Selective IKK-beta Inhibitor)"),
            "BMS-345541"
        );
        assert_eq!(
            short_perturbagen_label("TPCA-1 (IKK-2 / NF-kB Inhibitor)"),
            "TPCA-1"
        );
        assert_eq!(short_perturbagen_label("Dexamethasone"), "Dexamethasone");
    }

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
            let out = ctx.run(input, |c| app.render_ui_frame(c));
            assert!(!out.shapes.is_empty(), "View {view} must emit shapes");
        }

        // Render all 5 tabs and subtabs headless (including Volcano Plot & CNA table & hidden minimap)
        for tab in WorkbenchTab::ALL {
            app.active_tab = tab;
            if tab == WorkbenchTab::PathwayAtlas {
                app.pathway_minimap_visible = false;
            }
            let input = egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                ..Default::default()
            };
            let out = ctx.run(input, |c| app.render_ui_frame(c));
            assert!(!out.shapes.is_empty(), "Tab {:?} must emit shapes", tab);
        }
        app.pathway_minimap_visible = true;
        app.active_tab = WorkbenchTab::CausalRegulators;
        app.causal_subtab = CausalRegulatorsSubTab::MultiHopNetworks;
        let out_cna = ctx.run(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                ..Default::default()
            },
            |c| app.render_ui_frame(c),
        );
        assert!(!out_cna.shapes.is_empty());

        // Verify Copyright & About modal toggle
        assert!(!app.about_open);
        app.apply_action(&ctx, GuiAction::ToggleAbout);
        assert!(app.about_open);
        let input_about = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
            ..Default::default()
        };
        let about_out = ctx.run(input_about, |c| app.render_ui_frame(c));
        assert!(!about_out.shapes.is_empty());
        app.apply_action(&ctx, GuiAction::ToggleAbout);
        assert!(!app.about_open);

        // Verify 1-click Drug Rescue & Active Clamps management
        app.configure_screenshot_view(&ctx, 1);
        app.apply_action(
            &ctx,
            GuiAction::ApplyDrugRescue {
                hit_idx: 0,
                switch_tab: true,
            },
        );
        assert_eq!(app.active_tab, WorkbenchTab::PathwayAtlas);
        assert!(!app.map_user_clamps.is_empty());
        assert!(app.active_drug_rescue_name.is_some());

        let first_clamp = app.map_user_clamps.keys().next().unwrap().clone();
        let initial_val = app.map_user_clamps[&first_clamp];
        app.apply_action(
            &ctx,
            GuiAction::SetNodeClamp(first_clamp.clone(), Some(-initial_val)),
        );
        assert_eq!(app.map_user_clamps[&first_clamp], -initial_val);
        app.apply_action(&ctx, GuiAction::SetNodeClamp(first_clamp.clone(), None));
        assert!(!app.map_user_clamps.contains_key(&first_clamp));

        let input_clamped = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
            ..Default::default()
        };
        let _ = ctx.run(input_clamped, |c| app.render_ui_frame(c));

        app.apply_action(&ctx, GuiAction::ClearAllClamps);
        assert!(app.map_user_clamps.is_empty());
        assert!(app.active_drug_rescue_name.is_none());

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }

    #[test]
    fn batch_finished_surfaces_errors_in_status_banner_and_sync_logs() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "causal_atlas_batch_err_test_{}.sqlite3",
            std::process::id()
        ));
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());
        app.worker_busy = true;

        let summary = BatchSyncSummary {
            streams_checked: 5,
            cache_hits: 3,
            streams_downloaded: 1,
            pathways_ingested: 12,
            causal_edges_ingested: 40,
            errors: vec!["WikiPathways GPML index timed out".to_string()],
        };
        app.tx
            .send(WorkerMessage::BatchFinished(summary))
            .expect("send BatchFinished");
        app.poll_worker();

        assert!(!app.worker_busy);
        assert!(app.status_banner.contains("WikiPathways GPML index timed out"));
        assert!(app
            .sync_logs
            .iter()
            .any(|line| line.contains("[ERROR] WikiPathways GPML index timed out")));

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }

    #[test]
    fn headless_f1_and_layered_escape_dismissal_priority() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "causal_atlas_shortcut_test_{}.sqlite3",
            std::process::id()
        ));
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());

        let step_key = |app: &mut CausalAtlasApp, key: egui::Key| {
            let input = egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                events: vec![egui::Event::Key {
                    key,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::NONE,
                }],
                ..Default::default()
            };
            let _ = ctx.run(input, |c| app.render_ui_frame(c));
        };

        app.omnibox_dropdown_open = true;
        app.pathway_fullscreen = true;
        app.pathway_fit_pending = false;
        app.pathway_pan = Vec2::new(123.0, -77.0);
        app.pathway_zoom = 2.75;
        assert!(!app.about_open);

        // F1 toggles About open
        step_key(&mut app, egui::Key::F1);
        assert!(app.about_open);

        // 1st Escape closes omnibox_dropdown_open first
        step_key(&mut app, egui::Key::Escape);
        assert!(!app.omnibox_dropdown_open);
        assert!(app.about_open);
        assert!(app.pathway_fullscreen);
        assert_eq!(app.pathway_pan, Vec2::new(123.0, -77.0));

        // 2nd Escape closes about_open next
        step_key(&mut app, egui::Key::Escape);
        assert!(!app.about_open);
        assert!(app.pathway_fullscreen);
        assert_eq!(app.pathway_pan, Vec2::new(123.0, -77.0));

        // 3rd Escape exits pathway_fullscreen and triggers pathway_fit_pending refit
        step_key(&mut app, egui::Key::Escape);
        assert!(!app.pathway_fullscreen);
        assert_eq!(app.pathway_pan, Vec2::ZERO);
        assert!((app.pathway_zoom - 2.75).abs() > 1e-3);

        // On another tab, exiting pathway_fullscreen leaves pathway_fit_pending = true until returning
        app.active_tab = WorkbenchTab::CohortVolcano;
        app.pathway_fullscreen = true;
        app.pathway_fit_pending = false;
        step_key(&mut app, egui::Key::Escape);
        assert!(!app.pathway_fullscreen);
        assert!(app.pathway_fit_pending);

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }

    #[test]
    fn keyboard_shortcuts_and_text_focus_guard_work_headless() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "causal_atlas_hotkeys_test_{}.sqlite3",
            std::process::id()
        ));
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());
        app.map_user_clamps.insert("TNF".to_string(), 1.0);
        app.pathway_fit_pending = false;

        let run_key = |app: &mut CausalAtlasApp, key: egui::Key, modifiers: egui::Modifiers| {
            let input = egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                events: vec![egui::Event::Key {
                    key,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers,
                }],
                ..Default::default()
            };
            let _ = ctx.run(input, |c| app.render_ui_frame(c));
        };

        // Keys 1..=5 switch across the 5 workbench pillars
        run_key(&mut app, egui::Key::Num1, egui::Modifiers::NONE);
        assert_eq!(app.active_tab, WorkbenchTab::CohortVolcano);
        run_key(&mut app, egui::Key::Num3, egui::Modifiers::NONE);
        assert_eq!(app.active_tab, WorkbenchTab::CausalRegulators);
        run_key(&mut app, egui::Key::Num4, egui::Modifiers::NONE);
        assert_eq!(app.active_tab, WorkbenchTab::PhenotypeDrugMatch);
        run_key(&mut app, egui::Key::Num5, egui::Modifiers::NONE);
        assert_eq!(app.active_tab, WorkbenchTab::KnowledgeDatabases);
        run_key(&mut app, egui::Key::Num2, egui::Modifiers::NONE);
        assert_eq!(app.active_tab, WorkbenchTab::PathwayAtlas);

        // R clears MAP clamps
        assert!(!app.map_user_clamps.is_empty());
        run_key(&mut app, egui::Key::R, egui::Modifiers::NONE);
        assert!(app.map_user_clamps.is_empty());

        // F and 0 trigger pathway camera fit
        app.pathway_zoom = 2.5;
        app.pathway_pan = Vec2::new(90.0, -45.0);
        run_key(&mut app, egui::Key::F, egui::Modifiers::NONE);
        assert!(app.pathway_fit_pending);
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                ..Default::default()
            },
            |c| app.render_ui_frame(c),
        );
        assert_eq!(app.pathway_pan, Vec2::ZERO);

        app.pathway_pan = Vec2::new(60.0, 30.0);
        run_key(&mut app, egui::Key::Num0, egui::Modifiers::NONE);
        assert!(app.pathway_fit_pending);
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                ..Default::default()
            },
            |c| app.render_ui_frame(c),
        );
        assert_eq!(app.pathway_pan, Vec2::ZERO);

        // Text-input guard: when a TextEdit has real focus (wants_keyboard_input() == true),
        // 1-5/F/0/R/Ctrl+E must NOT trigger tab switch, clamp reset, or export
        app.active_tab = WorkbenchTab::PathwayAtlas;
        app.map_user_clamps.insert("IKBKB".to_string(), -1.0);
        app.pathway_fit_pending = false;
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                modifiers: egui::Modifiers::CTRL,
                events: vec![
                    egui::Event::Key {
                        key: egui::Key::Num1,
                        physical_key: None,
                        pressed: true,
                        repeat: false,
                        modifiers: egui::Modifiers::NONE,
                    },
                    egui::Event::Key {
                        key: egui::Key::R,
                        physical_key: None,
                        pressed: true,
                        repeat: false,
                        modifiers: egui::Modifiers::NONE,
                    },
                    egui::Event::Key {
                        key: egui::Key::E,
                        physical_key: None,
                        pressed: true,
                        repeat: false,
                        modifiers: egui::Modifiers::CTRL,
                    },
                ],
                ..Default::default()
            },
            |c| {
                let mut guarded_actions = Vec::new();
                app.process_keyboard_shortcuts(c, true, &mut guarded_actions);
                assert!(guarded_actions.is_empty());
            },
        );
        assert_eq!(app.active_tab, WorkbenchTab::PathwayAtlas);
        assert!(!app.map_user_clamps.is_empty());

        // Ctrl+E when NOT in text input queues ExportPathwaySvgDialog
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1600.0, 1000.0))),
                modifiers: egui::Modifiers::CTRL,
                events: vec![egui::Event::Key {
                    key: egui::Key::E,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::CTRL,
                }],
                ..Default::default()
            },
            |c| {
                let mut acts = Vec::new();
                app.process_keyboard_shortcuts(c, false, &mut acts);
                assert!(
                    matches!(acts.as_slice(), [GuiAction::ExportPathwaySvgDialog]),
                    "Unexpected acts: {}",
                    acts.len()
                );
            },
        );

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }

    #[test]
    fn tables_sort_filter_csv_export_and_minimap_toggle_work() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "causal_atlas_tables_test_{}.sqlite3",
            std::process::id()
        ));
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());

        // Minimap toggle action flips visibility cleanly
        assert!(app.pathway_minimap_visible);
        app.apply_action(&ctx, GuiAction::TogglePathwayMinimap);
        assert!(!app.pathway_minimap_visible);
        app.apply_action(&ctx, GuiAction::TogglePathwayMinimap);
        assert!(app.pathway_minimap_visible);

        // Gene filtering & sorting + CSV export
        app.gene_dir_filter = GeneDirectionFilter::DownOnly;
        app.gene_sort = Some(TableSort {
            col: GeneSortCol::Log2Fc,
            dir: SortDir::Asc,
        });
        let genes = app.filtered_sorted_genes();
        assert!(!genes.is_empty());
        assert!(genes.iter().all(|g| g.log2_fc < 0.0));
        assert!(genes.windows(2).all(|w| w[0].log2_fc <= w[1].log2_fc));
        let gene_csv = app.export_genes_csv_string();
        assert!(gene_csv.starts_with("symbol,log2_fc,p_value,bh_fdr,passes_cutoffs\n"));
        assert!(gene_csv.contains("KLF2"));

        // URA filtering & sorting + CSV export
        app.ura_state_filter = RegStateFilter::Inhibited;
        app.ura_sort = Some(TableSort {
            col: UraSortCol::ActZ,
            dir: SortDir::Asc,
        });
        let rep = app.reports[0].clone();
        let ura_rows = app.filtered_sorted_ura(&rep);
        assert!(!ura_rows.is_empty());
        assert!(ura_rows
            .iter()
            .all(|(_, h)| h.predicted_state == PredictedState::Inhibited));
        let ura_csv = app.export_ura_csv_string();
        assert!(ura_csv.contains("DEXAMETHASONE"));
        assert!(!ura_csv.contains(",TNF,"));

        // CNA filtering & sorting + CSV export
        app.cna_state_filter = RegStateFilter::Activated;
        app.cna_sort = Some(TableSort {
            col: CnaSortCol::NetBiasP,
            dir: SortDir::Asc,
        });
        let cna_rows = app.filtered_sorted_cna(&rep);
        assert!(!cna_rows.is_empty());
        assert!(cna_rows
            .windows(2)
            .all(|w| w[0].1.network_bias_p_value <= w[1].1.network_bias_p_value));
        let cna_csv = app.export_cna_csv_string();
        assert!(cna_csv.starts_with("rank,root_regulator,root_name,kind,max_depth_k,"));

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }

    #[test]
    fn ctrl_k_command_palette_filters_switches_tabs_and_computes_rank_diff() {
        let ctx = egui::Context::default();
        let db_path = std::env::temp_dir().join(format!(
            "ca_palette_test_{}.sqlite3",
            std::process::id()
        ));
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
        let mut app = CausalAtlasApp::new_headless(&ctx, db_path.clone());
        assert!(!app.show_command_palette);
        assert!(ca_fuzzy_palette_score("diff", "Export Regulator Rank Diff (CSV)", "diff").is_some());
        assert!(ca_fuzzy_palette_score("zzzz", "Fit Pathway Camera", "fit").is_none());

        let send_key = |app: &mut CausalAtlasApp, key: egui::Key, ctrl: bool| {
            let mods = if ctrl {
                egui::Modifiers::CTRL
            } else {
                egui::Modifiers::NONE
            };
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                Pos2::ZERO,
                Vec2::new(1280.0, 800.0),
            ));
            raw.modifiers = mods;
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: mods,
            });
            let _ = ctx.run(raw, |c| app.render_ui_frame(c));
        };

        send_key(&mut app, egui::Key::K, true);
        assert!(app.show_command_palette);

        app.command_palette_query = "multi-cohort comparison".to_string();
        send_key(&mut app, egui::Key::Enter, false);
        assert!(!app.show_command_palette);
        assert_eq!(app.active_tab, WorkbenchTab::PhenotypeDrugMatch);
        assert_eq!(
            app.phenotype_subtab,
            PhenotypeDrugSubTab::ComparisonAndReversers
        );

        let (ra, rb) = app.diff_report_pair().expect("built-in cohorts pair");
        let digest = compute_regulator_rank_diff(ra, rb);
        assert!(!digest.entries.is_empty());

        send_key(&mut app, egui::Key::K, true);
        assert!(app.show_command_palette);
        send_key(&mut app, egui::Key::Escape, false);
        assert!(!app.show_command_palette);

        drop(app);
        for suffix in ["", "-wal", "-shm"] {
            let _ = std::fs::remove_file(format!("{}{}", db_path.display(), suffix));
        }
    }
}
