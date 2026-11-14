#![forbid(unsafe_code)]

mod gui;
mod svg_export;
mod theme;
mod tui;

use ca_core::{
    builtin_cohorts, parse_expression_table, run_core_analysis, seed_kb::build_seed_knowledge_graph,
    simulate_map, CoreAnalysisConfig,
};
use ca_streams::{
    client::import_local_stream_file, run_batch_download_and_update, SqliteKnowledgeStore,
};
use clap::{Parser, Subcommand};
use gui::CausalAtlasApp;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "causal-atlas",
    version,
    about = "Independent native Rust causal pathway, upstream regulator, and network workbench with SQLite WAL caching"
)]
struct Cli {
    /// Path to the persistent SQLite 3 WAL knowledge & analysis cache database
    #[arg(long, global = true, default_value = "causal_atlas.sqlite3")]
    db: PathBuf,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Launch the native cross-platform desktop GUI (default when no subcommand is given)
    Gui,

    /// Batch-download & cache all human WikiPathways GPML diagrams, Reactome GMT pathways, and regulatory streams into SQLite
    DownloadAll {
        /// Bypass ETag / TTL cache freshness checks and force full remote re-download
        #[arg(long)]
        force: bool,
        /// Only download Canonical Pathway bundles (WikiPathways GPML zip + Reactome GMT zip)
        #[arg(long)]
        pathways_only: bool,
        /// Optional local file or `.zip` archive (`.gpml`, `.gmt`, `.tsv`, `.gaf.gz`) to import into SQLite
        #[arg(long)]
        local_file: Option<PathBuf>,
    },

    /// Launch the headless terminal UI (SSH/workstation compatible) backed by SQLite cache
    Tui {
        /// Optional path to a differential expression CSV/TSV file (defaults to built-in TNF benchmark)
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Execute headless Core Analysis using the SQLite knowledge & result cache
    Analyze {
        /// Path to input CSV/TSV differential expression file (`symbol`, `log2FC`, `pvalue`, `padj`)
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output JSON file path
        #[arg(short, long, default_value = "causal_atlas_report.json")]
        output: PathBuf,
        /// Minimum |log2FC| cutoff
        #[arg(long, default_value_t = 1.0)]
        fc_cutoff: f64,
        /// Maximum BH-FDR cutoff
        #[arg(long, default_value_t = 0.05)]
        fdr_cutoff: f64,
    },

    /// Export vector SVG diagrams (Canonical Pathway MAP & Downstream TreeMap) for a cohort
    ExportSvg {
        /// Input expression CSV/TSV (defaults to built-in TNF benchmark)
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output directory for generated `.svg` diagrams
        #[arg(short, long, default_value = ".")]
        out_dir: PathBuf,
    },
}

fn open_seeded_db(path: &PathBuf) -> Result<SqliteKnowledgeStore, Box<dyn std::error::Error>> {
    let mut store = SqliteKnowledgeStore::open(path)?;
    store.ensure_seeded(build_seed_knowledge_graph)?;
    Ok(store)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let db_path = cli.db;

    match cli.command.unwrap_or(Commands::Gui) {
        Commands::Gui => {
            let native_options = eframe::NativeOptions {
                viewport: egui::ViewportBuilder::default()
                    .with_title("Causal Atlas — Open-Source Pathway & Causal Network Workbench (SQLite WAL)")
                    .with_inner_size([1500.0, 920.0])
                    .with_min_inner_size([1024.0, 680.0]),
                ..Default::default()
            };
            eframe::run_native(
                "Causal Atlas",
                native_options,
                Box::new(move |cc| Ok(Box::new(CausalAtlasApp::new(cc, db_path)))),
            )
            .map_err(|e| format!("GUI error: {e}"))?;
        }
        Commands::DownloadAll {
            force,
            pathways_only,
            local_file,
        } => {
            let mut store = open_seeded_db(&db_path)?;
            if let Some(lf) = local_file {
                let (delta, cnt, sha) = import_local_stream_file(&lf, None)?;
                store.upsert_knowledge_graph(&delta)?;
                println!(
                    "Imported {cnt} records from {} into {} (SHA-256: {})",
                    lf.display(),
                    db_path.display(),
                    &sha[..16]
                );
            } else {
                println!(
                    "Starting SQLite WAL batch downloader -> {} (force={force})",
                    db_path.display()
                );
                let summary = run_batch_download_and_update(
                    &mut store,
                    true,
                    !pathways_only,
                    force,
                    None,
                    |ev| {
                        println!(
                            "[{}/{}] {} — {}",
                            ev.step_index, ev.total_steps, ev.source_name, ev.message
                        );
                    },
                );
                let stats = store.stats()?;
                println!(
                    "Batch sync complete: {} checked, {} cache hits, {} downloaded | SQLite totals: {} pathways ({} with 2D coords), {} causal edges ({:.1} KB)",
                    summary.streams_checked,
                    summary.cache_hits,
                    summary.streams_downloaded,
                    stats.cached_pathways,
                    stats.cached_diagrams_with_coords,
                    stats.cached_edges,
                    stats.file_size_bytes as f64 / 1024.0
                );
            }
        }
        Commands::Tui { input } => {
            let store = open_seeded_db(&db_path)?;
            let graph = store.load_knowledge_graph()?;
            let dataset = if let Some(in_path) = input {
                let raw = std::fs::read_to_string(&in_path)?;
                let name = in_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("Input");
                parse_expression_table(name, &raw)?
            } else {
                builtin_cohorts().into_iter().next().unwrap()
            };
            let cfg = CoreAnalysisConfig::default();
            let report = if let Some(cached) = store.get_cached_analysis(&dataset, &cfg)? {
                cached
            } else {
                let r = run_core_analysis(&dataset, &graph, &cfg);
                store.put_cached_analysis(&dataset, &cfg, &r)?;
                r
            };
            tui::run_headless_tui(&report)?;
        }
        Commands::Analyze {
            input,
            output,
            fc_cutoff,
            fdr_cutoff,
        } => {
            let store = open_seeded_db(&db_path)?;
            let graph = store.load_knowledge_graph()?;
            let dataset = if let Some(in_path) = input {
                let raw = std::fs::read_to_string(&in_path)?;
                let name = in_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("Input");
                parse_expression_table(name, &raw)?
            } else {
                builtin_cohorts().into_iter().next().unwrap()
            };
            let cfg = CoreAnalysisConfig {
                abs_log2fc_cutoff: fc_cutoff,
                fdr_cutoff,
                ..CoreAnalysisConfig::default()
            };
            let (report, was_hit) = if let Some(cached) = store.get_cached_analysis(&dataset, &cfg)? {
                (cached, true)
            } else {
                let r = run_core_analysis(&dataset, &graph, &cfg);
                store.put_cached_analysis(&dataset, &cfg, &r)?;
                (r, false)
            };
            let json = serde_json::to_string_pretty(&report)?;
            std::fs::write(&output, json)?;
            println!(
                "Saved Core Analysis report [{}] ({} pathways, {} upstream regulators, {} multi-hop causal networks) to {}",
                if was_hit { "SQLite CACHE HIT" } else { "COMPUTED & CACHED" },
                report.canonical_pathways.len(),
                report.upstream_regulators.len(),
                report.causal_networks.len(),
                output.display()
            );
        }
        Commands::ExportSvg { input, out_dir } => {
            std::fs::create_dir_all(&out_dir)?;
            let store = open_seeded_db(&db_path)?;
            let graph = store.load_knowledge_graph()?;
            let dataset = if let Some(in_path) = input {
                let raw = std::fs::read_to_string(&in_path)?;
                let name = in_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("Input");
                parse_expression_table(name, &raw)?
            } else {
                builtin_cohorts().into_iter().next().unwrap()
            };
            let cfg = CoreAnalysisConfig::default();
            let report = if let Some(cached) = store.get_cached_analysis(&dataset, &cfg)? {
                cached
            } else {
                let r = run_core_analysis(&dataset, &graph, &cfg);
                store.put_cached_analysis(&dataset, &cfg, &r)?;
                r
            };

            let tm_svg = svg_export::render_treemap_svg(&report.dataset_name, &report.treemap_rects);
            let tm_path = out_dir.join("downstream_effects_treemap.svg");
            std::fs::write(&tm_path, tm_svg)?;
            println!("Wrote {}", tm_path.display());

            if let Some(top_pw) = report.canonical_pathways.first() {
                if let Some(pw_def) = graph.pathways.iter().find(|p| p.id == top_pw.pathway_id) {
                    let active =
                        dataset.filter_active(cfg.abs_log2fc_cutoff, cfg.fdr_cutoff, &graph);
                    let map_res = simulate_map(
                        pw_def,
                        &active,
                        &BTreeMap::new(),
                        cfg.map_damping,
                        cfg.map_iterations,
                    );
                    let pw_svg = svg_export::render_pathway_map_svg(pw_def, &map_res, &graph);
                    let pw_path = out_dir.join(format!("{}_pathway_map.svg", pw_def.id));
                    std::fs::write(&pw_path, pw_svg)?;
                    println!("Wrote {}", pw_path.display());
                }
            }
        }
    }

    Ok(())
}
