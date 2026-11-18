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

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#![deny(unsafe_code)]

mod gui;
mod icon;
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

    /// Render a 1600x1000 offscreen GPU screenshot to a PNG file without opening a GUI window
    #[arg(long, global = true, value_name = "OUT_PNG")]
    screenshot: Option<PathBuf>,

    /// Select which workbench view to capture with --screenshot (1, 2, or 3)
    #[arg(long, global = true, default_value_t = 1, value_name = "VIEW")]
    screenshot_view: u8,

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

fn block_on<F: std::future::Future>(future: F) -> F::Output {
    use std::sync::Arc;
    use std::task::{Context, Poll, Wake, Waker};

    struct ThreadWaker(std::thread::Thread);
    impl Wake for ThreadWaker {
        fn wake(self: Arc<Self>) {
            self.0.unpark();
        }
        fn wake_by_ref(self: &Arc<Self>) {
            self.0.unpark();
        }
    }

    let waker = Waker::from(Arc::new(ThreadWaker(std::thread::current())));
    let mut cx = Context::from_waker(&waker);
    let mut future = std::pin::pin!(future);
    loop {
        match future.as_mut().poll(&mut cx) {
            Poll::Ready(out) => return out,
            Poll::Pending => std::thread::park(),
        }
    }
}

fn write_rgba_png(
    path: &std::path::Path,
    width: u32,
    height: u32,
    rgba: &[u8],
) -> std::io::Result<()> {
    fn crc32(data: &[u8]) -> u32 {
        let mut crc = 0xFFFF_FFFFu32;
        for &b in data {
            crc ^= b as u32;
            for _ in 0..8 {
                crc = if (crc & 1) != 0 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
            }
        }
        !crc
    }
    fn push_chunk(out: &mut Vec<u8>, tag: &[u8; 4], payload: &[u8]) {
        out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
        let start = out.len();
        out.extend_from_slice(tag);
        out.extend_from_slice(payload);
        out.extend_from_slice(&crc32(&out[start..]).to_be_bytes());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let row_bytes = width as usize * 4;
    let mut raw = Vec::with_capacity(height as usize * (1 + row_bytes));
    for row in rgba.chunks_exact(row_bytes) {
        raw.push(1); // PNG Sub filter
        for (i, &b) in row.iter().enumerate() {
            let left = if i >= 4 { row[i - 4] } else { 0 };
            raw.push(b.wrapping_sub(left));
        }
    }
    let zlib = {
        use std::io::Write;
        let mut enc =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::best());
        enc.write_all(&raw)?;
        enc.finish()?
    };
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.extend_from_slice(&[8, 6, 0, 0, 0]);
    let mut png = Vec::with_capacity(zlib.len() + 64);
    png.extend_from_slice(b"\x89PNG\r\n\x1a\n");
    push_chunk(&mut png, b"IHDR", &ihdr);
    push_chunk(&mut png, b"IDAT", &zlib);
    push_chunk(&mut png, b"IEND", &[]);
    std::fs::write(path, png)
}

fn capture_headless_screenshot(
    out_png: &std::path::Path,
    db_path: &std::path::Path,
    view: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    use eframe::{egui_wgpu, wgpu};

    const W: u32 = 1600;
    const H: u32 = 1000;

    let use_temp_db = db_path == std::path::Path::new("causal_atlas.sqlite3") || !db_path.exists();
    let effective_db = if db_path == std::path::Path::new("causal_atlas.sqlite3") {
        std::env::temp_dir().join(format!(
            "causal_atlas_screenshot_{}.sqlite3",
            std::process::id()
        ))
    } else {
        db_path.to_path_buf()
    };

    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::all(),
        ..Default::default()
    });
    let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .ok_or("No suitable wgpu adapter found for offscreen capture")?;

    let (device, queue) = block_on(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("causal_atlas_headless_device"),
            required_features: wgpu::Features::empty(),
            required_limits: wgpu::Limits::default(),
            memory_hints: wgpu::MemoryHints::default(),
        },
        None,
    ))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let mut renderer = egui_wgpu::Renderer::new(&device, format, None, 1, false);
    let tex_extent = wgpu::Extent3d {
        width: W,
        height: H,
        depth_or_array_layers: 1,
    };
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("causal_atlas_offscreen_target"),
        size: tex_extent,
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let tex_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let bytes_per_row = W * 4;
    let readback_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("causal_atlas_offscreen_readback"),
        size: (bytes_per_row * H) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let ctx = egui::Context::default();
    ctx.set_pixels_per_point(1.0);
    let mut app = CausalAtlasApp::new_headless(&ctx, effective_db.clone());
    app.configure_screenshot_view(&ctx, view);

    for frame_idx in 0..4u32 {
        let mut raw_input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(W as f32, H as f32),
            )),
            time: Some(frame_idx as f64 * 0.016),
            predicted_dt: 1.0 / 60.0,
            ..Default::default()
        };
        raw_input
            .viewports
            .entry(egui::ViewportId::ROOT)
            .or_default()
            .native_pixels_per_point = Some(1.0);

        let full_output = ctx.run(raw_input, |ctx| app.render_ui_frame(ctx));
        let clipped_primitives =
            ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [W, H],
            pixels_per_point: 1.0,
        };

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("causal_atlas_offscreen_encoder"),
        });
        for (id, image_delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, image_delta);
        }
        let user_cmds = renderer.update_buffers(
            &device,
            &queue,
            &mut encoder,
            &clipped_primitives,
            &screen_descriptor,
        );

        let bg = ctx.style().visuals.panel_fill.to_normalized_gamma_f32();
        {
            let mut rpass = encoder
                .begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("causal_atlas_offscreen_pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &tex_view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color {
                                r: bg[0] as f64,
                                g: bg[1] as f64,
                                b: bg[2] as f64,
                                a: 1.0,
                            }),
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    occlusion_query_set: None,
                    timestamp_writes: None,
                })
                .forget_lifetime();
            renderer.render(&mut rpass, &clipped_primitives, &screen_descriptor);
        }

        if frame_idx == 3 {
            encoder.copy_texture_to_buffer(
                texture.as_image_copy(),
                wgpu::ImageCopyBuffer {
                    buffer: &readback_buf,
                    layout: wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(bytes_per_row),
                        rows_per_image: Some(H),
                    },
                },
                tex_extent,
            );
        }

        queue.submit(user_cmds.into_iter().chain(std::iter::once(encoder.finish())));
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
    }

    let slice = readback_buf.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()??;
    let rgba = slice.get_mapped_range().to_vec();
    readback_buf.unmap();

    drop(app);
    if use_temp_db {
        for suffix in ["", "-wal", "-shm"] {
            let p = PathBuf::from(format!("{}{}", effective_db.display(), suffix));
            let _ = std::fs::remove_file(p);
        }
    }

    write_rgba_png(out_png, W, H, &rgba)?;
    println!(
        "Saved offscreen GPU screenshot ({}x{} RGBA) to {}",
        W,
        H,
        out_png.display()
    );
    Ok(())
}

#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_parent_console_if_cli() {
    if std::env::args_os().len() > 1 {
        const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
        extern "system" {
            fn AttachConsole(dw_process_id: u32) -> i32;
        }
        // SAFETY: Attaches stdout/stderr to an existing parent console only when CLI
        // arguments are supplied; failure when no parent console exists is harmless.
        unsafe {
            let _ = AttachConsole(ATTACH_PARENT_PROCESS);
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(windows)]
    attach_parent_console_if_cli();

    let cli = Cli::parse();
    let db_path = cli.db;

    if let Some(out_png) = cli.screenshot {
        return capture_headless_screenshot(&out_png, &db_path, cli.screenshot_view);
    }

    match cli.command.unwrap_or(Commands::Gui) {
        Commands::Gui => {
            let native_options = eframe::NativeOptions {
                viewport: egui::ViewportBuilder::default()
                    .with_title("Causal Atlas — Open-Source Pathway & Causal Network Workbench (SQLite WAL)")
                    .with_icon(std::sync::Arc::new(icon::app_icon_data()))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_rgba_png_uses_deflate_compression_under_1mb() {
        let w = 1600u32;
        let h = 1000u32;
        let mut rgba = vec![0u8; (w * h * 4) as usize];
        for y in 0..h {
            for x in 0..w {
                let idx = ((y * w + x) * 4) as usize;
                rgba[idx] = ((x / 32) as u8).wrapping_mul(17);
                rgba[idx + 1] = ((y / 25) as u8).wrapping_mul(19);
                rgba[idx + 2] = 48;
                rgba[idx + 3] = 255;
            }
        }
        let out_path = std::env::temp_dir().join(format!(
            "causal_atlas_deflate_test_{}.png",
            std::process::id()
        ));
        write_rgba_png(&out_path, w, h, &rgba).expect("write PNG");
        let bytes = std::fs::read(&out_path).expect("read PNG");
        let _ = std::fs::remove_file(&out_path);

        assert!(bytes.starts_with(b"\x89PNG\r\n\x1a\n"));
        assert!(bytes.windows(4).any(|w| w == b"IDAT"));
        assert!(bytes.windows(4).any(|w| w == b"IEND"));
        assert!(
            bytes.len() < 1_048_576,
            "DEFLATE PNG should be < 1.0 MB, got {} bytes",
            bytes.len()
        );
    }
}
