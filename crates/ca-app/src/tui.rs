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

use ca_core::CoreAnalysisReport;
use crossterm::{
    cursor,
    event::{self, Event, KeyCode, KeyEventKind},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs},
    Terminal,
};
use std::io;
use std::time::Duration;

fn restore_terminal_buffer<W: io::Write>(out: &mut W) -> io::Result<()> {
    if execute!(out, LeaveAlternateScreen, cursor::Show).is_err() {
        out.write_all(b"\x1b[?1049l\x1b[?25h")?;
        out.flush()?;
    }
    Ok(())
}

struct TuiSessionGuard {
    active: bool,
    cleanup: Option<Box<dyn FnOnce()>>,
}

impl TuiSessionGuard {
    fn enter() -> io::Result<Self> {
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        if let Err(err) = execute!(stdout, EnterAlternateScreen) {
            let _ = disable_raw_mode();
            return Err(err);
        }
        Ok(Self::from_cleanup(Box::new(|| {
            let _ = disable_raw_mode();
            let mut out = io::stdout();
            let _ = restore_terminal_buffer(&mut out);
        })))
    }

    fn from_cleanup(cleanup: Box<dyn FnOnce()>) -> Self {
        Self {
            active: true,
            cleanup: Some(cleanup),
        }
    }

    fn release(&mut self) {
        if self.active {
            self.active = false;
            if let Some(cleanup) = self.cleanup.take() {
                cleanup();
            }
        }
    }
}

impl Drop for TuiSessionGuard {
    fn drop(&mut self) {
        self.release();
    }
}

/// Launches the headless SSH-compatible terminal workbench (`causal-atlas tui`).
pub fn run_headless_tui(report: &CoreAnalysisReport) -> io::Result<()> {
    let mut guard = TuiSessionGuard::enter()?;
    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;

    let mut tab_idx: usize = 0;
    let tab_titles = [
        "1: Canonical Pathways",
        "2: Upstream Regulators",
        "3: Multi-Hop Causal Networks",
        "4: Perturbation Reversers",
    ];

    let loop_res = loop {
        terminal.draw(|f| {
            let chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(3),
                    Constraint::Length(3),
                    Constraint::Min(10),
                    Constraint::Length(2),
                ])
                .split(f.area());

            let header = Paragraph::new(Line::from(vec![
                Span::styled(
                    " CAUSAL ATLAS TUI ",
                    Style::default()
                        .fg(Color::Black)
                        .bg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw(format!(
                    "  Dataset: {}  |  Active DEGs: {} ({} up / {} down)  |  |log2FC| >= {:.2}, FDR <= {:.3}",
                    report.dataset_name,
                    report.summary.active_genes,
                    report.summary.upregulated_count,
                    report.summary.downregulated_count,
                    report.config.abs_log2fc_cutoff,
                    report.config.fdr_cutoff
                )),
            ]))
            .block(Block::default().borders(Borders::ALL));
            f.render_widget(header, chunks[0]);

            let tabs = Tabs::new(tab_titles.iter().copied().map(Line::from).collect::<Vec<_>>())
                .select(tab_idx)
                .highlight_style(
                    Style::default()
                        .fg(Color::Yellow)
                        .add_modifier(Modifier::BOLD),
                )
                .block(Block::default().borders(Borders::ALL).title("Views (Left/Right or 1-4)"));
            f.render_widget(tabs, chunks[1]);

            match tab_idx {
                0 => {
                    let rows = report.canonical_pathways.iter().take(25).map(|p| {
                        let z_str = p
                            .activation_z
                            .map(ca_core::format_z_2dp)
                            .unwrap_or_else(|| "n/a".to_string());
                        Row::new(vec![
                            Cell::from(p.pathway_id.as_str()),
                            Cell::from(p.pathway_name.as_str()),
                            Cell::from(p.source_db.as_str()),
                            Cell::from(format!("{}/{}", p.overlap_count, p.pathway_size)),
                            Cell::from(format!("{:.2e}", p.fdr)),
                            Cell::from(z_str),
                        ])
                    });
                    let tbl = Table::new(
                        rows,
                        [
                            Constraint::Length(14),
                            Constraint::Percentage(42),
                            Constraint::Length(14),
                            Constraint::Length(10),
                            Constraint::Length(12),
                            Constraint::Length(10),
                        ],
                    )
                    .header(
                        Row::new(["ID", "Pathway", "Source", "Overlap", "BH-FDR", "Act Z"])
                            .style(Style::default().add_modifier(Modifier::BOLD)),
                    )
                    .block(Block::default().borders(Borders::ALL).title("Canonical Pathways"));
                    f.render_widget(tbl, chunks[2]);
                }
                1 => {
                    let rows = report.upstream_regulators.iter().take(25).map(|r| {
                        Row::new(vec![
                            Cell::from(r.regulator_id.as_str()),
                            Cell::from(r.kind.label()),
                            Cell::from(r.predicted_state.label()),
                            Cell::from(ca_core::format_z_2dp(r.activation_z)),
                            Cell::from(ca_core::format_z_2dp(r.bias_corrected_z)),
                            Cell::from(format!("{:.2e}", r.p_value)),
                            Cell::from(format!("{}", r.overlap_count)),
                        ])
                    });
                    let tbl = Table::new(
                        rows,
                        [
                            Constraint::Length(14),
                            Constraint::Length(24),
                            Constraint::Length(12),
                            Constraint::Length(10),
                            Constraint::Length(10),
                            Constraint::Length(12),
                            Constraint::Length(10),
                        ],
                    )
                    .header(
                        Row::new([
                            "Regulator",
                            "Molecule Type",
                            "State",
                            "Act Z",
                            "Bias-Z",
                            "Overlap P",
                            "Targets",
                        ])
                        .style(Style::default().add_modifier(Modifier::BOLD)),
                    )
                    .block(
                        Block::default()
                            .borders(Borders::ALL)
                            .title("Upstream Regulator Analysis (Krämer et al., 2014)"),
                    );
                    f.render_widget(tbl, chunks[2]);
                }
                2 => {
                    let rows = report.causal_networks.iter().take(25).map(|c| {
                        Row::new(vec![
                            Cell::from(c.root_regulator.as_str()),
                            Cell::from(format!("K={}", c.max_depth_k)),
                            Cell::from(c.predicted_state.label()),
                            Cell::from(ca_core::format_z_2dp(c.activation_z)),
                            Cell::from(format!("{:.3e}", c.network_bias_p_value)),
                            Cell::from(c.intermediate_regulators.join(", ")),
                        ])
                    });
                    let tbl = Table::new(
                        rows,
                        [
                            Constraint::Length(14),
                            Constraint::Length(6),
                            Constraint::Length(12),
                            Constraint::Length(10),
                            Constraint::Length(14),
                            Constraint::Percentage(44),
                        ],
                    )
                    .header(
                        Row::new([
                            "Root Master",
                            "Depth",
                            "State",
                            "Act Z",
                            "Net-Bias P",
                            "Intermediates",
                        ])
                        .style(Style::default().add_modifier(Modifier::BOLD)),
                    )
                    .block(
                        Block::default()
                            .borders(Borders::ALL)
                            .title("Multi-Hop Causal Networks"),
                    );
                    f.render_widget(tbl, chunks[2]);
                }
                _ => {
                    let rows = report.perturbation_matches.iter().take(25).map(|p| {
                        Row::new(vec![
                            Cell::from(p.perturbagen.as_str()),
                            Cell::from(p.kind.label()),
                            Cell::from(format!("{:+.1}", p.connectivity_score)),
                            Cell::from(p.mode.as_str()),
                            Cell::from(p.mechanism_of_action.as_str()),
                        ])
                    });
                    let tbl = Table::new(
                        rows,
                        [
                            Constraint::Length(28),
                            Constraint::Length(18),
                            Constraint::Length(10),
                            Constraint::Length(28),
                            Constraint::Percentage(35),
                        ],
                    )
                    .header(
                        Row::new(["Perturbagen", "Class", "Score", "Direction", "Mechanism"])
                            .style(Style::default().add_modifier(Modifier::BOLD)),
                    )
                    .block(
                        Block::default()
                            .borders(Borders::ALL)
                            .title("SigCom LINCS Open Perturbation Connectivity"),
                    );
                    f.render_widget(tbl, chunks[2]);
                }
            }

            let footer = Paragraph::new(
                " [Tab / Left / Right / 1-4] Switch Module   [q / Esc] Quit Headless TUI",
            );
            f.render_widget(footer, chunks[3]);
        })?;

        if event::poll(Duration::from_millis(150))? {
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char('q') | KeyCode::Esc => break Ok(()),
                        KeyCode::Right | KeyCode::Tab => {
                            tab_idx = (tab_idx + 1) % tab_titles.len();
                        }
                        KeyCode::Left => {
                            tab_idx = (tab_idx + tab_titles.len() - 1) % tab_titles.len();
                        }
                        KeyCode::Char('1') => tab_idx = 0,
                        KeyCode::Char('2') => tab_idx = 1,
                        KeyCode::Char('3') => tab_idx = 2,
                        KeyCode::Char('4') => tab_idx = 3,
                        _ => {}
                    }
                }
            }
        }
    };

    drop(terminal);
    guard.release();
    loop_res
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    };

    #[test]
    fn tui_session_guard_restores_on_early_error_return_and_is_idempotent() {
        let calls = Arc::new(AtomicUsize::new(0));
        let calls_ref = Arc::clone(&calls);

        let early_err = (|| -> io::Result<()> {
            let _guard = TuiSessionGuard::from_cleanup(Box::new(move || {
                calls_ref.fetch_add(1, Ordering::SeqCst);
            }));
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "simulated draw failure"))?;
            Ok(())
        })();
        assert!(early_err.is_err());
        assert_eq!(calls.load(Ordering::SeqCst), 1);

        let calls2 = Arc::new(AtomicUsize::new(0));
        let calls2_ref = Arc::clone(&calls2);
        {
            let mut guard = TuiSessionGuard::from_cleanup(Box::new(move || {
                calls2_ref.fetch_add(1, Ordering::SeqCst);
            }));
            guard.release();
            assert_eq!(calls2.load(Ordering::SeqCst), 1);
            guard.release();
        }
        assert_eq!(calls2.load(Ordering::SeqCst), 1);

        let mut buf = Vec::new();
        restore_terminal_buffer(&mut buf).expect("restore terminal sequences");
        let seq = String::from_utf8_lossy(&buf);
        assert!(seq.contains("\x1b[?1049l"));
        assert!(seq.contains("\x1b[?25h"));
    }
}
