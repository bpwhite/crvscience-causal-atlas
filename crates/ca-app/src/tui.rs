#![forbid(unsafe_code)]

use ca_core::CoreAnalysisReport;
use crossterm::{
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

/// Launches the headless SSH-compatible terminal workbench (`causal-atlas tui`).
pub fn run_headless_tui(report: &CoreAnalysisReport) -> io::Result<()> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
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
                            .map(|z| format!("{z:+.2}"))
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
                            Cell::from(format!("{:+.2}", r.activation_z)),
                            Cell::from(format!("{:+.2}", r.bias_corrected_z)),
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
                            Cell::from(format!("{:+.2}", c.activation_z)),
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

    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;
    loop_res
}
