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

use ca_core::{EdgeSign, EntityKind};
use egui::{Color32, Context, Pos2, Rect, Shape, Stroke, Vec2, Visuals};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkbenchPaletteMode {
    AcademicDark,
    PublicationLight,
}

pub struct ScientificPalette {
    pub bg_canvas: Color32,
    pub bg_panel: Color32,
    pub grid_line: Color32,
    pub text_primary: Color32,
    pub text_muted: Color32,
    pub activated_strong: Color32,
    pub inhibited_strong: Color32,
    pub expr_up: Color32,
    pub expr_down: Color32,
    pub inconsistent_amber: Color32,
    pub neutral_node: Color32,
    pub border_default: Color32,
}

impl ScientificPalette {
    pub fn for_mode(mode: WorkbenchPaletteMode) -> Self {
        match mode {
            WorkbenchPaletteMode::AcademicDark => Self {
                bg_canvas: Color32::from_rgb(17, 21, 28),
                bg_panel: Color32::from_rgb(24, 29, 38),
                grid_line: Color32::from_rgb(42, 50, 64),
                text_primary: Color32::from_rgb(232, 237, 244),
                text_muted: Color32::from_rgb(148, 160, 178),
                // Colorblind-safe diverging orange (activated) vs steel-blue (inhibited)
                activated_strong: Color32::from_rgb(230, 105, 18),
                inhibited_strong: Color32::from_rgb(43, 131, 186),
                expr_up: Color32::from_rgb(217, 72, 48),
                expr_down: Color32::from_rgb(27, 158, 119),
                inconsistent_amber: Color32::from_rgb(228, 180, 41),
                neutral_node: Color32::from_rgb(54, 63, 79),
                border_default: Color32::from_rgb(112, 126, 148),
            },
            WorkbenchPaletteMode::PublicationLight => Self {
                bg_canvas: Color32::from_rgb(249, 250, 252),
                bg_panel: Color32::from_rgb(240, 243, 247),
                grid_line: Color32::from_rgb(214, 220, 230),
                text_primary: Color32::from_rgb(24, 30, 40),
                text_muted: Color32::from_rgb(94, 106, 124),
                activated_strong: Color32::from_rgb(214, 88, 12),
                inhibited_strong: Color32::from_rgb(33, 113, 181),
                expr_up: Color32::from_rgb(198, 52, 36),
                expr_down: Color32::from_rgb(21, 138, 102),
                inconsistent_amber: Color32::from_rgb(196, 142, 14),
                neutral_node: Color32::from_rgb(225, 230, 238),
                border_default: Color32::from_rgb(118, 130, 148),
            },
        }
    }

    /// Maps an Activation Z-score (typically in `[-4.0, +4.0]`) to a diverging blue-neutral-orange fill.
    pub fn color_for_z(&self, z: f64) -> Color32 {
        let t = (z / 3.5).clamp(-1.0, 1.0) as f32;
        if t >= 0.0 {
            lerp_color(self.neutral_node, self.activated_strong, t)
        } else {
            lerp_color(self.neutral_node, self.inhibited_strong, -t)
        }
    }

    /// Maps normalized activity `[-1.0, +1.0]` (used by MAP propagation) to diverging blue-neutral-orange.
    pub fn color_for_activity(&self, act: f64) -> Color32 {
        let t = act.clamp(-1.0, 1.0) as f32;
        if t >= 0.0 {
            lerp_color(self.neutral_node, self.activated_strong, t)
        } else {
            lerp_color(self.neutral_node, self.inhibited_strong, -t)
        }
    }

    /// Maps observed `log2FC` to coral (upregulated) vs teal (downregulated).
    pub fn color_for_log2fc(&self, fc: f64) -> Color32 {
        let t = (fc / 4.0).clamp(-1.0, 1.0) as f32;
        if t >= 0.0 {
            lerp_color(self.neutral_node, self.expr_up, t)
        } else {
            lerp_color(self.neutral_node, self.expr_down, -t)
        }
    }
}

pub fn apply_workbench_visuals(ctx: &Context, mode: WorkbenchPaletteMode) {
    let pal = ScientificPalette::for_mode(mode);
    let mut visuals = match mode {
        WorkbenchPaletteMode::AcademicDark => Visuals::dark(),
        WorkbenchPaletteMode::PublicationLight => Visuals::light(),
    };
    visuals.panel_fill = pal.bg_panel;
    visuals.window_fill = pal.bg_panel;
    visuals.extreme_bg_color = pal.bg_canvas;
    ctx.set_visuals(visuals);
}

pub fn lerp_color(a: Color32, b: Color32, t: f32) -> Color32 {
    let u = t.clamp(0.0, 1.0);
    let r = (f32::from(a.r()) * (1.0 - u) + f32::from(b.r()) * u).round() as u8;
    let g = (f32::from(a.g()) * (1.0 - u) + f32::from(b.g()) * u).round() as u8;
    let bl = (f32::from(a.b()) * (1.0 - u) + f32::from(b.b()) * u).round() as u8;
    Color32::from_rgb(r, g, bl)
}

/// Renders distinctive biological entity glyphs (Kinase/Enzyme diamond, TF oval, Receptor notched shield,
/// Chemical hexagon, Cytokine/Gene rounded box) centered at `center`.
pub fn draw_entity_glyph(
    painter: &egui::Painter,
    center: Pos2,
    half_extents: Vec2,
    kind: EntityKind,
    fill: Color32,
    stroke: Stroke,
) {
    let hw = half_extents.x;
    let hh = half_extents.y;
    match kind {
        EntityKind::Kinase => {
            // Classic inverted kinase triangle with double top bar
            let pts = vec![
                Pos2::new(center.x - hw * 0.98, center.y - hh * 1.08),
                Pos2::new(center.x + hw * 0.98, center.y - hh * 1.08),
                Pos2::new(center.x, center.y + hh * 1.25),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
            painter.line_segment(
                [
                    Pos2::new(center.x - hw * 0.78, center.y - hh * 0.75),
                    Pos2::new(center.x + hw * 0.78, center.y - hh * 0.75),
                ],
                stroke,
            );
        }
        EntityKind::Phosphatase => {
            // Upward phosphatase triangle
            let pts = vec![
                Pos2::new(center.x, center.y - hh * 1.25),
                Pos2::new(center.x + hw * 0.98, center.y + hh * 1.05),
                Pos2::new(center.x - hw * 0.98, center.y + hh * 1.05),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        EntityKind::Enzyme => {
            // Clean vertical diamond
            let pts = vec![
                Pos2::new(center.x, center.y - hh * 1.2),
                Pos2::new(center.x + hw * 0.95, center.y),
                Pos2::new(center.x, center.y + hh * 1.2),
                Pos2::new(center.x - hw * 0.95, center.y),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        EntityKind::Cytokine | EntityKind::GrowthFactor => {
            // Clean square / compact ligand box
            let sz = (hw * 1.4).min(hh * 2.05);
            let rect = Rect::from_center_size(center, Vec2::new(sz * 1.15, sz));
            painter.rect_filled(rect, 1.5, fill);
            painter.rect_stroke(rect, 1.5, stroke);
        }
        EntityKind::TransmembraneReceptor | EntityKind::IonChannel => {
            // Vertical receptor capsule straddling the plasma membrane bar
            let rect = Rect::from_center_size(center, Vec2::new(hw * 1.38, hh * 2.75));
            painter.rect_filled(rect, hw * 0.69, fill);
            painter.rect_stroke(rect, hw * 0.69, stroke);
        }
        EntityKind::ChemicalDrug => {
            let pts = vec![
                Pos2::new(center.x - hw * 0.6, center.y - hh),
                Pos2::new(center.x + hw * 0.6, center.y - hh),
                Pos2::new(center.x + hw * 1.08, center.y),
                Pos2::new(center.x + hw * 0.6, center.y + hh),
                Pos2::new(center.x - hw * 0.6, center.y + hh),
                Pos2::new(center.x - hw * 1.08, center.y),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        EntityKind::TranscriptionRegulator => {
            let rect = Rect::from_center_size(center, Vec2::new(hw * 2.0, hh * 1.82));
            painter.rect_filled(rect, hh * 0.91, fill);
            painter.rect_stroke(rect, hh * 0.91, stroke);
        }
        EntityKind::BiologicalProcess | EntityKind::DiseasePhenotype => {
            let c = 0.45;
            let pts = vec![
                Pos2::new(center.x - hw * c, center.y - hh),
                Pos2::new(center.x + hw * c, center.y - hh),
                Pos2::new(center.x + hw, center.y - hh * c),
                Pos2::new(center.x + hw, center.y + hh * c),
                Pos2::new(center.x + hw * c, center.y + hh),
                Pos2::new(center.x - hw * c, center.y + hh),
                Pos2::new(center.x - hw, center.y + hh * c),
                Pos2::new(center.x - hw, center.y - hh * c),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        _ => {
            let rect = Rect::from_center_size(center, Vec2::new(hw * 1.8, hh * 1.75));
            painter.rect_filled(rect, 5.0, fill);
            painter.rect_stroke(rect, 5.0, stroke);
        }
    }
}

/// Draws a small circular `(P)` post-translational phosphorylation badge on the top-right rim of a molecule.
pub fn draw_phospho_badge(painter: &egui::Painter, pos: Pos2, bg: Color32, fg: Color32) {
    painter.circle_filled(pos, 6.2, bg);
    painter.circle_stroke(pos, 6.2, Stroke::new(1.1, fg));
    painter.text(
        pos,
        egui::Align2::CENTER_CENTER,
        "P",
        egui::FontId::proportional(8.5),
        fg,
    );
}

/// Draws a promoter double-helix DNA strand across `[x_start, x_end]` centered at `y_center`.
pub fn draw_dna_double_helix(
    painter: &egui::Painter,
    x_start: f32,
    x_end: f32,
    y_center: f32,
    amp: f32,
    color: Color32,
) {
    let span = (x_end - x_start).max(40.0);
    let steps = 72usize;
    let freq = std::f32::consts::PI * 5.0 / span;
    let stroke = Stroke::new(1.15, color);

    let mut gx = x_start + 10.0;
    while gx < x_end - 10.0 {
        let ph = (gx - x_start) * freq;
        let dy = ph.sin() * amp;
        painter.line_segment([Pos2::new(gx, y_center - dy), Pos2::new(gx, y_center + dy)], stroke);
        gx += 14.0;
    }

    for sign in [-1.0f32, 1.0f32] {
        let mut prev = Pos2::new(x_start, y_center);
        for i in 1..=steps {
            let px = x_start + (i as f32 / steps as f32) * span;
            let py = y_center + sign * ((px - x_start) * freq).sin() * amp;
            let cur = Pos2::new(px, py);
            painter.line_segment([prev, cur], stroke);
            prev = cur;
        }
    }
}

/// Draws a directed causal edge between `p1` and `p2` with pointed arrowhead (`Activating`),
/// perpendicular T-bar blunt (`Inhibiting`), or small circle (`Ambiguous`).
pub fn draw_directed_causal_edge(
    painter: &egui::Painter,
    p1: Pos2,
    p2: Pos2,
    target_backoff: f32,
    sign: EdgeSign,
    stroke: Stroke,
) {
    draw_bezier_causal_edge(
        painter,
        p1,
        p2,
        target_backoff * 0.75,
        target_backoff,
        14.0,
        sign,
        stroke,
    );
}

/// Draws a smooth cubic Bezier directed causal edge between `p1` and `p2` with tangent-aligned
/// arrowhead (`Activating`), T-bar blunt (`Inhibiting`), or circle (`Ambiguous`).
#[allow(clippy::too_many_arguments)]
pub fn draw_bezier_causal_edge(
    painter: &egui::Painter,
    p1: Pos2,
    p2: Pos2,
    source_backoff: f32,
    target_backoff: f32,
    curve_offset: f32,
    sign: EdgeSign,
    stroke: Stroke,
) {
    let delta = p2 - p1;
    let dist = delta.length();
    if dist < source_backoff + target_backoff + 4.0 {
        return;
    }
    let dir = delta / dist;
    let perp = Vec2::new(-dir.y, dir.x);
    let start = p1 + dir * source_backoff;
    let tip = p2 - dir * target_backoff;
    let chord = tip - start;
    let chord_len = chord.length().max(1.0);

    let (c1, c2) = if chord.y.abs() >= chord.x.abs() * 0.45 {
        let bow = perp * (curve_offset * 0.65);
        (
            Pos2::new(start.x + bow.x, start.y + chord.y * 0.44 + bow.y),
            Pos2::new(tip.x + bow.x, tip.y - chord.y * 0.44 + bow.y),
        )
    } else {
        let arch = if curve_offset.abs() < 1.0 {
            (chord_len * 0.16).clamp(12.0, 48.0)
        } else {
            curve_offset
        };
        (
            start + dir * (chord_len * 0.28) + perp * arch,
            tip - dir * (chord_len * 0.28) + perp * arch,
        )
    };

    const STEPS: usize = 18;
    let mut pts = Vec::with_capacity(STEPS + 1);
    for step in 0..=STEPS {
        let t = step as f32 / STEPS as f32;
        let u = 1.0 - t;
        let b0 = u * u * u;
        let b1 = 3.0 * u * u * t;
        let b2 = 3.0 * u * t * t;
        let b3 = t * t * t;
        let px = b0 * start.x + b1 * c1.x + b2 * c2.x + b3 * tip.x;
        let py = b0 * start.y + b1 * c1.y + b2 * c2.y + b3 * tip.y;
        pts.push(Pos2::new(px, py));
    }

    let end_delta = pts[STEPS] - pts[STEPS - 2];
    let end_dir = if end_delta.length_sq() > 1e-4 {
        end_delta.normalized()
    } else {
        dir
    };
    let end_perp = Vec2::new(-end_dir.y, end_dir.x);

    painter.add(Shape::line(pts, stroke));

    match sign {
        EdgeSign::Activating => {
            let arrow_len = 9.0f32;
            let arrow_half_w = 4.6f32;
            let base = tip - end_dir * arrow_len;
            let tri = vec![
                tip,
                base + end_perp * arrow_half_w,
                base - end_perp * arrow_half_w,
            ];
            painter.add(Shape::convex_polygon(tri, stroke.color, Stroke::NONE));
        }
        EdgeSign::Inhibiting => {
            let bar_half = 6.5f32;
            let t1 = tip + end_perp * bar_half;
            let t2 = tip - end_perp * bar_half;
            painter.line_segment([t1, t2], Stroke::new(stroke.width + 1.0, stroke.color));
        }
        EdgeSign::Ambiguous => {
            painter.circle_filled(tip, 3.5, stroke.color);
        }
    }
}

/// Draws a bright pulsing gold selection/search halo around a matched molecule node.
pub fn draw_gold_pulse_halo(
    painter: &egui::Painter,
    center: Pos2,
    half_extents: Vec2,
    phase: f32,
) {
    let pulse = 0.5 + 0.5 * (phase * 6.0).sin();
    let outer_pad = 9.0 + pulse * 5.0;
    let inner_pad = 5.0 + pulse * 2.0;

    let outer_rect = Rect::from_center_size(
        center,
        Vec2::new(
            half_extents.x * 2.0 + outer_pad * 2.0,
            half_extents.y * 2.0 + outer_pad * 2.0,
        ),
    );
    painter.rect_stroke(
        outer_rect,
        10.0,
        Stroke::new(2.2_f32, Color32::from_rgba_unmultiplied(250, 204, 21, 145)),
    );

    let inner_rect = Rect::from_center_size(
        center,
        Vec2::new(
            half_extents.x * 2.0 + inner_pad * 2.0,
            half_extents.y * 2.0 + inner_pad * 2.0,
        ),
    );
    painter.rect_stroke(
        inner_rect,
        8.0,
        Stroke::new(3.0_f32, Color32::from_rgb(253, 224, 71)),
    );
}
