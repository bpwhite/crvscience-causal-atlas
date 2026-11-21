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
    painter.circle_stroke(pos, 6.2, Stroke::new(1.1_f32, fg));
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
    let stroke = Stroke::new(1.15_f32, color);

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
#[allow(dead_code)]
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

/// Screen-space node bounding obstacle used to deflect Bezier edges around intermediate nodes.
#[derive(Debug, Clone)]
pub struct EdgeObstacle {
    pub id: String,
    pub center: Pos2,
    pub half_extents: Vec2,
}

fn baseline_bezier_controls(
    start: Pos2,
    tip: Pos2,
    dir: Vec2,
    perp: Vec2,
    curve_offset: f32,
) -> (Pos2, Pos2) {
    let chord = tip - start;
    let chord_len = chord.length().max(1.0);
    if chord.y.abs() >= chord.x.abs() * 0.45 {
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
    }
}

fn sample_cubic_bezier(start: Pos2, c1: Pos2, c2: Pos2, tip: Pos2, steps: usize) -> Vec<Pos2> {
    let mut pts = Vec::with_capacity(steps + 1);
    for step in 0..=steps {
        let t = step as f32 / steps as f32;
        let u = 1.0 - t;
        let b0 = u * u * u;
        let b1 = 3.0 * u * u * t;
        let b2 = 3.0 * u * t * t;
        let b3 = t * t * t;
        let px = b0 * start.x + b1 * c1.x + b2 * c2.x + b3 * tip.x;
        let py = b0 * start.y + b1 * c1.y + b2 * c2.y + b3 * tip.y;
        pts.push(Pos2::new(px, py));
    }
    pts
}

fn segment_ellipse_closest(a: Pos2, b: Pos2, o: &EdgeObstacle) -> (f32, f32, Pos2) {
    let rx = o.half_extents.x.max(1.0);
    let ry = o.half_extents.y.max(1.0);
    let ax = (a.x - o.center.x) / rx;
    let ay = (a.y - o.center.y) / ry;
    let bx = (b.x - o.center.x) / rx;
    let by = (b.y - o.center.y) / ry;
    let dx = bx - ax;
    let dy = by - ay;
    let denom = dx * dx + dy * dy;
    let s = if denom > 1e-8 {
        (-(ax * dx + ay * dy) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let qx = ax + s * dx;
    let qy = ay + s * dy;
    let r = qx.hypot(qy);
    let p = Pos2::new(a.x + (b.x - a.x) * s, a.y + (b.y - a.y) * s);
    (r, s, p)
}

/// Computes cubic Bezier polyline samples between `start` and `tip`, smoothly deflecting
/// control points around any intermediate `obstacles` (`id != from_id && id != to_id`).
pub fn solve_deflected_cubic_bezier(
    start: Pos2,
    tip: Pos2,
    curve_offset: f32,
    obstacles: &[EdgeObstacle],
    from_id: &str,
    to_id: &str,
) -> Vec<Pos2> {
    let chord = tip - start;
    let chord_len = chord.length().max(1.0);
    let dir = chord / chord_len;
    let perp = Vec2::new(-dir.y, dir.x);
    let (c1_base, c2_base) = baseline_bezier_controls(start, tip, dir, perp, curve_offset);

    const STEPS: usize = 24;
    if obstacles.is_empty() {
        return sample_cubic_bezier(start, c1_base, c2_base, tip, 18);
    }

    let pad = 360.0 + curve_offset.abs();
    let min_x = start.x.min(tip.x) - pad;
    let max_x = start.x.max(tip.x) + pad;
    let min_y = start.y.min(tip.y) - pad;
    let max_y = start.y.max(tip.y) + pad;

    let nearby: Vec<&EdgeObstacle> = obstacles
        .iter()
        .filter(|o| {
            o.id != from_id
                && o.id != to_id
                && o.center.x >= min_x
                && o.center.x <= max_x
                && o.center.y >= min_y
                && o.center.y <= max_y
        })
        .collect();

    if nearby.is_empty() {
        return sample_cubic_bezier(start, c1_base, c2_base, tip, 18);
    }

    // Check every polyline segment `[pts[i], pts[i + 1]]` against every nearby obstacle
    // ellipse so neither near-endpoint nodes (t = 1/24, 23/24) nor segment-straddling chords
    // can slice through intermediate glyphs. Only immovable exact endpoints (`start`/`tip`)
    // already overlapping an obstacle are ignored if the segment exits immediately.
    let eval_hits = |pts: &[Pos2]| -> (usize, f32, Option<(f32, usize, Pos2)>) {
        let mut count = 0usize;
        let mut min_r = f32::INFINITY;
        let mut worst: Option<(f32, usize, Pos2)> = None;
        let seg_count = pts.len().saturating_sub(1);
        for i in 0..seg_count {
            let a = pts[i];
            let b = pts[i + 1];
            for (oi, o) in nearby.iter().enumerate() {
                let (r, s, p) = segment_ellipse_closest(a, b, o);
                if (i == 0 && s <= 1e-4) || (i + 1 == seg_count && s >= 1.0 - 1e-4) {
                    continue;
                }
                if r < 1.0 {
                    count += 1;
                }
                if r < min_r {
                    min_r = r;
                    if r < 1.0 {
                        let t_hit = ((i as f32 + s) / seg_count as f32).clamp(0.02, 0.98);
                        worst = Some((t_hit, oi, p));
                    }
                }
            }
        }
        (count, min_r, worst)
    };

    let base_pts = sample_cubic_bezier(start, c1_base, c2_base, tip, STEPS);
    let (base_hits, _, base_worst) = eval_hits(&base_pts);
    if base_hits == 0 {
        return base_pts;
    }

    let preferred_side = if let Some((_, oi, _)) = base_worst {
        let lat = (nearby[oi].center - start).dot(perp);
        if lat.abs() >= 1e-2 {
            if lat >= 0.0 {
                -1.0f32
            } else {
                1.0f32
            }
        } else if curve_offset >= 0.0 {
            1.0f32
        } else {
            -1.0f32
        }
    } else if curve_offset >= 0.0 {
        1.0f32
    } else {
        -1.0f32
    };

    let solve_side = |side: f32| -> (Vec<Pos2>, usize, f32, f32) {
        let mut c1 = c1_base;
        let mut c2 = c2_base;
        let push_dir = perp * side;
        let mut best: Option<(Vec<Pos2>, usize, f32, f32)> = None;
        for _ in 0..14 {
            let pts = sample_cubic_bezier(start, c1, c2, tip, STEPS);
            let (hits, min_r, worst) = eval_hits(&pts);
            let exc = (c1 - c1_base).length() + (c2 - c2_base).length();
            if hits == 0 {
                return (pts, 0, min_r.max(1.0), exc);
            }
            let replace_best = match &best {
                None => true,
                Some(b) => {
                    hits < b.1
                        || (hits == b.1 && (min_r > b.2 || (min_r == b.2 && exc < b.3)))
                }
            };
            if replace_best {
                best = Some((pts, hits, min_r, exc));
            }
            let Some((t_hit, oi, p_hit)) = worst else {
                break;
            };
            let o = nearby[oi];
            let rx = o.half_extents.x.max(1.0);
            let ry = o.half_extents.y.max(1.0);
            let eff_r = 1.0 / ((push_dir.x / rx).hypot(push_dir.y / ry)).max(1e-4);
            let proj = (p_hit - o.center).dot(push_dir);
            let shift = (eff_r * 1.15 - proj).clamp(12.0, 180.0);
            let u = 1.0 - t_hit;
            let denom = (3.0 * u * u * t_hit * u + 3.0 * u * t_hit * t_hit * t_hit).max(0.04);
            c1 += push_dir * (shift * u / denom).clamp(0.0, 520.0);
            c2 += push_dir * (shift * t_hit / denom).clamp(0.0, 520.0);
        }
        let pts = sample_cubic_bezier(start, c1, c2, tip, STEPS);
        let (hits, min_r, _) = eval_hits(&pts);
        let exc = (c1 - c1_base).length() + (c2 - c2_base).length();
        match best {
            Some(b)
                if b.1 < hits
                    || (b.1 == hits && (b.2 > min_r || (b.2 == min_r && b.3 < exc))) =>
            {
                b
            }
            _ => (pts, hits, min_r, exc),
        }
    };

    let cand_a = solve_side(preferred_side);
    if cand_a.1 == 0 {
        return cand_a.0;
    }
    let cand_b = solve_side(-preferred_side);
    if cand_b.1 < cand_a.1
        || (cand_b.1 == cand_a.1
            && (cand_b.2 > cand_a.2 || (cand_b.2 == cand_a.2 && cand_b.3 < cand_a.3)))
    {
        cand_b.0
    } else {
        cand_a.0
    }
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
    draw_deflected_bezier_causal_edge(
        painter,
        p1,
        p2,
        source_backoff,
        target_backoff,
        curve_offset,
        sign,
        stroke,
        &[],
        "",
        "",
    );
}

/// Draws a smooth cubic Bezier directed causal edge deflected around intermediate `obstacles`.
#[allow(clippy::too_many_arguments)]
pub fn draw_deflected_bezier_causal_edge(
    painter: &egui::Painter,
    p1: Pos2,
    p2: Pos2,
    source_backoff: f32,
    target_backoff: f32,
    curve_offset: f32,
    sign: EdgeSign,
    stroke: Stroke,
    obstacles: &[EdgeObstacle],
    from_id: &str,
    to_id: &str,
) {
    let delta = p2 - p1;
    let dist = delta.length();
    if dist < source_backoff + target_backoff + 4.0 {
        return;
    }
    let dir = delta / dist;
    let start = p1 + dir * source_backoff;
    let tip = p2 - dir * target_backoff;

    let pts = solve_deflected_cubic_bezier(start, tip, curve_offset, obstacles, from_id, to_id);
    let end_delta = if pts.len() >= 3 {
        pts[pts.len() - 1] - pts[pts.len() - 3]
    } else {
        delta
    };
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deflected_cubic_bezier_routes_cleanly_around_intermediate_node() {
        // Collinear vertical cascade: TNF -> VCAM1 (intermediate) -> KLF2
        let start = Pos2::new(520.0, 80.0);
        let tip = Pos2::new(560.0, 560.0);
        let vcam1 = EdgeObstacle {
            id: "VCAM1".to_string(),
            center: Pos2::new(525.0, 175.0),
            half_extents: Vec2::new(38.0, 28.0),
        };
        let pts = solve_deflected_cubic_bezier(
            start,
            tip,
            12.0,
            std::slice::from_ref(&vcam1),
            "TNF",
            "KLF2",
        );
        assert!(pts.len() >= 19);
        for p in &pts[2..(pts.len() - 2)] {
            let r = ((p.x - vcam1.center.x) / vcam1.half_extents.x)
                .hypot((p.y - vcam1.center.y) / vcam1.half_extents.y);
            assert!(
                r >= 0.99,
                "deflected edge sample {p:?} penetrated VCAM1 obstacle (r={r:.3})"
            );
        }
    }

    #[test]
    fn deflected_cubic_bezier_clears_collinear_multi_node_column() {
        let start = Pos2::new(480.0, 70.0);
        let tip = Pos2::new(480.0, 580.0);
        let obstacles = vec![
            EdgeObstacle {
                id: "TNFRSF1A".to_string(),
                center: Pos2::new(480.0, 180.0),
                half_extents: Vec2::new(36.0, 26.0),
            },
            EdgeObstacle {
                id: "TRAF2".to_string(),
                center: Pos2::new(480.0, 310.0),
                half_extents: Vec2::new(36.0, 22.0),
            },
            EdgeObstacle {
                id: "MAP3K7".to_string(),
                center: Pos2::new(480.0, 440.0),
                half_extents: Vec2::new(36.0, 24.0),
            },
        ];
        let pts = solve_deflected_cubic_bezier(start, tip, 12.0, &obstacles, "TNF", "IKBKB");
        for p in &pts[2..(pts.len() - 2)] {
            for obs in &obstacles {
                let r = ((p.x - obs.center.x) / obs.half_extents.x)
                    .hypot((p.y - obs.center.y) / obs.half_extents.y);
                assert!(
                    r >= 0.98,
                    "sample {p:?} hit obstacle {} (r={r:.3})",
                    obs.id
                );
            }
        }
    }

    #[test]
    fn deflected_cubic_bezier_clears_near_start_adv_a_and_segment_straddle() {
        // Synthetic Adv A: intermediate obstacle right at t = 1/24 (pt[1])
        let start = Pos2::new(400.0, 100.0);
        let tip = Pos2::new(400.0, 700.0);
        let t_near = 1.0_f32 / 24.0;
        let adv_center = Pos2::new(400.0, start.y + (tip.y - start.y) * t_near);
        let adv_a = EdgeObstacle {
            id: "ADV_A".to_string(),
            center: adv_center,
            half_extents: Vec2::new(36.0, 22.0),
        };
        let pts = solve_deflected_cubic_bezier(
            start,
            tip,
            12.0,
            std::slice::from_ref(&adv_a),
            "SRC",
            "DST",
        );
        for seg in pts.windows(2) {
            let (r, _, _) = segment_ellipse_closest(seg[0], seg[1], &adv_a);
            assert!(
                r >= 0.99,
                "Adv A obstacle penetrated along segment {:?} -> {:?} (r={r:.3})",
                seg[0],
                seg[1]
            );
        }

        // Segment-straddling probe: obstacle sitting between two chord samples on deflected path
        let mid_obs = EdgeObstacle {
            id: "STRADDLE".to_string(),
            center: Pos2::new(440.0, 220.0),
            half_extents: Vec2::new(35.0, 22.0),
        };
        let pts2 = solve_deflected_cubic_bezier(
            Pos2::new(420.0, 140.0),
            Pos2::new(460.0, 420.0),
            12.0,
            std::slice::from_ref(&mid_obs),
            "A",
            "B",
        );
        for seg in pts2.windows(2) {
            for k in 0..=16 {
                let s = k as f32 / 16.0;
                let p = Pos2::new(
                    seg[0].x + (seg[1].x - seg[0].x) * s,
                    seg[0].y + (seg[1].y - seg[0].y) * s,
                );
                let r = ((p.x - mid_obs.center.x) / mid_obs.half_extents.x)
                    .hypot((p.y - mid_obs.center.y) / mid_obs.half_extents.y);
                assert!(
                    r >= 0.99,
                    "segment-straddling obstacle hit at {p:?} (r={r:.3})"
                );
            }
        }
    }

    #[test]
    fn deflected_cubic_bezier_clears_all_seed_pathways_at_all_zooms() {
        use ca_core::{
            compute_pathway_world_layout_cached, seed_kb::build_seed_knowledge_graph, EntityKind,
        };
        use std::collections::{HashMap, HashSet};

        let kg = build_seed_knowledge_graph();
        let canvas_rect = Rect::from_min_size(Pos2::new(340.0, 120.0), Vec2::new(1240.0, 840.0));

        for pw in &kg.pathways {
            let wl = compute_pathway_world_layout_cached(pw, &kg);
            let fit_zoom = ((canvas_rect.width() - 36.0) / wl.world_w)
                .min((canvas_rect.height() - 36.0) / wl.world_h)
                .clamp(0.15, 3.5);
            let world_center = Pos2::new(wl.world_w * 0.5, wl.world_h * 0.5);

            for zoom in [0.5_f32, fit_zoom, 1.0_f32, 1.5_f32] {
                let to_screen = |wx: f32, wy: f32| -> Pos2 {
                    Pos2::new(
                        canvas_rect.center().x + (wx - world_center.x) * zoom,
                        canvas_rect.center().y + (wy - world_center.y) * zoom,
                    )
                };

                let base_half = Vec2::new((32.0 * zoom).max(16.0), (16.0 * zoom).max(9.0));
                let obs_pad = (7.0 * zoom).clamp(4.0, 10.0);
                let mut obstacles = Vec::with_capacity(wl.positions.len());
                let mut raw_extents: HashMap<String, (Pos2, Vec2)> = HashMap::new();
                for (nid, &(wx, wy)) in &wl.positions {
                    let kind = kg
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
                    let center = to_screen(wx, wy);
                    let raw_half = Vec2::new(base_half.x * scale.x, base_half.y * scale.y);
                    raw_extents.insert(nid.clone(), (center, raw_half));
                    obstacles.push(EdgeObstacle {
                        id: nid.clone(),
                        center,
                        half_extents: Vec2::new(raw_half.x + obs_pad, raw_half.y + obs_pad),
                    });
                }

                let mut pair_counts: HashMap<(String, String), usize> = HashMap::new();
                let mut edge_seen_rev: HashSet<(String, String)> = HashSet::new();
                for e in &wl.edges {
                    edge_seen_rev.insert((e.from.clone(), e.to.clone()));
                }

                let source_backoff = (20.0 * zoom).clamp(8.0, 36.0);
                let target_backoff = (24.0 * zoom).clamp(10.0, 42.0);

                for edge in &wl.edges {
                    let Some(&(wx1, wy1)) = wl.positions.get(&edge.from) else {
                        continue;
                    };
                    let Some(&(wx2, wy2)) = wl.positions.get(&edge.to) else {
                        continue;
                    };
                    let p1 = to_screen(wx1, wy1);
                    let p2 = to_screen(wx2, wy2);
                    let delta = p2 - p1;
                    let dist = delta.length();
                    if dist < source_backoff + target_backoff + 4.0 {
                        continue;
                    }
                    let dir = delta / dist;
                    let start = p1 + dir * source_backoff;
                    let tip = p2 - dir * target_backoff;

                    let key = (edge.from.clone(), edge.to.clone());
                    let dup_i = pair_counts.entry(key).or_insert(0);
                    let has_rev = edge_seen_rev.contains(&(edge.to.clone(), edge.from.clone()));
                    let bow = if has_rev {
                        22.0 * zoom
                    } else {
                        (12.0 + (*dup_i as f32) * 14.0) * zoom
                    };
                    *dup_i += 1;

                    let pts = solve_deflected_cubic_bezier(
                        start,
                        tip,
                        bow,
                        &obstacles,
                        &edge.from,
                        &edge.to,
                    );

                    for (seg_idx, seg) in pts.windows(2).enumerate() {
                        for k in 0..=16 {
                            if (seg_idx == 0 && k == 0)
                                || (seg_idx + 2 == pts.len() && k == 16)
                            {
                                continue;
                            }
                            let s = k as f32 / 16.0;
                            let p = Pos2::new(
                                seg[0].x + (seg[1].x - seg[0].x) * s,
                                seg[0].y + (seg[1].y - seg[0].y) * s,
                            );
                            for (oid, &(oc, oh)) in &raw_extents {
                                if oid == &edge.from || oid == &edge.to {
                                    continue;
                                }
                                let r_raw = ((p.x - oc.x) / oh.x).hypot((p.y - oc.y) / oh.y);
                                assert!(
                                    r_raw >= 1.0,
                                    "Pathway {} edge {} -> {} at zoom={zoom:.3} penetrated raw glyph {} (r_raw={r_raw:.3}, p={p:?})",
                                    pw.id,
                                    edge.from,
                                    edge.to,
                                    oid
                                );
                            }
                            for obs in &obstacles {
                                if obs.id == edge.from || obs.id == edge.to {
                                    continue;
                                }
                                let r_pad = ((p.x - obs.center.x) / obs.half_extents.x)
                                    .hypot((p.y - obs.center.y) / obs.half_extents.y);
                                assert!(
                                    r_pad >= 0.995,
                                    "Pathway {} edge {} -> {} at zoom={zoom:.3} penetrated padded obstacle {} (r_pad={r_pad:.3}, p={p:?})",
                                    pw.id,
                                    edge.from,
                                    edge.to,
                                    obs.id
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
