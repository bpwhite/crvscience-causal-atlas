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

use ca_core::{
    compute_pathway_world_layout_cached, EdgeSign, KnowledgeGraph, MapSimulationResult, PathwayDef,
    SubcellularLocation, TreeMapRect,
};
use std::collections::BTreeMap;
use std::fmt::Write;

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn z_hex_color(v: f64, limit: f64) -> String {
    let t = (v / limit).clamp(-1.0, 1.0);
    // Neutral light slate (226, 232, 240) -> Orange (230, 105, 18) or Blue (43, 131, 186)
    let (r, g, b) = if t >= 0.0 {
        (
            (226.0 * (1.0 - t) + 230.0 * t).round() as u8,
            (232.0 * (1.0 - t) + 105.0 * t).round() as u8,
            (240.0 * (1.0 - t) + 18.0 * t).round() as u8,
        )
    } else {
        let u = -t;
        (
            (226.0 * (1.0 - u) + 43.0 * u).round() as u8,
            (232.0 * (1.0 - u) + 131.0 * u).round() as u8,
            (240.0 * (1.0 - u) + 186.0 * u).round() as u8,
        )
    };
    format!("#{r:02X}{g:02X}{b:02X}")
}

const NODE_WIDTH: f32 = 74.0;
const NODE_HEIGHT: f32 = 30.0;
const NODE_HALF_W: f32 = NODE_WIDTH * 0.5;
const NODE_HALF_H: f32 = NODE_HEIGHT * 0.5;
const EDGE_CLIP_CLEARANCE: f32 = 1.2;

/// Clips a ray from rectangle center `(cx, cy)` toward `(toward_x, toward_y)` onto the
/// perimeter of `[cx ± half_w, cy ± half_h]`.
pub fn clip_ray_to_rect_perimeter(
    cx: f32,
    cy: f32,
    toward_x: f32,
    toward_y: f32,
    half_w: f32,
    half_h: f32,
) -> (f32, f32) {
    let dx = toward_x - cx;
    let dy = toward_y - cy;
    if dx.abs() < 1e-5 && dy.abs() < 1e-5 {
        return (cx, cy);
    }
    let tx = if dx.abs() > 1e-6 {
        half_w / dx.abs()
    } else {
        f32::INFINITY
    };
    let ty = if dy.abs() > 1e-6 {
        half_h / dy.abs()
    } else {
        f32::INFINITY
    };
    let t = tx.min(ty);
    (cx + dx * t, cy + dy * t)
}

/// Computes Quadratic Bezier endpoints `(start, control, end)` clipped outside both
/// `74x30` node rectangles (`[cx±37, cy±15]`) so SVG markers render outside node borders.
pub fn clip_quadratic_edge_to_node_boxes(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
) -> Option<((f32, f32), (f32, f32), (f32, f32))> {
    let chord_len = (x2 - x1).hypot(y2 - y1);
    if chord_len < 4.0 {
        return None;
    }
    let mx = (x1 + x2) * 0.5 + (y2 - y1) * 0.08;
    let my = (y1 + y2) * 0.5 - (x2 - x1) * 0.08;
    let hw = NODE_HALF_W + EDGE_CLIP_CLEARANCE;
    let hh = NODE_HALF_H + EDGE_CLIP_CLEARANCE;
    let (sx, sy) = clip_ray_to_rect_perimeter(x1, y1, mx, my, hw, hh);
    let (ex, ey) = clip_ray_to_rect_perimeter(x2, y2, mx, my, hw, hh);
    if (ex - sx).hypot(ey - sy) < 2.0 {
        return None;
    }
    Some(((sx, sy), (mx, my), (ex, ey)))
}

/// Renders a Canonical Pathway subcellular topology + Molecule Activity Predictor (MAP) state to vector SVG.
pub fn render_pathway_map_svg(
    pathway: &PathwayDef,
    map_res: &MapSimulationResult,
    kg: &KnowledgeGraph,
) -> String {
    let wl = compute_pathway_world_layout_cached(pathway, kg);
    let width = wl.world_w.max(1100.0);
    let plot_h = wl.world_h.max(700.0);
    let legend_h = 48.0f32;
    let height = plot_h + legend_h;
    let mut svg = String::with_capacity(16_384);

    let _ = writeln!(
        &mut svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" font-family="Inter, Helvetica, Arial, sans-serif">"##
    );
    let _ = writeln!(
        &mut svg,
        r##"<defs>
  <marker id="arrow-act" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 1 L 10 5 L 0 9 z" fill="#D95F02"/>
  </marker>
  <marker id="arrow-inh" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto">
    <line x1="8" y1="1" x2="8" y2="9" stroke="#2B83BA" stroke-width="2.5"/>
  </marker>
  <linearGradient id="map-activity-grad" x1="0%" y1="0%" x2="100%" y2="0%">
    <stop offset="0%" stop-color="#2B83BA"/>
    <stop offset="50%" stop-color="#E2E8F0"/>
    <stop offset="100%" stop-color="#E66912"/>
  </linearGradient>
</defs>"##
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="0" y="0" width="{width}" height="{height}" fill="#F8FAFC"/>"##
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="28" y="34" font-size="18" font-weight="700" fill="#0F172A">{} ({}) — Subcellular MAP Diagram</text>"##,
        xml_escape(&pathway.name),
        xml_escape(&pathway.id)
    );

    let scale_x = width / wl.world_w.max(1.0);
    let scale_y = (plot_h - 40.0) / wl.world_h.max(1.0);
    let pm_y = 40.0 + wl.bands.plasma_membrane_y * scale_y;
    let org_y = 40.0 + wl.bands.organelle_divider_y * scale_y;
    let nm_y = 40.0 + wl.bands.nuclear_envelope_y * scale_y;

    // Subtle subcellular compartment band fills
    let _ = writeln!(
        &mut svg,
        r##"<rect x="28" y="44.0" width="{:.1}" height="{:.1}" fill="#EFF6FF" opacity="0.65"/>"##,
        width - 56.0,
        (pm_y - 44.0).max(0.0)
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="28" y="{:.1}" width="{:.1}" height="{:.1}" fill="#F1F5F9" opacity="0.72"/>"##,
        pm_y,
        width - 56.0,
        (org_y - pm_y).max(0.0)
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="48" y="{:.1}" width="{:.1}" height="{:.1}" fill="#EEF2FF" opacity="0.60"/>"##,
        org_y,
        width - 96.0,
        (nm_y - org_y).max(0.0)
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="60" y="{:.1}" width="{:.1}" height="{:.1}" fill="#F5F3FF" opacity="0.65"/>"##,
        nm_y,
        width - 120.0,
        (plot_h - nm_y - 4.0).max(0.0)
    );

    let _ = writeln!(
        &mut svg,
        r##"<rect x="36" y="{:.1}" width="{:.1}" height="6" fill="#73777D"/>"##,
        pm_y - 3.0,
        width - 72.0
    );
    let _ = writeln!(
        &mut svg,
        r##"<line x1="56" y1="{org_y:.1}" x2="{:.1}" y2="{org_y:.1}" stroke="#CBD5E1" stroke-width="1.5" stroke-dasharray="6,4"/>"##,
        width - 56.0
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="68" y="{:.1}" width="{:.1}" height="6" fill="#73777D"/>"##,
        nm_y - 3.0,
        width - 136.0
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="38" y="{:.1}" font-size="12" font-weight="600" fill="#374151">{}</text>"##,
        pm_y - 10.0,
        SubcellularLocation::Extracellular.label()
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="38" y="{:.1}" font-size="12" font-weight="600" fill="#374151">{}</text>"##,
        pm_y + 20.0,
        SubcellularLocation::Cytoplasm.label()
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="58" y="{:.1}" font-size="11" font-weight="600" fill="#64748B">{}</text>"##,
        org_y + 15.0,
        SubcellularLocation::Organelle.label()
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="70" y="{:.1}" font-size="12" font-weight="600" fill="#374151">{}</text>"##,
        nm_y + 20.0,
        SubcellularLocation::Nucleus.label()
    );

    let coords: BTreeMap<String, (f32, f32)> = wl
        .positions
        .into_iter()
        .map(|(nid, (wx, wy))| (nid, (wx * scale_x, 40.0 + wy * scale_y)))
        .collect();

    for edge in &wl.edges {
        if let (Some(&(x1, y1)), Some(&(x2, y2))) = (coords.get(&edge.from), coords.get(&edge.to)) {
            let Some(((sx, sy), (mx, my), (ex, ey))) =
                clip_quadratic_edge_to_node_boxes(x1, y1, x2, y2)
            else {
                continue;
            };
            let (stroke, marker) = match edge.sign {
                EdgeSign::Activating => ("#D95F02", "url(#arrow-act)"),
                EdgeSign::Inhibiting => ("#2B83BA", "url(#arrow-inh)"),
                EdgeSign::Ambiguous => ("#64748B", ""),
            };
            let marker_attr = if marker.is_empty() {
                String::new()
            } else {
                format!(r#" marker-end="{marker}""#)
            };
            let _ = writeln!(
                &mut svg,
                r##"<path d="M {sx:.1} {sy:.1} Q {mx:.1} {my:.1} {ex:.1} {ey:.1}" fill="none" stroke="{stroke}" stroke-width="2.0" opacity="0.85"{marker_attr}/>"##
            );
        }
    }

    for (nid, st) in &map_res.states {
        if let Some(&(cx, cy)) = coords.get(nid) {
            let fill = z_hex_color(st.predicted_activity, 1.0);
            let border = if st.user_clamp.is_some() {
                "#7C3AED"
            } else if st.observed_log2fc.is_some() {
                "#0F172A"
            } else {
                "#64748B"
            };
            let bw = if st.observed_log2fc.is_some() || st.user_clamp.is_some() {
                2.4
            } else {
                1.2
            };
            let _ = writeln!(
                &mut svg,
                r##"<rect x="{:.1}" y="{:.1}" width="{NODE_WIDTH:.0}" height="{NODE_HEIGHT:.0}" rx="6" fill="{fill}" stroke="{border}" stroke-width="{bw}"/>"##,
                cx - NODE_HALF_W,
                cy - NODE_HALF_H
            );
            let _ = writeln!(
                &mut svg,
                r##"<text x="{cx:.1}" y="{:.1}" text-anchor="middle" font-size="11" font-weight="600" fill="#0F172A">{}</text>"##,
                cy + 4.0,
                xml_escape(nid)
            );
        }
    }

    // Activity Colorbar & Edge-Sign Legend
    let leg_y = plot_h + 6.0;
    let leg_mid_y = leg_y + 18.0;
    let _ = writeln!(
        &mut svg,
        r##"<g id="activity-colorbar-legend">
  <rect x="28" y="{leg_y:.1}" width="{:.1}" height="34" rx="6" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="1"/>
  <text x="42" y="{:.1}" font-size="11" font-weight="700" fill="#0F172A">MAP Activity:</text>
  <text x="128" y="{:.1}" font-size="10" fill="#2B83BA">-1.0 (Inhibited)</text>
  <rect x="212" y="{:.1}" width="170" height="12" rx="3" fill="url(#map-activity-grad)" stroke="#94A3B8" stroke-width="0.8"/>
  <text x="297" y="{:.1}" text-anchor="middle" font-size="9.5" fill="#475569">0.0</text>
  <text x="390" y="{:.1}" font-size="10" fill="#D95F02">+1.0 (Activated)</text>
  <line x1="500" y1="{:.1}" x2="528" y2="{:.1}" stroke="#D95F02" stroke-width="2" marker-end="url(#arrow-act)"/>
  <text x="536" y="{:.1}" font-size="10" fill="#334155">Activating</text>
  <line x1="608" y1="{:.1}" x2="636" y2="{:.1}" stroke="#2B83BA" stroke-width="2" marker-end="url(#arrow-inh)"/>
  <text x="644" y="{:.1}" font-size="10" fill="#334155">Inhibiting</text>
  <rect x="716" y="{:.1}" width="14" height="10" rx="2" fill="#E2E8F0" stroke="#0F172A" stroke-width="2"/>
  <text x="736" y="{:.1}" font-size="10" fill="#334155">Measured DEG</text>
  <rect x="828" y="{:.1}" width="14" height="10" rx="2" fill="#E2E8F0" stroke="#7C3AED" stroke-width="2"/>
  <text x="848" y="{:.1}" font-size="10" fill="#334155">MAP Clamp</text>
</g>"##,
        width - 56.0,
        leg_mid_y + 3.0,
        leg_mid_y + 3.0,
        leg_y + 11.0,
        leg_y + 9.0,
        leg_mid_y + 3.0,
        leg_mid_y,
        leg_mid_y,
        leg_mid_y + 3.0,
        leg_mid_y,
        leg_mid_y,
        leg_mid_y + 3.0,
        leg_y + 12.0,
        leg_mid_y + 3.0,
        leg_y + 12.0,
        leg_mid_y + 3.0,
    );

    svg.push_str("</svg>\n");
    svg
}

/// Renders the Downstream Effects hierarchical TreeMap to vector SVG.
pub fn render_treemap_svg(dataset_name: &str, rects: &[TreeMapRect]) -> String {
    let width = 1000.0f32;
    let height = 640.0f32;
    let pad_top = 52.0f32;
    let plot_w = width - 32.0;
    let plot_h = height - pad_top - 20.0;

    let mut svg = String::with_capacity(12_288);
    let _ = writeln!(
        &mut svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" font-family="Inter, Helvetica, Arial, sans-serif">"##
    );
    let _ = writeln!(
        &mut svg,
        r##"<rect x="0" y="0" width="{width}" height="{height}" fill="#F8FAFC"/>"##
    );
    let _ = writeln!(
        &mut svg,
        r##"<text x="16" y="32" font-size="17" font-weight="700" fill="#0F172A">Downstream Effects TreeMap — {}</text>"##,
        xml_escape(dataset_name)
    );

    for r in rects {
        let rx = 16.0 + r.x * plot_w;
        let ry = pad_top + r.y * plot_h;
        let rw = (r.w * plot_w - 2.0).max(2.0);
        let rh = (r.h * plot_h - 2.0).max(2.0);
        let fill = z_hex_color(r.activation_z, 3.5);
        let _ = writeln!(
            &mut svg,
            r##"<rect x="{rx:.1}" y="{ry:.1}" width="{rw:.1}" height="{rh:.1}" rx="3" fill="{fill}" stroke="#FFFFFF" stroke-width="1.5"/>"##
        );
        if rw > 70.0 && rh > 28.0 {
            let _ = writeln!(
                &mut svg,
                r##"<text x="{:.1}" y="{:.1}" font-size="11" font-weight="600" fill="#0F172A">{}</text>"##,
                rx + 6.0,
                ry + 16.0,
                xml_escape(&r.term_name)
            );
            let _ = writeln!(
                &mut svg,
                r##"<text x="{:.1}" y="{:.1}" font-size="10" fill="#1E293B">z = {} | -log10(p) = {:.1}</text>"##,
                rx + 6.0,
                ry + 30.0,
                ca_core::format_z_2dp(r.activation_z),
                r.neg_log10_p
            );
        }
    }

    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use ca_core::{builtin_cohorts, seed_kb::build_seed_knowledge_graph, simulate_map_in_kg};

    #[test]
    fn pathway_svg_clips_edges_outside_node_boxes_and_includes_activity_colorbar() {
        let kg = build_seed_knowledge_graph();
        let ds = &builtin_cohorts()[0];
        let active = ds.filter_active(1.0, 0.05, &kg);
        for pw in &kg.pathways {
            let map_res = simulate_map_in_kg(pw, &active, &kg, &BTreeMap::new(), 0.82, 16);
            let svg = render_pathway_map_svg(pw, &map_res, &kg);

            assert!(svg.contains("<linearGradient id=\"map-activity-grad\""));
            assert!(svg.contains("id=\"activity-colorbar-legend\""));
            assert!(svg.contains("MAP Activity:"));

            // Collect all 74x30 node centers `(cx, cy)` from `<rect ... width="74" height="30" ...>`
            let mut node_centers = Vec::new();
            for line in svg.lines() {
                if line.contains(r#"width="74" height="30""#) {
                    let x_str = line.split("x=\"").nth(1).unwrap().split('"').next().unwrap();
                    let y_str = line.split("y=\"").nth(1).unwrap().split('"').next().unwrap();
                    let cx = x_str.parse::<f32>().unwrap() + 37.0;
                    let cy = y_str.parse::<f32>().unwrap() + 15.0;
                    node_centers.push((cx, cy));
                }
            }
            assert!(!node_centers.is_empty());

            // Verify every causal `<path d="M sx sy Q mx my ex ey"` endpoint lies strictly outside
            // all 74x30 node rect interiors `[cx±37, cy±15]`.
            let mut checked_paths = 0usize;
            for line in svg.lines() {
                if line.starts_with("<path d=\"M ") {
                    let d = line.split("d=\"").nth(1).unwrap().split('"').next().unwrap();
                    let toks: Vec<&str> = d.split_whitespace().collect();
                    // M sx sy Q mx my ex ey
                    assert_eq!(toks.len(), 8);
                    let sx: f32 = toks[1].parse().unwrap();
                    let sy: f32 = toks[2].parse().unwrap();
                    let ex: f32 = toks[6].parse().unwrap();
                    let ey: f32 = toks[7].parse().unwrap();
                    checked_paths += 1;

                    for &(cx, cy) in &node_centers {
                        let s_inside = (sx - cx).abs() < 37.0 && (sy - cy).abs() < 15.0;
                        let e_inside = (ex - cx).abs() < 37.0 && (ey - cy).abs() < 15.0;
                        assert!(
                            !s_inside,
                            "Pathway {}: edge start ({sx}, {sy}) penetrated node box at ({cx}, {cy})",
                            pw.id
                        );
                        assert!(
                            !e_inside,
                            "Pathway {}: edge end ({ex}, {ey}) penetrated node box at ({cx}, {cy})",
                            pw.id
                        );
                    }
                }
            }
            assert!(checked_paths > 0);
        }
    }
}
