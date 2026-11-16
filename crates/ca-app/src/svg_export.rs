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

/// Renders a Canonical Pathway subcellular topology + Molecule Activity Predictor (MAP) state to vector SVG.
pub fn render_pathway_map_svg(
    pathway: &PathwayDef,
    map_res: &MapSimulationResult,
    kg: &KnowledgeGraph,
) -> String {
    let wl = compute_pathway_world_layout_cached(pathway, kg);
    let width = wl.world_w.max(1100.0);
    let height = wl.world_h.max(700.0);
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
    let scale_y = (height - 40.0) / wl.world_h.max(1.0);
    let pm_y = 40.0 + wl.bands.plasma_membrane_y * scale_y;
    let org_y = 40.0 + wl.bands.organelle_divider_y * scale_y;
    let nm_y = 40.0 + wl.bands.nuclear_envelope_y * scale_y;

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
            let mx = (x1 + x2) * 0.5 + (y2 - y1) * 0.08;
            let my = (y1 + y2) * 0.5 - (x2 - x1) * 0.08;
            let _ = writeln!(
                &mut svg,
                r##"<path d="M {x1:.1} {y1:.1} Q {mx:.1} {my:.1} {x2:.1} {y2:.1}" fill="none" stroke="{stroke}" stroke-width="2.0" opacity="0.82"{marker_attr}/>"##
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
                r##"<rect x="{:.1}" y="{:.1}" width="74" height="30" rx="6" fill="{fill}" stroke="{border}" stroke-width="{bw}"/>"##,
                cx - 37.0,
                cy - 15.0
            );
            let _ = writeln!(
                &mut svg,
                r##"<text x="{cx:.1}" y="{:.1}" text-anchor="middle" font-size="11" font-weight="600" fill="#0F172A">{}</text>"##,
                cy + 4.0,
                xml_escape(nid)
            );
        }
    }

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
                r##"<text x="{:.1}" y="{:.1}" font-size="10" fill="#1E293B">z = {:+.2} | -log10(p) = {:.1}</text>"##,
                rx + 6.0,
                ry + 30.0,
                r.activation_z,
                r.neg_log10_p
            );
        }
    }

    svg.push_str("</svg>\n");
    svg
}
