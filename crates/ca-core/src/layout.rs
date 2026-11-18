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

use crate::graph::{
    EdgeClass, KnowledgeGraph, PathwayDef, PathwayReactionEdge, SubcellularLocation,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::{Mutex, OnceLock};

/// Default anisotropic inverse-square repulsion radii (`rx`, `ry`) in virtual world pixels.
pub const REPULSION_RX: f32 = 92.0;
pub const REPULSION_RY: f32 = 48.0;

/// Ideal Hooke spring rest length in virtual world pixels.
pub const HOOKE_IDEAL_LEN: f32 = 155.0;

/// Hard minimum bounding-box center-to-center separation (`min_dx`, `min_dy`) guaranteeing zero node overlap.
pub const MIN_NODE_DX: f32 = 88.0;
pub const MIN_NODE_DY: f32 = 46.0;

/// Number of deterministic physics iterations per layout solve.
pub const LAYOUT_ITERATIONS: usize = 120;

/// Normalized Y centers for the 5 subcellular compartment bands + optional phenotype band.
#[derive(Debug, Clone, Copy)]
pub struct CompartmentBands {
    pub plasma_membrane_y: f32,
    pub organelle_divider_y: f32,
    pub nuclear_envelope_y: f32,
    pub dna_helix_y: f32,
}

impl CompartmentBands {
    pub fn for_world(world_h: f32) -> Self {
        Self {
            plasma_membrane_y: 0.21 * world_h,
            organelle_divider_y: 0.59 * world_h,
            nuclear_envelope_y: 0.73 * world_h,
            dna_helix_y: 0.895 * world_h,
        }
    }
}

/// Complete deterministic world-space pathway layout (virtual canvas, non-overlapping node coords, induced edges).
#[derive(Debug, Clone)]
pub struct PathwayWorldLayout {
    pub world_w: f32,
    pub world_h: f32,
    pub positions: BTreeMap<String, (f32, f32)>,
    pub edges: Vec<PathwayReactionEdge>,
    pub bands: CompartmentBands,
}

/// Complete deterministic world-space tiered causal layout (for URA Mechanistic & Regulator Effects).
#[derive(Debug, Clone)]
pub struct TieredWorldLayout {
    pub world_w: f32,
    pub world_h: f32,
    pub positions: BTreeMap<String, (f32, f32)>,
}

/// Sizes a dynamic virtual world canvas proportional to node count `N`:
/// `world_w = max(960, sqrt(N) * 210)`, `world_h = max(760, sqrt(N) * 165)`.
pub fn compute_world_canvas_size(node_count: usize) -> (f32, f32) {
    let n_sqrt = (node_count.max(1) as f32).sqrt();
    let world_w = (n_sqrt * 210.0).max(960.0);
    let world_h = (n_sqrt * 165.0).max(760.0);
    (world_w, world_h)
}

/// Induces matching `KnowledgeGraph` causal edges among pathway members so gene-set pathways gain real topology.
pub fn induce_pathway_topology_edges(
    pathway: &PathwayDef,
    kg: &KnowledgeGraph,
) -> Vec<PathwayReactionEdge> {
    let mut member_set: BTreeSet<String> = pathway.members.iter().cloned().collect();
    for e in &pathway.topology_edges {
        member_set.insert(e.from.clone());
        member_set.insert(e.to.clone());
    }

    let mut edges = pathway.topology_edges.clone();
    let mut seen: BTreeSet<(String, String)> = edges
        .iter()
        .map(|e| (e.from.clone(), e.to.clone()))
        .collect();

    for kg_edge in &kg.edges {
        if kg_edge.edge_class == EdgeClass::Phenotype {
            continue;
        }
        if kg_edge.source != kg_edge.target
            && member_set.contains(&kg_edge.source)
            && member_set.contains(&kg_edge.target)
            && seen.insert((kg_edge.source.clone(), kg_edge.target.clone()))
        {
            edges.push(PathwayReactionEdge {
                from: kg_edge.source.clone(),
                to: kg_edge.target.clone(),
                sign: kg_edge.sign,
                interaction_label: kg_edge.mechanism.clone(),
            });
        }
    }

    edges
}

/// Maps `SubcellularLocation` onto the 5 molecular compartment band indices (`0..=4`), plus `5` for phenotypes.
pub fn subcellular_band_5(loc: SubcellularLocation) -> usize {
    match loc {
        SubcellularLocation::Extracellular => 0,
        SubcellularLocation::PlasmaMembrane => 1,
        SubcellularLocation::Cytoplasm => 2,
        SubcellularLocation::Organelle => 3,
        SubcellularLocation::Nucleus => 4,
        SubcellularLocation::PhenotypeLayer => 5,
    }
}

fn band_corridor_y(band: usize, world_h: f32) -> (f32, f32, f32) {
    // Returns (target_y, min_corridor_y, max_corridor_y)
    match band {
        0 => (0.095 * world_h, 42.0, 0.165 * world_h),
        1 => (0.210 * world_h, 0.180 * world_h, 0.240 * world_h),
        2 => (0.425 * world_h, 0.275 * world_h, 0.565 * world_h),
        3 => (0.650 * world_h, 0.605 * world_h, 0.700 * world_h),
        4 => (0.835 * world_h, 0.755 * world_h, world_h - 46.0),
        _ => (0.935 * world_h, 0.890 * world_h, world_h - 36.0),
    }
}

/// Runs 120 deterministic iterations of anisotropic inverse-square repulsion,
/// Hooke spring edge attraction, compartment Y-band restoring forces, and hard minimum
/// bounding-box overlap projection (`min_dx`, `min_dy`).
fn solve_compartment_force_layout(
    node_ids: &[String],
    node_bands: &[usize],
    node_target_y: &[f32],
    node_corridors: &[(f32, f32)],
    edge_pairs: &[(usize, usize)],
    initial_positions: &[(f32, f32)],
    world_w: f32,
    world_h: f32,
    min_dx: f32,
    min_dy: f32,
    band_stiffness: f32,
    run_force_iterations: bool,
) -> Vec<(f32, f32)> {
    let n = node_ids.len();
    if n == 0 {
        return Vec::new();
    }

    let margin_x = (min_dx * 0.58).max(52.0);
    let margin_y = (min_dy * 0.65).max(36.0);
    let max_x = (world_w - margin_x).max(margin_x + 1.0);
    let max_y = (world_h - margin_y).max(margin_y + 1.0);

    let mut pos: Vec<(f32, f32)> = initial_positions
        .iter()
        .map(|&(x, y)| (x.clamp(margin_x, max_x), y.clamp(margin_y, max_y)))
        .collect();

    if run_force_iterations && n > 1 {
        let mut disp_x = vec![0.0f32; n];
        let mut disp_y = vec![0.0f32; n];
        let center_x = world_w * 0.5;
        let rx = REPULSION_RX.max(min_dx);
        let ry = REPULSION_RY.max(min_dy);
        let inv_rx2 = 1.0 / (rx * rx);
        let inv_ry2 = 1.0 / (ry * ry);
        let cutoff_r2 = 16.0f32; // ~(4x repulsion radius)^2

        for iter in 0..LAYOUT_ITERATIONS {
            disp_x.fill(0.0);
            disp_y.fill(0.0);

            let progress = iter as f32 / LAYOUT_ITERATIONS as f32;
            let max_step = (28.0 * (1.0 - 0.72 * progress)).max(5.0);

            // 1. Anisotropic inverse-square repulsion across pairs
            for i in 0..n {
                let (xi, yi) = pos[i];
                for j in (i + 1)..n {
                    let (xj, yj) = pos[j];
                    let mut dx = xi - xj;
                    let mut dy = yi - yj;
                    if dx.abs() < 0.05 && dy.abs() < 0.05 {
                        // Deterministic tie-breaker when two nodes start nearly coincident
                        let phase = ((i * 7 + j * 13) % 16) as f32 * 0.392_699_1;
                        dx = 0.5 * phase.cos();
                        dy = 0.5 * phase.sin();
                    }
                    let nx = dx / rx;
                    let ny = dy / ry;
                    let r2 = (nx * nx + ny * ny).max(0.04);
                    if r2 > cutoff_r2 {
                        continue;
                    }
                    // Smooth falloff near cutoff so distant nodes stay cohesive
                    let cutoff_fade = (1.0 - r2 / cutoff_r2).max(0.0);
                    let inv_r3 = 1.0 / (r2 * r2.sqrt());
                    let mag = 520.0 * inv_r3 * cutoff_fade;
                    let fx = mag * dx * inv_rx2;
                    let fy = mag * dy * inv_ry2;
                    disp_x[i] += fx;
                    disp_y[i] += fy;
                    disp_x[j] -= fx;
                    disp_y[j] -= fy;
                }
            }

            // 2. Hooke spring attraction along causal edges (L_ideal = 155 px),
            //    degree-normalized and gentle across different compartment tiers so hubs do not collapse child rows.
            let mut deg = vec![1.0f32; n];
            for &(u, v) in edge_pairs {
                if u < n && v < n && u != v {
                    deg[u] += 1.0;
                    deg[v] += 1.0;
                }
            }
            for &(u, v) in edge_pairs {
                if u >= n || v >= n || u == v {
                    continue;
                }
                let dx = pos[v].0 - pos[u].0;
                let dy = pos[v].1 - pos[u].1;
                let dist = (dx * dx + dy * dy).sqrt().max(1.0);
                let cross_band = node_bands[u] != node_bands[v];
                let ideal = if cross_band {
                    (node_target_y[v] - node_target_y[u]).abs().max(HOOKE_IDEAL_LEN)
                } else {
                    HOOKE_IDEAL_LEN
                };
                let stretch = dist - ideal;
                let hub_damp = 1.0 / deg[u].max(deg[v]).sqrt();
                let k_spring = 0.085f32 * hub_damp;
                let fx = k_spring * stretch * (dx / dist);
                let fy = k_spring * stretch * (dy / dist) * if cross_band { 0.15 } else { 0.50 };
                disp_x[u] += fx;
                disp_y[u] += fy;
                disp_x[v] -= fx;
                disp_y[v] -= fy;
            }

            // 3. 5-band subcellular Y restoring force + gentle horizontal centering
            for i in 0..n {
                let (ymin, ymax) = node_corridors[i];
                let y = pos[i].1;
                let target_y = node_target_y[i];
                let mut fy = (target_y - y) * band_stiffness;
                if y < ymin {
                    fy += (ymin - y) * (band_stiffness * 3.2);
                } else if y > ymax {
                    fy += (ymax - y) * (band_stiffness * 3.2);
                }
                disp_y[i] += fy;
                disp_x[i] += (center_x - pos[i].0) * 0.035;
            }

            // 4. Integrate capped step
            for i in 0..n {
                let dx = disp_x[i];
                let dy = disp_y[i];
                let len = (dx * dx + dy * dy).sqrt();
                let (sx, sy) = if len > max_step && len > 1e-6 {
                    let s = max_step / len;
                    (dx * s, dy * s)
                } else {
                    (dx, dy)
                };
                pos[i].0 = (pos[i].0 + sx).clamp(margin_x, max_x);
                pos[i].1 = (pos[i].1 + sy).clamp(margin_y, max_y);
            }

            // Light mid-loop overlap relaxation every 15 iterations to keep stratification clean
            if iter % 15 == 14 {
                project_hard_non_overlap(
                    &mut pos,
                    node_bands,
                    node_target_y,
                    world_w,
                    min_dx,
                    min_dy,
                    margin_x,
                    max_x,
                    margin_y,
                    max_y,
                    8,
                );
            }
        }
    }

    // Final hard minimum bounding-box overlap projection guaranteeing zero overlaps
    project_hard_non_overlap(
        &mut pos,
        node_bands,
        node_target_y,
        world_w,
        min_dx,
        min_dy,
        margin_x,
        max_x,
        margin_y,
        max_y,
        320,
    );

    pos
}

/// Iterative Gauss-Seidel AABB separation ensuring every pair `(i, j)` satisfies
/// `|x_i - x_j| >= min_dx || |y_i - y_j| >= min_dy`.
#[allow(clippy::too_many_arguments)]
fn project_hard_non_overlap(
    pos: &mut [(f32, f32)],
    node_bands: &[usize],
    node_target_y: &[f32],
    world_w: f32,
    min_dx: f32,
    min_dy: f32,
    margin_x: f32,
    max_x: f32,
    margin_y: f32,
    max_y: f32,
    max_sweeps: usize,
) {
    let n = pos.len();
    if n < 2 {
        return;
    }
    let eps = 0.08f32;
    let req_dx = min_dx + eps;
    let req_dy = min_dy + eps;

    for sweep in 0..max_sweeps {
        let mut any_overlap = false;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos[i].0 - pos[j].0;
                let dy = pos[i].1 - pos[j].1;
                let adx = dx.abs();
                let ady = dy.abs();
                if adx < req_dx && ady < req_dy {
                    any_overlap = true;
                    let ox = req_dx - adx;
                    let oy = req_dy - ady;
                    let same_band = node_bands[i] == node_bands[j];
                    let norm_ox = ox / req_dx;
                    let norm_oy = oy / req_dy;

                    // Prefer horizontal separation within the same band when penetration is comparable
                    let separate_x = if same_band && sweep < max_sweeps / 2 {
                        norm_ox <= norm_oy * 1.35
                    } else {
                        norm_ox <= norm_oy
                    };

                    if separate_x {
                        let sgn = if dx > 0.0 || (dx == 0.0 && (i + j) % 2 == 0) {
                            1.0
                        } else {
                            -1.0
                        };
                        let push = ox * 0.54;
                        let xi_new = (pos[i].0 + sgn * push).clamp(margin_x, max_x);
                        let xj_new = (pos[j].0 - sgn * push).clamp(margin_x, max_x);
                        let applied_i = (xi_new - pos[i].0).abs();
                        let applied_j = (xj_new - pos[j].0).abs();
                        pos[i].0 = xi_new;
                        pos[j].0 = xj_new;
                        if applied_i + applied_j < ox {
                            let rem = ox - (applied_i + applied_j);
                            if applied_i < applied_j {
                                pos[j].0 = (pos[j].0 - sgn * rem).clamp(margin_x, max_x);
                            } else {
                                pos[i].0 = (pos[i].0 + sgn * rem).clamp(margin_x, max_x);
                            }
                        }
                    } else {
                        let sgn = if dy > 0.0 || (dy == 0.0 && node_bands[i] >= node_bands[j]) {
                            1.0
                        } else {
                            -1.0
                        };
                        let push = oy * 0.54;
                        let yi_new = (pos[i].1 + sgn * push).clamp(margin_y, max_y);
                        let yj_new = (pos[j].1 - sgn * push).clamp(margin_y, max_y);
                        let applied_i = (yi_new - pos[i].1).abs();
                        let applied_j = (yj_new - pos[j].1).abs();
                        pos[i].1 = yi_new;
                        pos[j].1 = yj_new;
                        if applied_i + applied_j < oy {
                            let rem = oy - (applied_i + applied_j);
                            if applied_i < applied_j {
                                pos[j].1 = (pos[j].1 - sgn * rem).clamp(margin_y, max_y);
                            } else {
                                pos[i].1 = (pos[i].1 + sgn * rem).clamp(margin_y, max_y);
                            }
                        }
                    }
                }
            }
        }
        if !any_overlap {
            break;
        }
    }

    // Per-band centered fallback grid repair ONLY if extreme density still leaves any overlapping pair
    if has_any_overlap(pos, min_dx, min_dy) {
        let usable_w = (max_x - margin_x).max(req_dx);
        let cols = (usable_w / (req_dx + 8.0)).floor().max(1.0) as usize;
        let mut by_band: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (idx, &b) in node_bands.iter().enumerate() {
            by_band.entry(b).or_default().push(idx);
        }
        let mut next_band_min_y = margin_y;
        for indices in by_band.values_mut() {
            indices.sort_by(|&a, &b| pos[a].0.total_cmp(&pos[b].0).then_with(|| pos[a].1.total_cmp(&pos[b].1)));
            let rows = (indices.len() + cols - 1) / cols;
            let avg_target_y = indices.iter().map(|&i| node_target_y[i]).sum::<f32>()
                / (indices.len().max(1) as f32);
            let band_y = avg_target_y.max(next_band_min_y).min(max_y - (rows.saturating_sub(1) as f32) * req_dy);
            let step_x = req_dx + 10.0;
            for (slot, &idx) in indices.iter().enumerate() {
                let row = slot / cols;
                let col = slot % cols;
                let row_len = if row + 1 == rows && indices.len() % cols != 0 {
                    indices.len() % cols
                } else {
                    cols.min(indices.len())
                };
                let row_w = (row_len.saturating_sub(1) as f32) * step_x;
                let x0 = (world_w * 0.5 - row_w * 0.5).clamp(margin_x, max_x);
                pos[idx].0 = (x0 + (col as f32) * step_x).min(max_x);
                pos[idx].1 = (band_y + (row as f32) * (req_dy + 4.0)).min(max_y);
            }
            next_band_min_y = band_y + (rows as f32) * (req_dy + 8.0) + 20.0;
        }
    }
}

fn has_any_overlap(pos: &[(f32, f32)], min_dx: f32, min_dy: f32) -> bool {
    for i in 0..pos.len() {
        for j in (i + 1)..pos.len() {
            if (pos[i].0 - pos[j].0).abs() < min_dx - 1e-3
                && (pos[i].1 - pos[j].1).abs() < min_dy - 1e-3
            {
                return true;
            }
        }
    }
    false
}

/// Computes deterministic multi-row initial coordinates inside each subcellular compartment band.
fn initial_compartment_grid(
    node_ids: &[String],
    node_bands: &[usize],
    edges: &[PathwayReactionEdge],
    world_w: f32,
    world_h: f32,
) -> (Vec<(f32, f32)>, Vec<f32>, Vec<(f32, f32)>) {
    let n = node_ids.len();
    let id_to_idx: HashMap<&str, usize> = node_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.as_str(), i))
        .collect();

    // Compute normalized topological depth within Cytoplasm (band 2) for vertical cascade seeding
    let mut cyto_depth = vec![0usize; n];
    for _ in 0..6 {
        for e in edges {
            if let (Some(&u), Some(&v)) = (id_to_idx.get(e.from.as_str()), id_to_idx.get(e.to.as_str()))
            {
                if node_bands[u] == 2 && node_bands[v] == 2 && cyto_depth[v] <= cyto_depth[u] {
                    cyto_depth[v] = (cyto_depth[u] + 1).min(4);
                }
            }
        }
    }
    let max_cyto_depth = cyto_depth.iter().copied().max().unwrap_or(0).max(1);

    let mut by_band: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, &b) in node_bands.iter().enumerate() {
        by_band.entry(b).or_default().push(idx);
    }

    let mut init_pos = vec![(world_w * 0.5, world_h * 0.5); n];
    let mut target_y = vec![world_h * 0.5; n];
    let mut corridors = vec![(40.0, world_h - 40.0); n];

    let margin_x = 74.0f32;
    let usable_w = (world_w - 2.0 * margin_x).max(120.0);

    for (band, mut indices) in by_band {
        if band == 2 {
            indices.sort_by(|&a, &b| {
                cyto_depth[a]
                    .cmp(&cyto_depth[b])
                    .then_with(|| node_ids[a].cmp(&node_ids[b]))
            });
        } else {
            indices.sort_by(|&a, &b| node_ids[a].cmp(&node_ids[b]));
        }

        let (band_ty, cmin, cmax) = band_corridor_y(band, world_h);
        let max_cols = ((usable_w / 108.0).floor() as usize).max(1);
        let m = indices.len();
        let rows = (m + max_cols - 1) / max_cols;
        let row_step = if rows <= 1 {
            0.0
        } else {
            ((cmax - cmin) / (rows as f32)).clamp(50.0, 96.0)
        };
        let total_span_y = (rows.saturating_sub(1) as f32) * row_step;
        let start_y = (band_ty - total_span_y * 0.5).clamp(cmin, (cmax - total_span_y).max(cmin));

        for (k, &idx) in indices.iter().enumerate() {
            let r = k / max_cols;
            let c = k % max_cols;
            let row_items = if r + 1 == rows && m % max_cols != 0 {
                m % max_cols
            } else {
                max_cols.min(m)
            };
            let stagger = if r % 2 == 1 { 22.0 } else { 0.0 };
            let x_span = (row_items.saturating_sub(1) as f32) * ((usable_w / row_items.max(1) as f32).min(165.0));
            let x0 = world_w * 0.5 - x_span * 0.5 + stagger;
            let step_x = if row_items <= 1 {
                0.0
            } else {
                x_span / (row_items - 1) as f32
            };
            let px = (x0 + (c as f32) * step_x).clamp(margin_x, world_w - margin_x);
            let py = if band == 2 && max_cyto_depth > 0 && rows == 1 {
                let frac = cyto_depth[idx] as f32 / max_cyto_depth as f32;
                cmin + (0.15 + 0.70 * frac) * (cmax - cmin)
            } else {
                start_y + (r as f32) * row_step
            };
            init_pos[idx] = (px, py);
            target_y[idx] = if band == 2 { py } else { band_ty };
            corridors[idx] = (cmin, cmax);
        }
    }

    (init_pos, target_y, corridors)
}

/// Computes the full world-space layout for a `PathwayDef`, either uncached or forced re-layout.
pub fn compute_pathway_world_layout_uncached(
    pathway: &PathwayDef,
    kg: &KnowledgeGraph,
    force_relayout: bool,
) -> PathwayWorldLayout {
    let edges = induce_pathway_topology_edges(pathway, kg);
    let mut all_ids_set: BTreeSet<String> = pathway.members.iter().cloned().collect();
    for e in &edges {
        all_ids_set.insert(e.from.clone());
        all_ids_set.insert(e.to.clone());
    }
    let node_ids: Vec<String> = all_ids_set.into_iter().collect();
    let n = node_ids.len();
    let (world_w, world_h) = compute_world_canvas_size(n);
    let bands = CompartmentBands::for_world(world_h);

    let node_bands: Vec<usize> = node_ids
        .iter()
        .map(|nid| {
            kg.nodes
                .get(nid)
                .map(|node| subcellular_band_5(node.location))
                .unwrap_or(2)
        })
        .collect();

    let id_to_idx: HashMap<&str, usize> = node_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.as_str(), i))
        .collect();

    let edge_pairs: Vec<(usize, usize)> = edges
        .iter()
        .filter_map(|e| {
            let u = *id_to_idx.get(e.from.as_str())?;
            let v = *id_to_idx.get(e.to.as_str())?;
            Some((u, v))
        })
        .collect();

    let has_full_curated_coords = !force_relayout
        && !pathway.node_coords.is_empty()
        && node_ids.iter().all(|id| pathway.node_coords.contains_key(id));

    let (mut init_pos, target_y, corridors) =
        initial_compartment_grid(&node_ids, &node_bands, &edges, world_w, world_h);

    if has_full_curated_coords {
        for (i, nid) in node_ids.iter().enumerate() {
            if let Some(&(nx, ny)) = pathway.node_coords.get(nid) {
                init_pos[i] = (
                    (0.05 + nx.clamp(0.0, 1.0) * 0.90) * world_w,
                    (0.04 + ny.clamp(0.0, 1.0) * 0.92) * world_h,
                );
            }
        }
    }

    // Preserve hand-curated GPML / seed coordinates with strict non-overlap projection;
    // run full 120-iteration Inverse-Square Repulsion + Hooke Spring + 5-band physics otherwise.
    let solved = solve_compartment_force_layout(
        &node_ids,
        &node_bands,
        &target_y,
        &corridors,
        &edge_pairs,
        &init_pos,
        world_w,
        world_h,
        MIN_NODE_DX,
        MIN_NODE_DY,
        0.16,
        !has_full_curated_coords,
    );

    let positions: BTreeMap<String, (f32, f32)> = node_ids.into_iter().zip(solved).collect();

    PathwayWorldLayout {
        world_w,
        world_h,
        positions,
        edges,
        bands,
    }
}

#[derive(Debug, Clone)]
struct CachedLayoutEntry {
    kg_edge_count: usize,
    node_coords_len: usize,
    layout: PathwayWorldLayout,
}

fn layout_cache() -> &'static Mutex<HashMap<(String, usize), CachedLayoutEntry>> {
    static CACHE: OnceLock<Mutex<HashMap<(String, usize), CachedLayoutEntry>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Returns a deterministic cached `PathwayWorldLayout` keyed by `(pathway_id, N)`
/// (automatically refreshing if `kg.edges` or curated `node_coords` change).
pub fn compute_pathway_world_layout_cached(
    pathway: &PathwayDef,
    kg: &KnowledgeGraph,
) -> PathwayWorldLayout {
    let mut all_ids: BTreeSet<&str> = pathway.members.iter().map(String::as_str).collect();
    for e in &pathway.topology_edges {
        all_ids.insert(e.from.as_str());
        all_ids.insert(e.to.as_str());
    }
    let key = (pathway.id.clone(), all_ids.len());
    let kg_edges_len = kg.edges.len();
    let coords_len = pathway.node_coords.len();

    if let Ok(guard) = layout_cache().lock() {
        if let Some(entry) = guard.get(&key) {
            if entry.kg_edge_count == kg_edges_len && entry.node_coords_len == coords_len {
                return entry.layout.clone();
            }
        }
    }

    let computed = compute_pathway_world_layout_uncached(pathway, kg, false);
    if let Ok(mut guard) = layout_cache().lock() {
        if guard.len() > 256 {
            guard.clear();
        }
        guard.insert(
            key,
            CachedLayoutEntry {
                kg_edge_count: kg_edges_len,
                node_coords_len: coords_len,
                layout: computed.clone(),
            },
        );
    }
    computed
}

/// Clears or refreshes a single pathway's cached layout with a forced 120-iteration physics re-layout.
pub fn recompute_and_cache_pathway_layout(
    pathway: &PathwayDef,
    kg: &KnowledgeGraph,
) -> PathwayWorldLayout {
    let mut all_ids: BTreeSet<&str> = pathway.members.iter().map(String::as_str).collect();
    for e in &pathway.topology_edges {
        all_ids.insert(e.from.as_str());
        all_ids.insert(e.to.as_str());
    }
    let key = (pathway.id.clone(), all_ids.len());
    let computed = compute_pathway_world_layout_uncached(pathway, kg, true);
    if let Ok(mut guard) = layout_cache().lock() {
        guard.insert(
            key,
            CachedLayoutEntry {
                kg_edge_count: kg.edges.len(),
                node_coords_len: pathway.node_coords.len(),
                layout: computed.clone(),
            },
        );
    }
    computed
}

/// Normalized `[0.0, 1.0]^2` wrapper used by `PathwayDef::resolve_layout_coords` and SVG exporters.
pub fn resolve_normalized_pathway_coords(
    pathway: &PathwayDef,
    kg: &KnowledgeGraph,
) -> BTreeMap<String, (f32, f32)> {
    let wl = compute_pathway_world_layout_cached(pathway, kg);
    let w = wl.world_w.max(1.0);
    let h = wl.world_h.max(1.0);
    wl.positions
        .into_iter()
        .map(|(k, (x, y))| (k, ((x / w).clamp(0.02, 0.98), (y / h).clamp(0.02, 0.98))))
        .collect()
}

/// Input specification for a hierarchical/tiered causal network node (URA Mechanistic or Regulator Effects).
#[derive(Debug, Clone)]
pub struct TieredNodeSpec {
    pub id: String,
    pub tier: usize,
}

/// Computes a deterministic, zero-overlap world-space layout for tiered causal graphs
/// (used by URA Mechanistic Networks and 3-Tier Regulator Effects Explorer).
pub fn compute_tiered_causal_layout(
    cache_key: &str,
    nodes: &[TieredNodeSpec],
    edges: &[(String, String)],
    min_dx: f32,
    min_dy: f32,
) -> TieredWorldLayout {
    let mut unique_map: BTreeMap<String, usize> = BTreeMap::new();
    for node in nodes {
        unique_map.entry(node.id.clone()).or_insert(node.tier);
    }
    let node_ids: Vec<String> = unique_map.keys().cloned().collect();
    let n = node_ids.len();
    if n == 0 {
        return TieredWorldLayout {
            world_w: 960.0,
            world_h: 760.0,
            positions: BTreeMap::new(),
        };
    }

    let key = (format!("__TIERED__{cache_key}"), n);
    if let Ok(guard) = layout_cache().lock() {
        if let Some(entry) = guard.get(&key) {
            if entry.kg_edge_count == edges.len() {
                return TieredWorldLayout {
                    world_w: entry.layout.world_w,
                    world_h: entry.layout.world_h,
                    positions: entry.layout.positions.clone(),
                };
            }
        }
    }

    let (base_w, base_h) = compute_world_canvas_size(n);
    // Ensure wide nodes (such as BiologicalProcess octagons in Regulator Effects) always fit cleanly
    let tier_count = unique_map.values().copied().max().unwrap_or(0) + 1;
    let mut max_in_tier = 1usize;
    let mut tier_buckets: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (idx, nid) in node_ids.iter().enumerate() {
        let t = unique_map[nid];
        tier_buckets.entry(t).or_default().push(idx);
    }
    for b in tier_buckets.values() {
        max_in_tier = max_in_tier.max(b.len());
    }
    let min_width_for_tier = ((max_in_tier.min(10) as f32) * (min_dx + 18.0) + 140.0).min(3600.0);
    let world_w = base_w.max(min_width_for_tier);
    let world_h = base_h.max((tier_count as f32) * (min_dy * 3.2)).max(760.0);

    let node_bands: Vec<usize> = node_ids.iter().map(|id| unique_map[id]).collect();
    let max_tier = tier_count.saturating_sub(1).max(1) as f32;

    let mut init_pos = vec![(world_w * 0.5, world_h * 0.5); n];
    let mut target_y = vec![world_h * 0.5; n];
    let mut corridors = vec![(48.0, world_h - 48.0); n];

    let margin_x = (min_dx * 0.58).max(68.0);
    let usable_w = (world_w - 2.0 * margin_x).max(120.0);

    for (tier, indices) in tier_buckets {
        let frac = if tier_count <= 1 {
            0.5
        } else {
            0.12 + (tier as f32 / max_tier) * 0.76
        };
        let ty = frac * world_h;
        let band_half = (world_h / (tier_count.max(2) as f32)) * 0.34;
        let cmin = (ty - band_half).max(48.0);
        let cmax = (ty + band_half).min(world_h - 48.0).max(cmin + 10.0);

        let max_cols = ((usable_w / (min_dx + 22.0)).floor() as usize).max(1);
        let m = indices.len();
        let rows = (m + max_cols - 1) / max_cols;
        let row_step = if rows <= 1 {
            0.0
        } else {
            ((cmax - cmin) / (rows as f32)).clamp(min_dy + 6.0, 96.0)
        };
        let start_y = (ty - (rows.saturating_sub(1) as f32) * row_step * 0.5).clamp(cmin, cmax);

        for (k, &idx) in indices.iter().enumerate() {
            let r = k / max_cols;
            let c = k % max_cols;
            let row_items = if r + 1 == rows && m % max_cols != 0 {
                m % max_cols
            } else {
                max_cols.min(m)
            };
            let span_x = (row_items.saturating_sub(1) as f32) * ((usable_w / row_items.max(1) as f32).min(min_dx + 75.0));
            let x0 = world_w * 0.5 - span_x * 0.5;
            let step_x = if row_items <= 1 {
                0.0
            } else {
                span_x / (row_items - 1) as f32
            };
            init_pos[idx] = (
                (x0 + (c as f32) * step_x).clamp(margin_x, world_w - margin_x),
                (start_y + (r as f32) * row_step).clamp(48.0, world_h - 48.0),
            );
            target_y[idx] = ty;
            corridors[idx] = (cmin, cmax);
        }
    }

    let id_to_idx: HashMap<&str, usize> = node_ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.as_str(), i))
        .collect();
    let edge_pairs: Vec<(usize, usize)> = edges
        .iter()
        .filter_map(|(a, b)| {
            Some((*id_to_idx.get(a.as_str())?, *id_to_idx.get(b.as_str())?))
        })
        .collect();

    let solved = solve_compartment_force_layout(
        &node_ids,
        &node_bands,
        &target_y,
        &corridors,
        &edge_pairs,
        &init_pos,
        world_w,
        world_h,
        min_dx.max(MIN_NODE_DX),
        min_dy.max(MIN_NODE_DY),
        0.24,
        true,
    );

    let positions: BTreeMap<String, (f32, f32)> = node_ids.into_iter().zip(solved).collect();

    if let Ok(mut guard) = layout_cache().lock() {
        if guard.len() > 256 {
            guard.clear();
        }
        guard.insert(
            key,
            CachedLayoutEntry {
                kg_edge_count: edges.len(),
                node_coords_len: 0,
                layout: PathwayWorldLayout {
                    world_w,
                    world_h,
                    positions: positions.clone(),
                    edges: Vec::new(),
                    bands: CompartmentBands::for_world(world_h),
                },
            },
        );
    }

    TieredWorldLayout {
        world_w,
        world_h,
        positions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{CausalEdge, EdgeSign, EntityKind, EntityNode};

    #[test]
    fn world_canvas_scales_with_sqrt_n_exactly() {
        assert_eq!(compute_world_canvas_size(1), (960.0, 760.0));
        assert_eq!(compute_world_canvas_size(16), (960.0, 760.0));
        let (w144, h144) = compute_world_canvas_size(144);
        assert!((w144 - 2520.0).abs() < 1e-3, "got {w144}");
        assert!((h144 - 1980.0).abs() < 1e-3, "got {h144}");
    }

    #[test]
    fn compartment_force_layout_144_nodes_has_zero_overlap_and_induces_edges() {
        let mut kg = KnowledgeGraph::new();
        let mut members = Vec::with_capacity(144);

        for i in 0..144usize {
            let sym = format!("G{:03}", i);
            let loc = if i < 10 {
                SubcellularLocation::Extracellular
            } else if i < 22 {
                SubcellularLocation::PlasmaMembrane
            } else if i < 124 {
                SubcellularLocation::Cytoplasm
            } else if i < 132 {
                SubcellularLocation::Organelle
            } else {
                SubcellularLocation::Nucleus
            };
            kg.upsert_node(EntityNode {
                id: sym.clone(),
                display_name: sym.clone(),
                kind: EntityKind::GeneOther,
                location: loc,
                family_category: "Test".to_string(),
                external_id: None,
                provenance_source: "Test".to_string(),
            });
            members.push(sym);
        }

        // Add causal KG edges without putting them into PathwayDef.topology_edges
        for i in 0..80usize {
            kg.add_edge(CausalEdge {
                source: format!("G{:03}", i),
                target: format!("G{:03}", i + 12),
                edge_class: EdgeClass::Activation,
                sign: EdgeSign::Activating,
                weight: 0.9,
                mechanism: "induced signal".to_string(),
                source_db: "SIGNOR".to_string(),
                pmids: vec![],
            });
        }
        kg.rebuild_indices();

        let pw = PathwayDef {
            id: "TEST-144".to_string(),
            name: "Dense 144-Node Cytoplasmic Cascade".to_string(),
            source_db: "Test".to_string(),
            category: "Test".to_string(),
            url: String::new(),
            members,
            member_expected_signs: BTreeMap::new(),
            topology_edges: Vec::new(),
            node_coords: BTreeMap::new(),
            phospho_nodes: Vec::new(),
            cell_context_label: String::new(),
        };

        let layout = compute_pathway_world_layout_cached(&pw, &kg);
        assert!((layout.world_w - 2520.0).abs() < 1e-3);
        assert!((layout.world_h - 1980.0).abs() < 1e-3);
        assert_eq!(layout.positions.len(), 144);
        assert!(
            layout.edges.len() >= 80,
            "Expected induced KG edges, got {}",
            layout.edges.len()
        );

        let coords: Vec<(String, (f32, f32))> = layout.positions.clone().into_iter().collect();
        for i in 0..coords.len() {
            for j in (i + 1)..coords.len() {
                let dx = (coords[i].1 .0 - coords[j].1 .0).abs();
                let dy = (coords[i].1 .1 - coords[j].1 .1).abs();
                assert!(
                    dx >= MIN_NODE_DX - 1e-2 || dy >= MIN_NODE_DY - 1e-2,
                    "Overlap detected between {} ({:?}) and {} ({:?}): dx={dx:.2}, dy={dy:.2}",
                    coords[i].0,
                    coords[i].1,
                    coords[j].0,
                    coords[j].1
                );
            }
        }

        // Verify deterministic caching returns bit-for-bit identical positions
        let layout2 = compute_pathway_world_layout_cached(&pw, &kg);
        assert_eq!(layout.positions, layout2.positions);

        // Verify subcellular band vertical stratification
        let mean_y = |lo: usize, hi: usize| -> f32 {
            let sum: f32 = (lo..hi)
                .map(|i| layout.positions[&format!("G{:03}", i)].1)
                .sum();
            sum / (hi - lo) as f32
        };
        let y_ext = mean_y(0, 10);
        let y_pm = mean_y(10, 22);
        let y_cyto = mean_y(22, 124);
        let y_org = mean_y(124, 132);
        let y_nuc = mean_y(132, 144);
        assert!(
            y_ext < y_pm && y_pm < y_cyto && y_cyto < y_org && y_org < y_nuc,
            "Subcellular bands out of order: ext={y_ext:.1}, pm={y_pm:.1}, cyto={y_cyto:.1}, org={y_org:.1}, nuc={y_nuc:.1}"
        );
    }

    #[test]
    fn tiered_causal_layout_prevents_overlap_on_dense_targets() {
        let mut nodes = vec![TieredNodeSpec {
            id: "TNF".to_string(),
            tier: 0,
        }];
        let mut edges = Vec::new();
        for i in 0..4 {
            let r = format!("REG{i}");
            edges.push(("TNF".to_string(), r.clone()));
            nodes.push(TieredNodeSpec { id: r, tier: 1 });
        }
        for i in 0..16 {
            let g = format!("GENE{i:02}");
            edges.push(("TNF".to_string(), g.clone()));
            nodes.push(TieredNodeSpec { id: g, tier: 2 });
        }

        let res = compute_tiered_causal_layout("test-ura", &nodes, &edges, MIN_NODE_DX, MIN_NODE_DY);
        let coords: Vec<(f32, f32)> = res.positions.values().copied().collect();
        for i in 0..coords.len() {
            for j in (i + 1)..coords.len() {
                let dx = (coords[i].0 - coords[j].0).abs();
                let dy = (coords[i].1 - coords[j].1).abs();
                assert!(
                    dx >= MIN_NODE_DX - 1e-2 || dy >= MIN_NODE_DY - 1e-2,
                    "Tiered layout overlap: dx={dx:.2}, dy={dy:.2}"
                );
            }
        }
    }
}
