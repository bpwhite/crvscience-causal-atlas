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

fn blend_px(buf: &mut [u8], width: u32, x: i32, y: i32, rgba: [u8; 4], alpha_scale: f32) {
    if x < 0 || y < 0 || x >= width as i32 || y >= width as i32 {
        return;
    }
    let a = ((rgba[3] as f32 / 255.0) * alpha_scale).clamp(0.0, 1.0);
    if a <= 0.001 {
        return;
    }
    let idx = ((y as u32 * width + x as u32) * 4) as usize;
    let inv = 1.0 - a;
    buf[idx] = (rgba[0] as f32 * a + buf[idx] as f32 * inv).round() as u8;
    buf[idx + 1] = (rgba[1] as f32 * a + buf[idx + 1] as f32 * inv).round() as u8;
    buf[idx + 2] = (rgba[2] as f32 * a + buf[idx + 2] as f32 * inv).round() as u8;
    let dst_a = buf[idx + 3] as f32 / 255.0;
    let out_a = (a + dst_a * inv).clamp(0.0, 1.0);
    buf[idx + 3] = (out_a * 255.0).round() as u8;
}

fn draw_disc(buf: &mut [u8], size: u32, cx: f32, cy: f32, r: f32, rgba: [u8; 4]) {
    let min_x = (cx - r - 1.5).floor() as i32;
    let max_x = (cx + r + 1.5).ceil() as i32;
    let min_y = (cy - r - 1.5).floor() as i32;
    let max_y = (cy + r + 1.5).ceil() as i32;
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let d = ((x as f32 + 0.5 - cx).hypot(y as f32 + 0.5 - cy)) - r;
            let cov = (0.5 - d).clamp(0.0, 1.0);
            if cov > 0.0 {
                blend_px(buf, size, x, y, rgba, cov);
            }
        }
    }
}

fn draw_segment(
    buf: &mut [u8],
    size: u32,
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    half_w: f32,
    rgba: [u8; 4],
) {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let len2 = dx * dx + dy * dy;
    if len2 <= 1e-6 {
        return;
    }
    let min_x = (x0.min(x1) - half_w - 1.5).floor() as i32;
    let max_x = (x0.max(x1) + half_w + 1.5).ceil() as i32;
    let min_y = (y0.min(y1) - half_w - 1.5).floor() as i32;
    let max_y = (y0.max(y1) + half_w + 1.5).ceil() as i32;
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let px = x as f32 + 0.5 - x0;
            let py = y as f32 + 0.5 - y0;
            let t = ((px * dx + py * dy) / len2).clamp(0.0, 1.0);
            let dist = (px - dx * t).hypot(py - dy * t) - half_w;
            let cov = (0.5 - dist).clamp(0.0, 1.0);
            if cov > 0.0 {
                blend_px(buf, size, x, y, rgba, cov);
            }
        }
    }
}

/// Renders the Causal-Atlas RGBA icon (`size x size`) depicting a signed hierarchical causal network
/// over a dark slate rounded badge.
pub fn render_icon_rgba(size: u32) -> Vec<u8> {
    let mut buf = vec![0u8; (size * size * 4) as usize];
    let s = size as f32;

    // Rounded badge background (#11151C) with subtle slate rim (#3A4659)
    let pad = s * 0.04;
    let corner_r = s * 0.20;
    let half = s * 0.5 - pad;
    let center = s * 0.5;
    let border_w = (s * 0.045).max(1.0);

    for y in 0..size {
        for x in 0..size {
            let px = (x as f32 + 0.5 - center).abs();
            let py = (y as f32 + 0.5 - center).abs();
            let qx = (px - (half - corner_r)).max(0.0);
            let qy = (py - (half - corner_r)).max(0.0);
            let sdf = qx.hypot(qy) - corner_r;
            let outer_cov = (0.5 - sdf).clamp(0.0, 1.0);
            if outer_cov > 0.0 {
                let is_rim = (sdf + border_w).clamp(0.0, 1.0);
                let bg = [17, 21, 28, 255];
                let rim = [58, 70, 89, 255];
                let r = (bg[0] as f32 * (1.0 - is_rim) + rim[0] as f32 * is_rim).round() as u8;
                let g = (bg[1] as f32 * (1.0 - is_rim) + rim[1] as f32 * is_rim).round() as u8;
                let b = (bg[2] as f32 * (1.0 - is_rim) + rim[2] as f32 * is_rim).round() as u8;
                blend_px(&mut buf, size, x as i32, y as i32, [r, g, b, 255], outer_cov);
            }
        }
    }

    // Subtle subcellular membrane divider lines
    let mem_col = [42, 50, 64, 210];
    draw_segment(
        &mut buf,
        size,
        s * 0.14,
        s * 0.34,
        s * 0.86,
        s * 0.34,
        (s * 0.016).max(0.6),
        mem_col,
    );
    draw_segment(
        &mut buf,
        size,
        s * 0.14,
        s * 0.67,
        s * 0.86,
        s * 0.67,
        (s * 0.016).max(0.6),
        mem_col,
    );

    // Causal network coordinates
    let root = (s * 0.50, s * 0.22);
    let mid_l = (s * 0.30, s * 0.49);
    let mid_r = (s * 0.70, s * 0.49);
    let bot_l = (s * 0.24, s * 0.77);
    let bot_c = (s * 0.50, s * 0.78);
    let bot_r = (s * 0.76, s * 0.77);

    let line_w = (s * 0.036).max(0.9);
    let orange = [230, 105, 18, 255];
    let blue = [43, 131, 186, 255];
    let teal = [27, 158, 119, 255];
    let coral = [217, 72, 48, 255];
    let gold = [250, 204, 21, 255];

    // Directed causal edges
    draw_segment(&mut buf, size, root.0, root.1, mid_l.0, mid_l.1, line_w, orange);
    draw_segment(&mut buf, size, root.0, root.1, mid_r.0, mid_r.1, line_w, blue);
    draw_segment(&mut buf, size, mid_l.0, mid_l.1, bot_l.0, bot_l.1, line_w, orange);
    draw_segment(&mut buf, size, mid_l.0, mid_l.1, bot_c.0, bot_c.1, line_w, orange);
    draw_segment(&mut buf, size, mid_r.0, mid_r.1, bot_c.0, bot_c.1, line_w, blue);
    draw_segment(&mut buf, size, mid_r.0, mid_r.1, bot_r.0, bot_r.1, line_w, teal);

    // Node glows & discs
    let nr = (s * 0.092).max(1.8);
    let root_r = (s * 0.115).max(2.2);
    let rim_r = (s * 0.025).max(0.7);

    // Gold selection halo around root & central hub
    draw_disc(&mut buf, size, root.0, root.1, root_r + rim_r * 1.4, [250, 204, 21, 160]);
    draw_disc(&mut buf, size, root.0, root.1, root_r, orange);

    // Mid tier
    draw_disc(&mut buf, size, mid_l.0, mid_l.1, nr + rim_r, [232, 237, 244, 220]);
    draw_disc(&mut buf, size, mid_l.0, mid_l.1, nr, gold);

    draw_disc(&mut buf, size, mid_r.0, mid_r.1, nr + rim_r, [232, 237, 244, 220]);
    draw_disc(&mut buf, size, mid_r.0, mid_r.1, nr, blue);

    // Bottom target tier
    draw_disc(&mut buf, size, bot_l.0, bot_l.1, nr + rim_r * 0.8, [232, 237, 244, 200]);
    draw_disc(&mut buf, size, bot_l.0, bot_l.1, nr * 0.94, coral);

    draw_disc(&mut buf, size, bot_c.0, bot_c.1, nr + rim_r * 0.8, [232, 237, 244, 200]);
    draw_disc(&mut buf, size, bot_c.0, bot_c.1, nr * 0.94, orange);

    draw_disc(&mut buf, size, bot_r.0, bot_r.1, nr + rim_r * 0.8, [232, 237, 244, 200]);
    draw_disc(&mut buf, size, bot_r.0, bot_r.1, nr * 0.94, teal);

    buf
}

/// Returns a 64x64 RGBA `egui::IconData` for `ViewportBuilder::with_icon`.
pub fn app_icon_data() -> egui::IconData {
    egui::IconData {
        rgba: render_icon_rgba(64),
        width: 64,
        height: 64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn app_icon_has_expected_rgba_dimensions_and_non_uniform_pixels() {
        let icon = app_icon_data();
        assert_eq!(icon.width, 64);
        assert_eq!(icon.height, 64);
        assert_eq!(icon.rgba.len(), 64 * 64 * 4);
        let opaque_count = icon.rgba.chunks_exact(4).filter(|px| px[3] > 200).count();
        assert!(opaque_count > 2000, "icon must contain visible opaque badge");
        let first = &icon.rgba[0..4];
        assert!(icon.rgba.chunks_exact(4).any(|px| px != first));
    }
}
